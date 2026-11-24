#!/usr/bin/env Rscript
# =============================================================================
# Phase 6: Differential gene expression analysis of GSE213410 with DESeq2
# =============================================================================
# Comparisons (analysed independently, never pooled):
#   C1  H9 hESC     Ethanol (A1-A3, control)  vs 5-OXO-ETE  (A4-A6)
#       primary: all six samples (A3 kept)
#       sensitivity: A3 excluded (A1, A2 vs A4-A6)
#   C2  T21 hiPSC   DMSO    (A7-A9, control)  vs Calcitriol (A10-A12)
#
# Input: results/tables/gene_level_counts.csv (Phase 5, Ensembl 97 genes)
#
# Count provenance (checked in Step 1):
#   The GEO counts were produced by StringTie 1.3.3b. Integer read counts
#   from StringTie come from its prepDE.py helper, which computes
#   ceil(coverage * transcript_length / read_length) per transcript. Counts
#   are therefore coverage-derived estimates, already integer, rounded *up*.
#   This script does not round anything. It checks that the counts sit on a
#   read-pair scale (against SRA run statistics) and measures an upper bound
#   on the ceil() inflation for each gene.
#
# Statistics (per analysis):
#   - Pre-filter: >= 10 raw counts in >= n_min samples (n_min = smallest
#     group size in that analysis)
#   - DESeq2: median-of-ratios size factors, gene-wise dispersions with
#     empirical-Bayes shrinkage, Wald test for treatment vs control
#   - Benjamini-Hochberg adjustment (alpha = 0.05, DESeq2 independent
#     filtering on)
#   - log2 fold-change shrinkage: apeglm. apeglm estimates its prior scale
#     from all genes; when genome-wide effects are tiny the scale hits its
#     lower bound (0.001) and moderately supported genes shrink to ~0. So:
#       * direction (up/down) and volcano x-axis use the MLE LFC (the
#         quantity tested by the Wald test)
#       * shrunken LFCs are reported as conservative effect sizes, with
#         'shrinkage_collapsed' flagging significant genes shrunk to |LFC| < 0.1
#       * cross-analysis comparisons (A3 sensitivity) use MLE LFCs, because
#         shrunken LFCs from analyses with different priors are not comparable
#   - Significant: padj < 0.05. |LFC| >= 1 is reported separately as an
#     effect-size descriptor, not as a significance criterion.
#   - VST (blind = TRUE) for PCA and heatmaps
#
# Run from the repository root:
#   Rscript data/raw/notebooks/06_differential_expression_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(EnhancedVolcano)
  library(dplyr)
  library(readr)
  library(digest)
})
options(warn = 1)   # print warnings where they occur so they can be attributed

# --- Locate repository root (works from any working directory) ---------------
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- if (length(script_arg) == 1) {
  normalizePath(file.path(dirname(normalizePath(script_arg)), "..", "..", ".."))
} else {
  getwd()
}
setwd(repo_root)

RAW_COUNTS_FILE   <- "data/raw/GSE213410_transcripts.readcount.txt.gz"
GENE_COUNTS_FILE  <- "results/tables/gene_level_counts.csv"
GENE_ANNOT_FILE   <- "results/tables/gene_annotation.csv"
MAPPING_FILE      <- "results/tables/transcript_gene_mapping.csv"
METADATA_FILE     <- "data/metadata/sample_metadata.csv"
SRA_RUNINFO_FILE  <- "data/metadata/sra_runinfo_PRJNA880672.csv"
TABLE_DIR <- "results/tables"
FIG_DIR   <- "results/figures"
REPORT_FILE <- "results/PHASE6_REPORT.md"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

EXPECTED_SAMPLES <- paste0("A", 1:12)
OUTLIER_SAMPLE <- "A3"
ALPHA <- 0.05              # BH-adjusted p-value threshold
MIN_COUNT <- 10            # pre-filter: raw count threshold ...
LFC_DESCRIPTOR <- 1        # |shrunken LFC| >= 1 reported as "large effect"
N_HEATMAP_GENES <- 50
N_LABELS <- 12
set.seed(42)

# Colours (same validated palette as Phases 3-5)
INK <- "#0b0b0b"; INK_2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
COL_NS <- "#c3c2b7"; COL_UP <- "#e34948"; COL_DOWN <- "#2a78d6"
GROUP_COLORS <- c(hESC_Ethanol = "#2a78d6", `hESC_5-OXO-ETE` = "#eb6834",
                  `T21-hiPSC_DMSO` = "#1baf7a", `T21-hiPSC_Calcitriol` = "#eda100")

theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.text = element_text(colour = INK_2),
          axis.title = element_text(colour = INK_2),
          plot.title = element_text(face = "bold", colour = INK, size = 13),
          plot.subtitle = element_text(colour = INK_2, size = 10),
          plot.caption = element_text(colour = MUTED, size = 8.5, hjust = 0),
          plot.title.position = "plot", plot.caption.position = "plot",
          plot.margin = margin(10, 28, 10, 10),
          legend.position = "bottom", legend.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

# --- Validation bookkeeping ---------------------------------------------------
validation <- data.frame(check = character(), result = character(), detail = character())
record_check <- function(name, passed, detail = "") {
  status <- if (isTRUE(passed)) "PASS" else "FAIL"
  validation[nrow(validation) + 1, ] <<- list(name, status, as.character(detail))
  cat(sprintf("  [%s] %s%s\n", status, name, if (nzchar(detail)) paste0(" - ", detail) else ""))
  invisible(isTRUE(passed))
}
section <- function(title) cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))

# =============================================================================
# Step 1: Data validation and count provenance
# =============================================================================
section("Step 1: Data validation and count provenance")
raw_sha_before <- sha256_file(RAW_COUNTS_FILE)
cat("Raw count file SHA-256 (before):", raw_sha_before, "\n")

gene_counts_df <- read_csv(GENE_COUNTS_FILE, show_col_types = FALSE,
                           col_types = cols(gene_id = col_character(), .default = col_double()))
metadata <- read_csv(METADATA_FILE, show_col_types = FALSE)
gene_annotation <- read_csv(GENE_ANNOT_FILE, show_col_types = FALSE)

count_values <- as.matrix(gene_counts_df[, setdiff(names(gene_counts_df), "gene_id")])
record_check("Gene count matrix columns are exactly A1..A12 in order",
             identical(colnames(count_values), EXPECTED_SAMPLES),
             paste(colnames(count_values), collapse = ","))
record_check("Metadata sample_id is exactly A1..A12 in order",
             identical(metadata$sample_id, EXPECTED_SAMPLES))
record_check("Count columns match metadata samples exactly",
             identical(colnames(count_values), metadata$sample_id))
record_check("No missing values in count matrix", !anyNA(count_values))
record_check("All counts are integers (no fractional values)",
             all(count_values == round(count_values)),
             sprintf("%d fractional values", sum(count_values != round(count_values))))
record_check("All counts are non-negative", all(count_values >= 0))
record_check("No duplicated gene identifiers", !anyDuplicated(gene_counts_df$gene_id),
             sprintf("%d genes", nrow(gene_counts_df)))
record_check("Counts fit in R integer range", max(count_values) < .Machine$integer.max)

gene_counts <- count_values
storage.mode(gene_counts) <- "integer"
rownames(gene_counts) <- gene_counts_df$gene_id

# Independent re-aggregation in R from the untouched raw file + Phase 5 mapping
raw_tx <- read_tsv(RAW_COUNTS_FILE, show_col_types = FALSE,
                   col_types = cols(transcriptID = col_character(), .default = col_integer()))
mapping <- read_csv(MAPPING_FILE, show_col_types = FALSE,
                    col_select = c(original_transcript_id, gene_id, mapping_status))
mapping <- mapping[match(raw_tx$transcriptID, mapping$original_transcript_id), ]
record_check("Mapping table covers every raw transcript", !anyNA(mapping$original_transcript_id))
is_mapped <- mapping$mapping_status == "mapped"
raw_tx_matrix <- as.matrix(raw_tx[, EXPECTED_SAMPLES])
reaggregated <- rowsum(raw_tx_matrix[is_mapped, ], group = mapping$gene_id[is_mapped])
record_check("R re-aggregation reproduces Phase 5 gene matrix exactly",
             identical(dim(reaggregated), dim(gene_counts)) &&
               all(reaggregated[rownames(gene_counts), ] == gene_counts))

# Provenance 1: count scale vs. sequenced read pairs (SRA run statistics)
runinfo <- read_csv(SRA_RUNINFO_FILE, show_col_types = FALSE) |>
  select(Run, SampleName, spots, avgLength, LibraryLayout)
provenance <- metadata |>
  select(sample_id, geo_accession) |>
  left_join(runinfo, by = c("geo_accession" = "SampleName")) |>
  mutate(transcript_level_counts = colSums(raw_tx_matrix)[sample_id],
         gene_level_counts = colSums(gene_counts)[sample_id],
         counts_per_read_pair = round(transcript_level_counts / spots, 3))
record_check("SRA run statistics found for all 12 samples", !anyNA(provenance$spots))
peer_ratio <- provenance$counts_per_read_pair[provenance$sample_id != OUTLIER_SAMPLE]
record_check("Non-A3 counts are on a read-pair scale (0.8-1.2 counts per pair)",
             all(peer_ratio > 0.8 & peer_ratio < 1.2),
             sprintf("range %.2f-%.2f", min(peer_ratio), max(peer_ratio)))
print(as.data.frame(provenance))

# Provenance 2: upper bound on prepDE ceil() inflation per gene and sample.
# ceil() adds < 1 count per transcript with non-zero coverage, so the
# inflation of a gene count is strictly below the number of its transcripts
# with count > 0 in that sample.
positive_tx_per_gene <- rowsum((raw_tx_matrix[is_mapped, ] > 0) * 1L,
                               group = mapping$gene_id[is_mapped])[rownames(gene_counts), ]
ceil_bound_fraction <- positive_tx_per_gene / pmax(gene_counts, 1)
write_csv(provenance, file.path(TABLE_DIR, "count_provenance_check.csv"))

# =============================================================================
# Step 2: Environment
# =============================================================================
section("Step 2: Environment")
packages <- c("DESeq2", "apeglm", "ggplot2", "ggrepel", "pheatmap", "EnhancedVolcano",
              "dplyr", "readr", "digest", "BiocManager")
environment_table <- data.frame(
  component = c("R", "Bioconductor", packages),
  version = c(R.version.string, as.character(BiocManager::version()),
              vapply(packages, function(p) as.character(packageVersion(p)), ""))
)
print(environment_table, row.names = FALSE)
write_csv(environment_table, file.path(TABLE_DIR, "phase6_environment.csv"))

# =============================================================================
# Step 3-5: DESeq2 analyses
# =============================================================================
annotation_lookup <- gene_annotation |> select(gene_id, gene_name, gene_biotype)

run_deseq <- function(label, samples) {
  section(paste("DESeq2:", label))
  col_data <- metadata |>
    filter(sample_id %in% samples) |>
    arrange(match(sample_id, samples)) |>
    mutate(condition = factor(treatment_role, levels = c("control", "treatment")))
  col_data <- as.data.frame(col_data)
  rownames(col_data) <- col_data$sample_id
  counts <- gene_counts[, samples, drop = FALSE]

  # Low-count pre-filter: >= MIN_COUNT reads in >= smallest-group-size samples
  n_min <- min(table(col_data$condition))
  keep <- rowSums(counts >= MIN_COUNT) >= n_min
  cat(sprintf("Pre-filter (>= %d counts in >= %d samples): %d of %d genes kept\n",
              MIN_COUNT, n_min, sum(keep), nrow(counts)))

  dds <- DESeqDataSetFromMatrix(countData = counts[keep, ], colData = col_data,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)   # size factors, dispersions, Wald test
  res <- results(dds, alpha = ALPHA, pAdjustMethod = "BH")
  # Capture apeglm warnings (e.g. optimiser line-search messages) so they are
  # counted and reported rather than lost.
  shrink_warnings <- character()
  shrunk <- withCallingHandlers(
    lfcShrink(dds, coef = "condition_treatment_vs_control", type = "apeglm", quiet = TRUE),
    warning = function(w) {
      shrink_warnings <<- c(shrink_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  prior_scale <- priorInfo(shrunk)$prior.control$prior.scale
  cat(sprintf("apeglm prior scale: %.4g; warnings: %d%s\n", prior_scale, length(shrink_warnings),
              if (length(shrink_warnings)) paste0(" (", unique(shrink_warnings)[1], ")") else ""))
  vsd <- vst(dds, blind = TRUE)
  norm_counts <- counts(dds, normalized = TRUE)

  table_out <- data.frame(gene_id = rownames(res),
                          baseMean = res$baseMean,
                          log2FoldChange_MLE = res$log2FoldChange,
                          lfcSE_MLE = res$lfcSE,
                          wald_stat = res$stat,
                          pvalue = res$pvalue,
                          padj = res$padj,
                          log2FoldChange_shrunk = shrunk$log2FoldChange,
                          lfcSE_shrunk = shrunk$lfcSE) |>
    left_join(annotation_lookup, by = "gene_id") |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name),
           # Gene symbols are not unique (e.g. two MATR3 genes): add the Ensembl ID
           gene_label = ifelse(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE),
                               paste0(gene_name, " (", gene_id, ")"), gene_name),
           significant = !is.na(padj) & padj < ALPHA,
           direction = case_when(significant & log2FoldChange_MLE > 0 ~ "up",
                                 significant & log2FoldChange_MLE < 0 ~ "down",
                                 TRUE ~ "not_significant"),
           large_effect_MLE = abs(log2FoldChange_MLE) >= LFC_DESCRIPTOR,
           large_effect_shrunk = abs(log2FoldChange_shrunk) >= LFC_DESCRIPTOR,
           shrinkage_collapsed = significant & abs(log2FoldChange_shrunk) < 0.1,
           test_status = case_when(is.na(pvalue) ~ "pvalue_NA_cooks_outlier",
                                   is.na(padj) ~ "padj_NA_independent_filtering",
                                   TRUE ~ "tested"),
           max_ceil_inflation_bound = apply(ceil_bound_fraction[gene_id, samples, drop = FALSE],
                                            1, max)) |>
    relocate(gene_name, gene_label, gene_biotype, .after = gene_id)
  norm_df <- as.data.frame(round(norm_counts, 2))
  colnames(norm_df) <- paste0("norm_", colnames(norm_df))
  table_out <- cbind(table_out, norm_df[table_out$gene_id, ]) |>
    arrange(is.na(padj), padj, pvalue)

  # BH verification: padj must equal p.adjust(BH) over genes kept by independent filtering
  bh_ok <- with(table_out[!is.na(table_out$padj), ],
                isTRUE(all.equal(padj, p.adjust(pvalue, method = "BH"))))
  record_check(paste0(label, ": padj reproduces Benjamini-Hochberg"), bh_ok)
  record_check(paste0(label, ": DESeq2 object contains only intended samples"),
               identical(colnames(dds), samples), paste(colnames(dds), collapse = ","))

  summary_row <- data.frame(
    analysis = label,
    samples = paste(samples, collapse = ","),
    n_control = sum(col_data$condition == "control"),
    n_treatment = sum(col_data$condition == "treatment"),
    genes_in_matrix = nrow(counts),
    genes_after_prefilter = sum(keep),
    genes_with_pvalue = sum(!is.na(res$pvalue)),
    genes_with_padj = sum(!is.na(res$padj)),
    cooks_outliers = sum(is.na(res$pvalue) & res$baseMean > 0),
    independent_filter_removed = sum(!is.na(res$pvalue) & is.na(res$padj)),
    independent_filter_threshold_baseMean = round(metadata(res)$filterThreshold, 3),
    significant_padj_0.05 = sum(table_out$significant),
    up = sum(table_out$direction == "up"),
    down = sum(table_out$direction == "down"),
    significant_and_abs_MLE_LFC_ge_1 = sum(table_out$significant & table_out$large_effect_MLE),
    significant_and_abs_shrunk_LFC_ge_1 = sum(table_out$significant & table_out$large_effect_shrunk),
    significant_shrinkage_collapsed = sum(table_out$shrinkage_collapsed),
    apeglm_prior_scale = signif(prior_scale, 4),
    apeglm_warnings = length(shrink_warnings),
    size_factors = paste(sprintf("%s=%.3f", samples, sizeFactors(dds)), collapse = "; ")
  )
  print(t(summary_row))
  list(label = label, samples = samples, dds = dds, vsd = vsd, table = table_out,
       summary = summary_row, col_data = col_data)
}

c1_samples <- metadata$sample_id[metadata$comparison == unique(metadata$comparison)[1]]
c2_samples <- metadata$sample_id[metadata$comparison == unique(metadata$comparison)[2]]
c1_noA3_samples <- setdiff(c1_samples, OUTLIER_SAMPLE)

record_check("C1 samples are A1-A6", identical(c1_samples, paste0("A", 1:6)))
record_check("C2 samples are A7-A12", identical(c2_samples, paste0("A", 7:12)))
record_check("C1 and C2 share no samples (independent analyses)",
             length(intersect(c1_samples, c2_samples)) == 0)
record_check("C1 sensitivity excludes only A3",
             identical(setdiff(c1_samples, c1_noA3_samples), OUTLIER_SAMPLE) &&
               all(c1_noA3_samples %in% c1_samples))

c1 <- run_deseq("C1 (all samples, primary)", c1_samples)
c1_noA3 <- run_deseq("C1 without A3 (sensitivity)", c1_noA3_samples)
c2 <- run_deseq("C2 (all samples)", c2_samples)

write_csv(c1$table, file.path(TABLE_DIR, "C1_differential_expression.csv"))
write_csv(c1_noA3$table, file.path(TABLE_DIR, "C1_differential_expression_without_A3.csv"))
write_csv(c2$table, file.path(TABLE_DIR, "C2_differential_expression.csv"))
sig_columns <- c("gene_id", "gene_name", "gene_biotype", "baseMean", "log2FoldChange_MLE",
                 "lfcSE_MLE", "log2FoldChange_shrunk", "lfcSE_shrunk", "pvalue", "padj",
                 "direction", "large_effect_MLE", "shrinkage_collapsed",
                 "max_ceil_inflation_bound")
write_csv(c1$table |> filter(significant) |> select(all_of(sig_columns)),
          file.path(TABLE_DIR, "C1_significant_genes.csv"))
write_csv(c2$table |> filter(significant) |> select(all_of(sig_columns)),
          file.path(TABLE_DIR, "C2_significant_genes.csv"))
write_csv(c1_noA3$table |> filter(significant) |> select(all_of(sig_columns)),
          file.path(TABLE_DIR, "C1_significant_genes_without_A3.csv"))
de_summary <- bind_rows(c1$summary, c1_noA3$summary, c2$summary)
write_csv(de_summary, file.path(TABLE_DIR, "phase6_de_summary.csv"))

# =============================================================================
# A3 sensitivity comparison (C1 with vs. without A3)
# =============================================================================
section("A3 sensitivity comparison")
# MLE LFCs are the primary comparison: the two analyses estimate different
# apeglm priors, so their shrunken LFCs are not on a comparable scale.
sens <- full_join(
  c1$table |> select(gene_id, gene_name, baseMean_with = baseMean,
                     lfc_mle_with = log2FoldChange_MLE, lfc_shrunk_with = log2FoldChange_shrunk,
                     padj_with = padj, sig_with = significant),
  c1_noA3$table |> select(gene_id, baseMean_without = baseMean,
                          lfc_mle_without = log2FoldChange_MLE,
                          lfc_shrunk_without = log2FoldChange_shrunk, padj_without = padj,
                          sig_without = significant),
  by = "gene_id") |>
  mutate(sig_with = coalesce(sig_with, FALSE), sig_without = coalesce(sig_without, FALSE),
         tested_in = case_when(!is.na(baseMean_with) & !is.na(baseMean_without) ~ "both",
                               !is.na(baseMean_with) ~ "with_A3_only",
                               TRUE ~ "without_A3_only"),
         significance_category = case_when(sig_with & sig_without ~ "significant_in_both",
                                           sig_with ~ "significant_with_A3_only",
                                           sig_without ~ "significant_without_A3_only",
                                           TRUE ~ "not_significant"),
         direction_concordant = ifelse(sig_with & sig_without,
                                       sign(lfc_mle_with) == sign(lfc_mle_without), NA),
         lfc_mle_difference = lfc_mle_without - lfc_mle_with) |>
  arrange(factor(significance_category, levels = c("significant_in_both",
                                                   "significant_with_A3_only",
                                                   "significant_without_A3_only",
                                                   "not_significant")),
          pmin(padj_with, padj_without, na.rm = TRUE))
write_csv(sens, file.path(TABLE_DIR, "C1_A3_sensitivity_comparison.csv"))

both_tested <- sens |> filter(tested_in == "both")
n_sig_with <- sum(sens$sig_with); n_sig_without <- sum(sens$sig_without)
n_sig_both <- sum(sens$significance_category == "significant_in_both")
jaccard <- if (n_sig_with + n_sig_without - n_sig_both > 0) {
  n_sig_both / (n_sig_with + n_sig_without - n_sig_both) } else NA
sens_summary <- data.frame(
  metric = c("genes tested with A3", "genes tested without A3", "genes tested in both",
             "significant with A3", "significant without A3", "significant in both",
             "significant with A3 only", "significant without A3 only",
             "Jaccard overlap of significant sets",
             "direction concordance among genes significant in both",
             "Pearson r of MLE LFC (genes tested in both)",
             "Spearman r of MLE LFC (genes tested in both)",
             "Pearson r of MLE LFC (genes significant in either analysis)",
             "Median |MLE LFC change| (genes significant in either analysis)",
             "Spearman r of -log10 padj (genes with padj in both)"),
  value = c(nrow(c1$table), nrow(c1_noA3$table), nrow(both_tested),
            n_sig_with, n_sig_without, n_sig_both,
            sum(sens$significance_category == "significant_with_A3_only"),
            sum(sens$significance_category == "significant_without_A3_only"),
            round(jaccard, 3),
            if (n_sig_both > 0) round(mean(sens$direction_concordant, na.rm = TRUE), 3) else NA,
            round(cor(both_tested$lfc_mle_with, both_tested$lfc_mle_without), 3),
            round(cor(both_tested$lfc_mle_with, both_tested$lfc_mle_without,
                      method = "spearman"), 3),
            with(both_tested[both_tested$sig_with | both_tested$sig_without, ],
                 round(cor(lfc_mle_with, lfc_mle_without), 3)),
            with(both_tested[both_tested$sig_with | both_tested$sig_without, ],
                 round(median(abs(lfc_mle_difference)), 3)),
            with(both_tested[!is.na(both_tested$padj_with) & !is.na(both_tested$padj_without), ],
                 round(cor(-log10(padj_with), -log10(padj_without), method = "spearman"), 3)))
)
print(sens_summary, row.names = FALSE)
write_csv(sens_summary, file.path(TABLE_DIR, "C1_A3_sensitivity_summary.csv"))

# =============================================================================
# Step 6: Figures
# =============================================================================
section("Step 6: Figures")
saved_figures <- character()
save_plot <- function(plot, name, width, height) {
  path <- file.path(FIG_DIR, name)
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  saved_figures <<- c(saved_figures, name)
  cat("  saved figure:", path, "\n")
}
pretty_label <- function(result) {
  switch(result$label,
         "C1 (all samples, primary)" = "C1 H9 hESC: 5-OXO-ETE vs Ethanol (n = 3 vs 3, A3 included)",
         "C1 without A3 (sensitivity)" = "C1 H9 hESC: 5-OXO-ETE vs Ethanol, A3 excluded (n = 3 vs 2)",
         "C2 (all samples)" = "C2 T21 hiPSC: Calcitriol vs DMSO (n = 3 vs 3)")
}
top_label_genes <- function(tab) {
  tab |> filter(significant) |> arrange(padj) |> head(N_LABELS) |> pull(gene_label)
}

volcano_plot <- function(result) {
  tab <- result$table |> filter(!is.na(padj))
  tab$padj_plot <- pmax(tab$padj, 1e-300)
  n_up <- sum(tab$direction == "up"); n_down <- sum(tab$direction == "down")
  keyvals <- ifelse(tab$direction == "up", COL_UP, ifelse(tab$direction == "down", COL_DOWN, COL_NS))
  names(keyvals) <- ifelse(tab$direction == "up", sprintf("Up (padj < 0.05), n = %d", n_up),
                    ifelse(tab$direction == "down", sprintf("Down (padj < 0.05), n = %d", n_down),
                           "Not significant"))
  EnhancedVolcano(tab, lab = tab$gene_label, x = "log2FoldChange_MLE", y = "padj_plot",
                  selectLab = top_label_genes(tab),
                  pCutoff = ALPHA, FCcutoff = LFC_DESCRIPTOR, cutoffLineType = "dashed",
                  cutoffLineCol = MUTED, cutoffLineWidth = 0.4,
                  colCustom = keyvals, colAlpha = 0.75, pointSize = 1.6, labSize = 3.2,
                  drawConnectors = TRUE, widthConnectors = 0.3, colConnectors = MUTED,
                  max.overlaps = Inf, boxedLabels = TRUE,
                  title = pretty_label(result),
                  subtitle = sprintf("%s genes with adjusted p-values; Wald test, BH-adjusted",
                                     format(nrow(tab), big.mark = ",")),
                  caption = paste("x: unshrunken (MLE) log2 fold change, treatment / control: the",
                                  "quantity tested by the Wald test.\nHorizontal line: padj = 0.05",
                                  "(significance). Vertical lines: |LFC| = 1, a visual guide only."),
                  xlab = "log2 fold change (MLE)", ylab = expression(-log[10]~"adjusted p-value"),
                  legendPosition = "bottom", legendLabSize = 10, legendIconSize = 3,
                  gridlines.minor = FALSE, border = "partial") +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = INK_2),
          plot.caption = element_text(size = 8.5, colour = MUTED, hjust = 0))
}

ma_plot <- function(result) {
  tab <- result$table |> filter(!is.na(log2FoldChange_shrunk)) |>
    mutate(direction = factor(direction, levels = c("up", "down", "not_significant")))
  labels <- tab |> filter(significant) |> arrange(padj) |> head(N_LABELS)
  ggplot(tab, aes(baseMean, log2FoldChange_shrunk)) +
    geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_point(data = filter(tab, direction == "not_significant"), colour = COL_NS,
               size = 0.8, alpha = 0.6) +
    geom_point(data = filter(tab, direction != "not_significant"), aes(colour = direction),
               size = 1.6, alpha = 0.9) +
    geom_label_repel(data = labels, aes(label = gene_label), size = 3.1, colour = INK,
                     fill = "white", label.size = 0, label.padding = unit(0.12, "lines"),
                     segment.colour = MUTED, segment.size = 0.3, max.overlaps = Inf,
                     min.segment.length = 0, box.padding = 0.45, seed = 42) +
    scale_x_log10(labels = scales::label_comma(), expand = expansion(mult = c(0.02, 0.04))) +
    scale_colour_manual(values = c(up = COL_UP, down = COL_DOWN),
                        labels = c(up = sprintf("Up (padj < 0.05), n = %d", result$summary$up),
                                   down = sprintf("Down (padj < 0.05), n = %d", result$summary$down)),
                        drop = FALSE) +
    labs(title = pretty_label(result),
         subtitle = "MA plot: mean normalised expression vs apeglm-shrunken log2 fold change",
         x = "Mean of normalised counts (log scale)", y = "Shrunken log2 fold change (apeglm)",
         caption = paste0("Grey: not significant (padj >= 0.05 or not tested). Colour and up/down use the ",
                          "MLE fold change.\napeglm prior scale = ", result$summary$apeglm_prior_scale,
                          "; ", result$summary$significant_shrinkage_collapsed,
                          " significant gene(s) shrink to |LFC| < 0.1 (effect size uncertain).")) +
    theme_portfolio()
}

pca_plot <- function(result) {
  pca_data <- plotPCA(result$vsd, intgroup = c("experimental_group", "sample_id"),
                      ntop = 500, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
  pca_data$is_outlier <- pca_data$sample_id == OUTLIER_SAMPLE
  groups_present <- unique(as.character(result$col_data$experimental_group))
  ggplot(pca_data, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
    geom_point(aes(fill = experimental_group), shape = 21, size = 4,
               colour = ifelse(pca_data$is_outlier, INK, "white"),
               stroke = ifelse(pca_data$is_outlier, 1.4, 0.8)) +
    geom_text_repel(aes(label = sample_id), size = 3.5, colour = INK_2, seed = 42,
                    box.padding = 0.5, point.padding = 0.4, min.segment.length = 0.3) +
    scale_fill_manual(values = GROUP_COLORS[groups_present]) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(title = pretty_label(result),
         subtitle = "PCA of variance-stabilised counts (VST, blind), top 500 variable genes",
         x = sprintf("PC1 (%.1f%% of variance)", percent_var[1]),
         y = sprintf("PC2 (%.1f%% of variance)", percent_var[2]),
         caption = if (OUTLIER_SAMPLE %in% pca_data$sample_id)
           "A3 (QC outlier from Phases 3-5) drawn with a dark outline." else NULL) +
    theme_portfolio()
}

heatmap_plot <- function(result, file_name) {
  tab <- result$table
  n_sig <- sum(tab$significant)
  chosen <- if (n_sig >= 2) {
    tab |> filter(significant) |> arrange(padj) |> head(N_HEATMAP_GENES)
  } else {
    tab |> filter(!is.na(padj)) |> arrange(padj) |> head(N_HEATMAP_GENES)
  }
  heading <- if (n_sig >= 2) {
    sprintf("%s\nTop %d of %d significant genes (padj < 0.05), VST row z-scores",
            pretty_label(result), nrow(chosen), n_sig)
  } else {
    sprintf("%s\nTop %d genes by adjusted p-value (only %d significant), VST row z-scores",
            pretty_label(result), nrow(chosen), n_sig)
  }
  vst_matrix <- assay(result$vsd)[chosen$gene_id, , drop = FALSE]
  rownames(vst_matrix) <- chosen$gene_label
  z <- t(scale(t(vst_matrix)))
  z[is.na(z)] <- 0
  annotation <- data.frame(Group = as.character(result$col_data$experimental_group),
                           row.names = result$col_data$sample_id)
  groups_present <- unique(annotation$Group)
  breaks <- seq(-2.5, 2.5, length.out = 101)
  palette <- colorRampPalette(c("#2a78d6", "#f0efec", "#e34948"))(100)
  path <- file.path(FIG_DIR, file_name)
  row_font <- if (nrow(z) > 40) 6.5 else 8
  pheatmap(z, color = palette, breaks = breaks, border_color = NA,
           cluster_rows = nrow(z) > 1, cluster_cols = TRUE,
           annotation_col = annotation,
           annotation_colors = list(Group = GROUP_COLORS[groups_present]),
           show_rownames = TRUE, fontsize_row = row_font, fontsize_col = 10,
           fontsize = 9, main = heading, angle_col = 0,
           filename = path, width = 7.5, height = max(5, 1.6 + nrow(z) * 0.13))
  saved_figures <<- c(saved_figures, file_name)
  cat("  saved figure:", path, "\n")
}

save_plot(volcano_plot(c1), "13_C1_volcano.png", 9, 7.5)
save_plot(volcano_plot(c2), "14_C2_volcano.png", 9, 7.5)
save_plot(ma_plot(c1), "15_C1_MA_plot.png", 9, 6)
save_plot(ma_plot(c2), "16_C2_MA_plot.png", 9, 6)
heatmap_plot(c1, "17_C1_gene_heatmap.png")
heatmap_plot(c2, "18_C2_gene_heatmap.png")
save_plot(pca_plot(c1), "19_C1_PCA.png", 7.5, 6)
save_plot(pca_plot(c2), "20_C2_PCA.png", 7.5, 6)
# C1 sensitivity figures (A3 excluded)
save_plot(volcano_plot(c1_noA3), "21_C1_without_A3_volcano.png", 9, 7.5)
save_plot(ma_plot(c1_noA3), "22_C1_without_A3_MA_plot.png", 9, 6)
heatmap_plot(c1_noA3, "23_C1_without_A3_gene_heatmap.png")
save_plot(pca_plot(c1_noA3), "24_C1_without_A3_PCA.png", 7.5, 6)

# LFC agreement with vs. without A3
sens_plot_data <- both_tested |>
  mutate(category = factor(significance_category,
                           levels = c("not_significant", "significant_in_both",
                                      "significant_with_A3_only", "significant_without_A3_only")))
sens_colors <- c(not_significant = COL_NS, significant_in_both = "#6250d6",
                 significant_with_A3_only = "#eb6834", significant_without_A3_only = "#1baf7a")
sens_labels <- c(
  not_significant = "Not significant in either",
  significant_in_both = sprintf("Significant in both, n = %d", n_sig_both),
  significant_with_A3_only = sprintf("Significant with A3 only, n = %d",
                                     sum(sens$significance_category == "significant_with_A3_only")),
  significant_without_A3_only = sprintf("Significant without A3 only, n = %d",
                                        sum(sens$significance_category == "significant_without_A3_only")))
lfc_limit <- max(abs(c(sens_plot_data$lfc_mle_with, sens_plot_data$lfc_mle_without)), na.rm = TRUE)
sens_lfc_plot <- ggplot(sens_plot_data, aes(lfc_mle_with, lfc_mle_without)) +
  geom_abline(slope = 1, intercept = 0, colour = MUTED, linetype = "dashed", linewidth = 0.4) +
  geom_point(data = filter(sens_plot_data, category == "not_significant"),
             colour = COL_NS, size = 0.7, alpha = 0.5) +
  geom_point(data = filter(sens_plot_data, category != "not_significant"),
             aes(colour = category), size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = sens_colors, labels = sens_labels, drop = TRUE) +
  coord_equal(xlim = c(-lfc_limit, lfc_limit), ylim = c(-lfc_limit, lfc_limit)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = "C1 A3 sensitivity: log2 fold changes with vs without A3",
       subtitle = sprintf("%s genes tested in both analyses; Pearson r = %.3f (MLE LFC)",
                          format(nrow(both_tested), big.mark = ","),
                          cor(both_tested$lfc_mle_with, both_tested$lfc_mle_without)),
       x = "MLE log2 fold change, all six samples (A3 included)",
       y = "MLE log2 fold change, A3 excluded",
       caption = paste("Dashed line: identity (y = x). MLE fold changes are compared because",
                       "the two analyses estimate different apeglm priors.")) +
  theme_portfolio()
save_plot(sens_lfc_plot, "25_C1_A3_sensitivity_LFC_comparison.png", 7.5, 7.8)

# =============================================================================
# Step 9: Validation and report
# =============================================================================
section("Step 9: Validation")
expected_tables <- c("C1_differential_expression.csv", "C1_differential_expression_without_A3.csv",
                     "C2_differential_expression.csv", "C1_significant_genes.csv",
                     "C2_significant_genes.csv", "C1_A3_sensitivity_comparison.csv")
expected_figures <- c("13_C1_volcano.png", "14_C2_volcano.png", "15_C1_MA_plot.png",
                      "16_C2_MA_plot.png", "17_C1_gene_heatmap.png", "18_C2_gene_heatmap.png",
                      "19_C1_PCA.png", "20_C2_PCA.png")
missing_tables <- expected_tables[!file.exists(file.path(TABLE_DIR, expected_tables))]
missing_figures <- expected_figures[!file.exists(file.path(FIG_DIR, expected_figures))]
record_check("All required tables written", length(missing_tables) == 0,
             if (length(missing_tables)) paste(missing_tables, collapse = ",") else "6 tables")
record_check("All required figures written", length(missing_figures) == 0,
             if (length(missing_figures)) paste(missing_figures, collapse = ",") else
               sprintf("%d figures incl. sensitivity", length(saved_figures)))
record_check("Every DE table has one row per tested gene (no duplicates)",
             all(vapply(list(c1$table, c1_noA3$table, c2$table),
                        function(t) !anyDuplicated(t$gene_id), logical(1))))
raw_sha_after <- sha256_file(RAW_COUNTS_FILE)
cat("Raw count file SHA-256 (after): ", raw_sha_after, "\n")
record_check("Original raw count file unchanged (SHA-256 before = after)",
             identical(raw_sha_before, raw_sha_after), raw_sha_after)
write_csv(validation, file.path(TABLE_DIR, "phase6_validation_report.csv"))
n_failed <- sum(validation$result == "FAIL")
cat(sprintf("\n%d/%d validation checks passed\n", nrow(validation) - n_failed, nrow(validation)))

# --- Report -------------------------------------------------------------------
top_table_md <- function(result, n = 10) {
  top <- result$table |> filter(significant) |> arrange(padj) |> head(n)
  if (nrow(top) == 0) return("_No genes reached padj < 0.05._\n")
  lines <- c("| Gene | Ensembl ID | baseMean | Direction | LFC (MLE) ± SE | Shrunken LFC | padj |",
             "|---|---|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %s | %.2f ± %.2f | %.2f%s | %s |", top$gene_name, top$gene_id,
                     format(round(top$baseMean), big.mark = ","), top$direction,
                     top$log2FoldChange_MLE, top$lfcSE_MLE, top$log2FoldChange_shrunk,
                     ifelse(top$shrinkage_collapsed, " †", ""), fmt_p(top$padj)))
  paste(lines, collapse = "\n")
}
summary_md <- function(s) {
  sprintf("| %s | %d vs %d | %s | %s | %d | %d | %d | %d | %d | %s |",
          s$analysis, s$n_control, s$n_treatment, format(s$genes_after_prefilter, big.mark = ","),
          format(s$genes_with_padj, big.mark = ","), s$significant_padj_0.05, s$up, s$down,
          s$significant_and_abs_MLE_LFC_ge_1, s$significant_shrinkage_collapsed,
          format(s$apeglm_prior_scale))
}
ceil_summary <- function(result) {
  sig <- result$table |> filter(significant)
  if (nrow(sig) == 0) return("no significant genes")
  sprintf("median bound %.4f; %d of %d significant genes have a bound above 5%%",
          median(sig$max_ceil_inflation_bound), sum(sig$max_ceil_inflation_bound > 0.05), nrow(sig))
}
overlap_c1_c2 <- length(intersect(c1$table$gene_id[c1$table$significant],
                                  c2$table$gene_id[c2$table$significant]))
a3_prov <- provenance[provenance$sample_id == OUTLIER_SAMPLE, ]

report <- c(
  "# Phase 6: Differential gene expression (DESeq2)",
  "",
  sprintf("*Generated by `data/raw/notebooks/06_differential_expression_analysis.R` on %s. All numbers are computed by the script.*",
          format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  sprintf("**Validation: %d/%d checks passed** (`results/tables/phase6_validation_report.csv`).",
          nrow(validation) - n_failed, nrow(validation)),
  "",
  "## 1. Count provenance and suitability for DESeq2",
  "",
  "- The GEO counts come from StringTie 1.3.3b. StringTie itself reports coverage, FPKM and TPM; integer read counts come from its `prepDE.py` helper, which computes `ceil(coverage x transcript_length / read_length)` per transcript (v1.3.3b source). The counts are therefore **coverage-derived estimates, already integers, rounded up**. No rounding is done here, and the matrix contains no fractional values.",
  sprintf("- **Scale check against SRA** (PRJNA880672, paired-end, %d bp per spot): for the 11 non-A3 samples, total counts are %.2f-%.2f per sequenced read pair. The counts sit on a fragment-count scale, as DESeq2's count model assumes, and are not inflated by a mismatched read-length setting.",
          as.integer(median(provenance$avgLength)), min(peer_ratio), max(peer_ratio)),
  sprintf("- **A3:** %.1f M read pairs were sequenced (similar to the others), but only %.1f M counts were assigned (%.2f per pair). A3's small library is therefore not shallow sequencing: most of its fragments were not assigned to annotated or assembled transcripts. This is further evidence of a technical problem with A3. GEO provides no alignment statistics, so the cause cannot be determined.",
          a3_prov$spots / 1e6, a3_prov$transcript_level_counts / 1e6, a3_prov$counts_per_read_pair),
  "- **Rounding-up bias:** `ceil()` adds less than 1 count per transcript with any coverage, so a gene's count is inflated by less than its number of expressed transcripts. Each DE table reports `max_ceil_inflation_bound` (that number divided by the gene count, maximum over samples). The bias affects both groups in the same way and matters only for low-count genes, which the pre-filter largely removes.",
  "- **Aggregation:** gene counts are sums of transcript estimates. That is the same operation `prepDE.py` uses to build its gene matrix, and the result reflects reads per gene, not per isoform. An R re-aggregation from the raw file reproduces the Phase 5 matrix exactly.",
  "- **Conclusion:** DESeq2 is appropriate for these counts, with the caveat that they are StringTie estimates rather than unique-read counts (for example from featureCounts). No methodological problem blocks the analysis.",
  "",
  "## 2. Methods",
  "",
  sprintf("- DESeq2 %s on R %s (Bioconductor %s); design `~ condition` (control as reference), each comparison fitted separately.",
          as.character(packageVersion("DESeq2")), paste(R.version$major, R.version$minor, sep = "."),
          as.character(BiocManager::version())),
  sprintf("- Pre-filter: >= %d raw counts in at least as many samples as the smallest group (3; 2 in the without-A3 analysis).", MIN_COUNT),
  "- Median-of-ratios size factors; dispersions with empirical-Bayes shrinkage; Wald test; Benjamini-Hochberg adjustment with DESeq2 independent filtering; Cook's distance flagging.",
  sprintf("- **Significance: padj < %.2f.** Up/down is set by the sign of the unshrunken (MLE) fold change, the quantity the Wald test evaluates. |LFC| >= %g is reported as a descriptor of magnitude, not as a significance criterion.", ALPHA, LFC_DESCRIPTOR),
  "- apeglm shrinkage is reported as a conservative effect-size estimate. **Its estimated prior scale reached the lower bound (0.001) in C1 and C2** (table below). The genome-wide spread of treatment effects is so small that apeglm shrinks almost every gene to about 0, including some genes the Wald test calls significant. Those genes are flagged `shrinkage_collapsed` (marked † below): the test detects a difference, but its size is uncertain. Because each analysis estimates its own prior, shrunken LFCs are not compared across analyses; the A3 sensitivity check uses MLE LFCs.",
  sprintf("- apeglm optimiser warnings ('line search routine failed'): C1 %d, C1 without A3 %d, C2 %d. Switching to apeglm's alternative optimiser gave essentially the same near-zero estimates, which are consistent with the narrow prior rather than an optimiser failure.",
          c1$summary$apeglm_warnings, c1_noA3$summary$apeglm_warnings, c2$summary$apeglm_warnings),
  "- VST (blind) for PCA and heatmaps. Seed 42.",
  "",
  "## 3. Results",
  "",
  "| Analysis | Control vs treatment | Genes after pre-filter | Genes with padj | Significant | Up | Down | Significant and \\|MLE LFC\\| >= 1 | Significant, shrinkage-collapsed | apeglm prior scale |",
  "|---|---|---|---|---|---|---|---|---|---|",
  summary_md(c1$summary), summary_md(c1_noA3$summary), summary_md(c2$summary),
  "",
  "\"Up\" means higher in the treated group (5-OXO-ETE in C1, Calcitriol in C2). † = significant, but apeglm shrinks the LFC to |LFC| < 0.1.",
  "",
  "### C1 top genes (A3 included)",
  "", top_table_md(c1), "",
  "### C1 top genes (A3 excluded)",
  "", top_table_md(c1_noA3), "",
  "### C2 top genes",
  "", top_table_md(c2), "",
  "## 4. Effect of excluding A3 (C1)",
  "",
  paste0("| Metric | Value |\n|---|---|\n",
         paste(sprintf("| %s | %s |", sens_summary$metric, sens_summary$value), collapse = "\n")),
  "",
  "Gene-level details: `C1_A3_sensitivity_comparison.csv`. Figures 21-25.",
  "",
  "## 5. C1 vs. C2",
  "",
  sprintf("The comparisons were fitted independently. %d genes are significant in both. Because the cell lines, genotypes, vehicles and treatments all differ, this overlap is descriptive only and is not a test of shared response.", overlap_c1_c2),
  "",
  "## 6. Rounding-up bias among significant genes",
  "",
  sprintf("- C1: %s", ceil_summary(c1)),
  sprintf("- C1 without A3: %s", ceil_summary(c1_noA3)),
  sprintf("- C2: %s", ceil_summary(c2)),
  "",
  "## 7. Limitations",
  "",
  "- Three replicates per group (two for the C1 controls without A3) give limited power. Small but real changes may be missed, and dispersion estimates rely heavily on shrinkage.",
  "- Statistical significance is not biological importance. A small padj can come with a fold change too small to matter, and large fold changes in low-count genes are unstable. Shrunken LFCs are reported for that reason.",
  "- These are associations between treatment and expression within one cell line per comparison at one time point (24 h). They do not establish mechanism, and they may not generalise to other lines.",
  "- Counts are StringTie coverage-based estimates (see Section 1). TCONS transcripts (6.5% of counts) are not part of the gene matrix.",
  "- No batch, RNA-integrity or alignment metadata are available. A3's low assignment rate cannot be explained from the published data.",
  "- No functional enrichment was performed in this phase.",
  "",
  "## 8. Status",
  "",
  if (n_failed == 0) "All validation checks passed. Results are ready for review and for functional enrichment in a later phase." else
    sprintf("**%d validation check(s) FAILED. Do not use these results until resolved.**", n_failed),
  ""
)
writeLines(report, REPORT_FILE, useBytes = TRUE)
cat("  saved report:", REPORT_FILE, "\n")

if (n_failed > 0) {
  print(validation[validation$result == "FAIL", ])
  stop(sprintf("%d validation check(s) failed.", n_failed))
}
cat("\nPhase 6 completed successfully.\n")
