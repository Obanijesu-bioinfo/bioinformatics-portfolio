#!/usr/bin/env Rscript
# =============================================================================
# Phase 7b: regenerate GSEA running-score figures 28, 31 and 34 from saved results
# =============================================================================
# Why this script exists
#   The first Phase 7 run (2026-11-24) produced running-score figures whose title
#   overlapped the plot and which could show a single unlabelled curve. Rerunning
#   the whole enrichment analysis only to redraw three figures is unnecessary and
#   would query live KEGG again. This script redraws the figures WITHOUT running
#   any enrichment test.
#
# What it does
#   1. Rebuilds each GSEA ranking exactly as Step 4 of 07_functional_enrichment.R:
#      Phase 6 genes with non-NA p-value and a unique Entrez ID (from the saved
#      phase7_gene_id_mapping.csv), ranked by decreasing Wald statistic.
#   2. Loads the gene sets of the plotted terms from the same annotation packages
#      clusterProfiler used (org.Hs.eg.db GO, reactome.db).
#   3. VERIFIES, for every plotted gene set, that the ranking + gene set reproduce
#      the saved table: set size, enrichment score, peak rank and leading-edge genes.
#      Any mismatch stops the script before a figure is written.
#   4. Builds a gseaResult object whose result slot is the saved table rows
#      (NES, p-values and adjusted p-values are copied, never recomputed) and
#      draws the figure with the shared helper 07_gsea_running_score_plot.R.
#
# It reads tables only; it never writes or modifies any results table.
#
# Run from the repository root:
#   Rscript data/raw/notebooks/07b_regenerate_running_score_figures.R [output_dir]
# output_dir defaults to results/figures.
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(DOSE)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(readr)
})
options(warn = 1)

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- if (length(script_arg) == 1) {
  normalizePath(file.path(dirname(normalizePath(script_arg)), "..", "..", ".."))
} else getwd()
setwd(repo_root)
source("data/raw/notebooks/07_gsea_running_score_plot.R")

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else "results/figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
TABLE_DIR <- "results/tables"

# Same definitions as 07_functional_enrichment.R
FIGURES <- list(
  C1 = list(file = "results/tables/C1_differential_expression.csv", fig = 28,
            label = "C1 H9 hESC: 5-OXO-ETE vs Ethanol (A3 included, primary)"),
  C2 = list(file = "results/tables/C2_differential_expression.csv", fig = 31,
            label = "C2 T21 hiPSC: Calcitriol vs DMSO"),
  C1_without_A3 = list(file = "results/tables/C1_differential_expression_without_A3.csv", fig = 34,
                       label = "C1 H9 hESC: 5-OXO-ETE vs Ethanol (A3 excluded, exploratory)")
)
GSEA_EXPONENT <- 1   # clusterProfiler default, used in the original run

mapping <- read_csv(file.path(TABLE_DIR, "phase7_gene_id_mapping.csv"),
                    col_types = cols(.default = col_character()))
mapping <- mapping |> filter(mapping_status == "unique") |> select(gene_id, entrez_id)

# Gene sets as clusterProfiler builds them (internal helpers of the installed versions)
gene_set_source <- list(
  GO_BP = function() clusterProfiler:::get_GO_data(org.Hs.eg.db, "BP", "ENTREZID"),
  GO_MF = function() clusterProfiler:::get_GO_data(org.Hs.eg.db, "MF", "ENTREZID"),
  GO_CC = function() clusterProfiler:::get_GO_data(org.Hs.eg.db, "CC", "ENTREZID"),
  Reactome = function() ReactomePA:::get_Reactome_DATA("human")
)

checks <- data.frame()
for (nm in names(FIGURES)) {
  cfg <- FIGURES[[nm]]
  cat(sprintf("\n== %s (figure %d) ==\n", nm, cfg$fig))

  # 1. Ranking, identical to Step 4 of the main script
  de <- read_csv(cfg$file, show_col_types = FALSE)
  ranking <- de |> filter(!is.na(pvalue)) |> inner_join(mapping, by = "gene_id") |>
    filter(is.finite(wald_stat)) |> arrange(desc(wald_stat))
  gene_list <- setNames(ranking$wald_stat, ranking$entrez_id)
  stopifnot(!anyDuplicated(names(gene_list)))
  cat(sprintf("ranked genes: %d\n", length(gene_list)))

  # 2. Saved results and the sets to plot
  saved <- read_csv(file.path(TABLE_DIR, sprintf("phase7_%s_GSEA_results.csv", nm)),
                    show_col_types = FALSE)
  rows <- select_running_score_sets(filter(saved, significant))
  db <- rows$database[1]
  if (!db %in% names(gene_set_source)) stop("No gene-set loader for database ", db)
  env <- gene_set_source[[db]]()
  sets <- get("PATHID2EXTID", envir = env)[rows$ID]
  if (any(vapply(sets, is.null, TRUE))) stop("Gene set(s) missing for ", nm, ": ", paste(rows$ID, collapse = ","))

  # 3. Verify the saved table is reproduced by ranking + gene set
  for (i in seq_len(nrow(rows))) {
    gs <- intersect(sets[[i]], names(gene_list))
    sc <- DOSE:::gseaScores(gene_list, gs, exponent = GSEA_EXPONENT, fortify = TRUE)
    peak <- which.max(abs(sc$runningScore))
    es <- sc$runningScore[peak]
    hits <- names(gene_list)[sc$position == 1]
    core_entrez <- if (es >= 0) hits[match(hits, names(gene_list)) <= peak] else
      hits[match(hits, names(gene_list)) >= peak]
    core_symbol <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, core_entrez, "SYMBOL", "ENTREZID",
                                                multiVals = "first"))
    core_symbol[is.na(core_symbol)] <- core_entrez[is.na(core_symbol)]
    saved_core <- strsplit(rows$core_enrichment[i], "/")[[1]]
    row_check <- data.frame(
      analysis = nm, figure = cfg$fig, database = db, ID = rows$ID[i],
      setSize_saved = rows$setSize[i], setSize_rebuilt = length(gs),
      ES_saved = rows$enrichmentScore[i], ES_rebuilt = es,
      ES_abs_diff = abs(rows$enrichmentScore[i] - es),
      # clusterProfiler reports the peak position from the top of the ranking for
      # positive ES and from the bottom (n - peak + 1) for negative ES
      rank_saved = rows$rank[i],
      rank_rebuilt = if (es >= 0) peak else length(gene_list) - peak + 1,
      core_genes_saved = length(saved_core), core_genes_rebuilt = length(core_symbol),
      core_identical_set = setequal(saved_core, core_symbol))
    row_check$pass <- with(row_check, setSize_saved == setSize_rebuilt & ES_abs_diff < 1e-12 &
                             rank_saved == rank_rebuilt & core_identical_set)
    checks <- rbind(checks, row_check)
  }
  print(checks[checks$analysis == nm, ], row.names = FALSE)
  if (!all(checks$pass[checks$analysis == nm])) {
    stop("Rebuilt running scores do not reproduce the saved GSEA table for ", nm,
         ". No figure written.")
  }

  # 4. gseaResult whose statistics are the saved rows, unchanged
  res_df <- as.data.frame(rows |> select(ID, Description, setSize, enrichmentScore, NES,
                                         pvalue, p.adjust, rank, leading_edge, core_enrichment))
  rownames(res_df) <- res_df$ID
  obj <- new("gseaResult", result = res_df, geneSets = sets, geneList = gene_list,
             keytype = "ENTREZID", setType = db, organism = "Homo sapiens",
             params = list(exponent = GSEA_EXPONENT, pAdjustMethod = "BH"), readable = FALSE)
  path <- file.path(OUT_DIR, sprintf("%d_%s_GSEA_running_scores.png", cfg$fig, nm))
  save_running_score_figure(obj, rows,
                            sprintf("GSEA running scores (%s): %s", gsub("_", " ", db), cfg$label), path)
  cat("saved figure:", path, "\n")
}

cat(sprintf("\nVerification: %d/%d plotted gene sets reproduce the saved table (setSize, ES, rank, leading edge).\n",
            sum(checks$pass), nrow(checks)))
cat("Max |ES saved - ES rebuilt|:", format(max(checks$ES_abs_diff), digits = 3), "\n")
