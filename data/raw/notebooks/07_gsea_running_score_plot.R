# =============================================================================
# Phase 7 helper: GSEA running-score figures (figures 28, 31, 34)
# =============================================================================
# Sourced by:
#   data/raw/notebooks/07_functional_enrichment.R        (full analysis run)
#   data/raw/notebooks/07b_regenerate_running_score_figures.R
#                                                        (figures only, from saved tables)
# Keeping the plotting code in one file guarantees both scripts draw the same figure.
#
# Selection rule (fixed after the first run, 2026-11-25):
#   1. Take the single most significant gene set (smallest BH-adjusted p) of the analysis.
#   2. Use that set's database.
#   3. Plot the top `n_sets` significant sets of THAT database (by adjusted p).
# The original code took the top 3 sets across all databases and then kept only
# those from the first set's database, which could leave a single, unlabelled curve.
#
# Labels: each curve is labelled in the legend with its term name, NES and adjusted p,
# taken unchanged from the saved results table. The title is passed to gseaplot2()
# itself (the original version drew it on top of the plot with grid.text(), which
# overlapped the panel).
#
# This file only draws figures. It does not compute or change any statistic.
# =============================================================================

select_running_score_sets <- function(sig_table, n_sets = 3) {
  # sig_table: significant GSEA rows of ONE analysis (columns database, ID, p.adjust, ...)
  sig_table <- sig_table[order(sig_table$p.adjust), , drop = FALSE]
  if (nrow(sig_table) == 0) return(sig_table)
  db <- sig_table$database[1]
  head(sig_table[sig_table$database == db, , drop = FALSE], n_sets)
}

running_score_label <- function(rows) {
  sprintf("%s (NES %.2f, adj. p %s)", rows$Description, rows$NES,
          formatC(rows$p.adjust, format = "e", digits = 1))
}

save_running_score_figure <- function(gsea_obj, rows, title, path,
                                      colors = c("#e34948", "#2a78d6", "#1baf7a")) {
  # gsea_obj: clusterProfiler gseaResult holding geneList, geneSets and a result slot
  #           whose rows include rows$ID. Only the Description of the plotted rows is
  #           replaced by a display label; the object is a local copy.
  ids <- rows$ID
  gsea_obj@result[ids, "Description"] <- running_score_label(rows)
  p <- enrichplot::gseaplot2(gsea_obj, geneSetID = ids, title = title, base_size = 10,
                             color = colors[seq_along(ids)], pvalue_table = FALSE)
  # Put the legend under the running-score panel instead of over the curves
  p[[1]] <- p[[1]] +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical",
                   legend.justification = "left",
                   legend.text = ggplot2::element_text(size = 8.5),
                   plot.title = ggplot2::element_text(face = "bold", size = 11.5),
                   plot.title.position = "plot") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
  png(path, width = 9, height = 7, units = "in", res = 300, bg = "white")
  on.exit(dev.off(), add = TRUE)
  print(p)
  invisible(path)
}
