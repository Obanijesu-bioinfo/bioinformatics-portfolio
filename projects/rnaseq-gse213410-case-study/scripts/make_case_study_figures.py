"""
Case-study presentation figures (improved versions; originals are kept).

Creates new figure files next to the originals in
projects/rnaseq-gse213410-case-study/figures/. No original figure, analysis
script or Phase 3-8 output is modified. Every plotted value is read from a
saved machine-readable output; nothing is re-estimated except the Phase 5
gene-level PCA, which is recomputed with the Phase 5 method and checked
against the variance percentages Phase 5 recorded.

Outputs (the script stops if any already exists, so nothing is overwritten):
  figures/11b_gene_level_pca_accessible.png
  figures/13b_C1_volcano_phase8_targets.png
  figures/14b_C2_volcano_phase8_targets.png
  figures/25b_C1_A3_sensitivity_LFC_comparison_denominator.png
  figures/27b_C1_GSEA_theme_summary.png
  figures/30b_C2_GSEA_theme_summary.png
  tables/gsea_theme_assignments.csv   one row per significant gene set
  tables/gsea_theme_summary.csv       counts per theme and direction
  tables/gsea_histone_threshold_sensitivity.csv
                                      histone-dominated set counts at 40/50/60%
To regenerate in place, delete only the files listed above (never the original
figures), or pass a separate output folder as the first argument.

Inputs (read only):
  results/tables/gene_level_counts.csv, gene_level_A3_sensitivity.csv
  data/metadata/sample_metadata.csv
  results/tables/C1_differential_expression.csv, C2_differential_expression.csv
  results/tables/C1_A3_sensitivity_comparison.csv, C1_A3_sensitivity_summary.csv
  results/tables/phase6_de_summary.csv
  results/tables/phase7_C1_GSEA_results.csv, phase7_C2_GSEA_results.csv
  results/tables/phase8_gene_interpretation.csv

GSEA theme summary method (descriptive only):
  1. Every significant gene set (BH-adjusted p < 0.05) is kept; none is
     dropped or hand-picked.
  2. "Histone-gene dominated": at least half of the set's leading-edge
     (core enrichment) genes are histone genes (symbol pattern HISTONE_RE).
     This flags the overlap between many chromatin, infection and disease
     sets that share the same histone genes.
  3. Otherwise the first matching keyword rule in THEME_RULES (applied to the
     set name, in the listed order) gives the theme; no match gives "Other".
     Short tokens are bounded so they cannot match inside other words (for
     example "RNA" is case-sensitive, so it does not match "external").
  Every assignment is written to tables/gsea_theme_assignments.csv. The themes
  are approximate groupings for presentation, not biological classifications.
  The 50% histone threshold is a choice; the sensitivity table shows how many
  sets it captures at 40% and 60%.

Run from anywhere:
    python projects/rnaseq-gse213410-case-study/scripts/make_case_study_figures.py
"""

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

REPO_ROOT = Path(__file__).resolve().parents[3]
CASE = REPO_ROOT / "projects/rnaseq-gse213410-case-study"
FIG = CASE / "figures"
TAB_OUT = CASE / "tables"
if len(sys.argv) > 1:                      # optional draft folder, e.g. for layout checks
    FIG = TAB_OUT = Path(sys.argv[1]).resolve()
    FIG.mkdir(parents=True, exist_ok=True)
T = REPO_ROOT / "results/tables"

OUTPUTS = {
    "pca": FIG / "11b_gene_level_pca_accessible.png",
    "c1_volcano": FIG / "13b_C1_volcano_phase8_targets.png",
    "c2_volcano": FIG / "14b_C2_volcano_phase8_targets.png",
    "sens": FIG / "25b_C1_A3_sensitivity_LFC_comparison_denominator.png",
    "c1_gsea": FIG / "27b_C1_GSEA_theme_summary.png",
    "c2_gsea": FIG / "30b_C2_GSEA_theme_summary.png",
    "theme_rows": TAB_OUT / "gsea_theme_assignments.csv",
    "theme_summary": TAB_OUT / "gsea_theme_summary.csv",
    "histone_sensitivity": TAB_OUT / "gsea_histone_threshold_sensitivity.csv",
}
existing = [str(p) for p in OUTPUTS.values() if p.exists()]
if existing:
    sys.exit("Refusing to overwrite existing outputs:\n  " + "\n  ".join(existing))
TAB_OUT.mkdir(exist_ok=True)

# ------------------------------------------------------------------ style
# Same palette as Phases 5-7 (validated for colour-vision deficiency there);
# larger type than the originals for readability when scaled in a README.
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#6b6a66"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
COL_UP, COL_DOWN, COL_NS = "#e34948", "#2a78d6", "#c3c2b7"
GROUP_COLORS = {"hESC_Ethanol": "#2a78d6", "hESC_5-OXO-ETE": "#eb6834",
                "T21-hiPSC_DMSO": "#1baf7a", "T21-hiPSC_Calcitriol": "#eda100"}
GROUP_LABELS = {"hESC_Ethanol": "Ethanol (control)", "hESC_5-OXO-ETE": "5-OXO-ETE",
                "T21-hiPSC_DMSO": "DMSO (control)", "T21-hiPSC_Calcitriol": "Calcitriol"}
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 12, "text.color": INK, "axes.labelcolor": INK_2, "axes.labelsize": 12.5,
    "axes.titlesize": 14, "axes.titleweight": "semibold", "axes.titlelocation": "left",
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "legend.frameon": False, "legend.fontsize": 11,
})


def save(fig, key):
    with open(OUTPUTS[key], "xb") as handle:          # exclusive: never overwrite
        fig.savefig(handle, format="png", facecolor="white")
    plt.close(fig)
    print("saved", OUTPUTS[key].name)


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail else ""))
    if not ok:
        sys.exit(f"Check failed: {name}")


def footnote(fig, text, y=-0.02):
    fig.text(0.01, y, text, ha="left", va="top", fontsize=10, color=MUTED, wrap=True)


# ------------------------------------------------------------------ 11b PCA
# Phase 5 method, unchanged: per comparison, DESeq2-style median-of-ratios size
# factors (genes > 0 in every sample), keep normalised count >= 10 in >= 3
# samples, log2(normalised + 1), PCA of the 500 most variable genes (centred,
# unscaled).
meta = pd.read_csv(REPO_ROOT / "data/metadata/sample_metadata.csv").set_index("sample_id")
counts = pd.read_csv(T / "gene_level_counts.csv", index_col="gene_id")
recorded = pd.read_csv(T / "gene_level_A3_sensitivity.csv")


def pca_scores(members):
    sub = counts[members]
    pos = sub[(sub > 0).all(axis=1)]
    logc = np.log(pos.astype(float))
    sf = np.exp(logc.sub(logc.mean(axis=1), axis=0).median())
    norm = sub / sf
    keep = (norm >= 10).sum(axis=1) >= 3
    logm = np.log2(norm[keep] + 1)
    top = logm.var(axis=1).sort_values(ascending=False).index[:500]
    x = logm.loc[top].T
    pca = PCA(n_components=min(x.shape))
    sc = pd.DataFrame(pca.fit_transform(x), index=x.index).iloc[:, :2]
    sc.columns = ["PC1", "PC2"]
    return sc, pca.explained_variance_ratio_[:2] * 100, int(keep.sum())


print("Figure 11b: gene-level PCA")
panels = []
for comp, rec_row, title in [("C1", 0, "C1: H9 hESC, 5-OXO-ETE vs ethanol"),
                             ("C2", 2, "C2: T21 hiPSC, Calcitriol vs DMSO")]:
    members = [s for s in meta.index if meta.at[s, "comparison"].startswith(comp + "_")]
    sc, var, n_genes = pca_scores(members)
    rec = recorded.iloc[rec_row]
    check(f"{comp} PCA matches Phase 5 record",
          abs(var[0] - rec.PC1_pct_variance) < 0.01 and abs(var[1] - rec.PC2_pct_variance) < 0.01
          and n_genes == rec.genes_after_filter,
          f"PC1 {var[0]:.2f}%, PC2 {var[1]:.2f}%, {n_genes} genes")
    panels.append((comp, title, members, sc, var))

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, (comp, title, members, sc, var) in zip(axes, panels):
    for s in members:
        group = meta.at[s, "experimental_group"]
        is_ctrl = meta.at[s, "treatment_role"] == "control"
        ax.scatter(sc.at[s, "PC1"], sc.at[s, "PC2"], s=130 if is_ctrl else 150,
                   marker="o" if is_ctrl else "^", color=GROUP_COLORS[group],
                   edgecolor=INK if s == "A3" else "white", linewidth=2.2 if s == "A3" else 1.2,
                   zorder=3)
        ax.annotate(s + (" (QC outlier)" if s == "A3" else ""), (sc.at[s, "PC1"], sc.at[s, "PC2"]),
                    xytext=(9, 6), textcoords="offset points", fontsize=11.5,
                    color=INK if s == "A3" else INK_2, fontweight="bold" if s == "A3" else "normal")
    ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    ax.axvline(0, color=AXIS, lw=0.8, zorder=1)
    ax.margins(0.22)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% of variance)")
    ax.set_title(title)
    groups = list(dict.fromkeys(meta.loc[members, "experimental_group"]))
    handles = [Line2D([], [], marker="o" if meta.loc[meta.experimental_group == g, "treatment_role"].iloc[0]
                      == "control" else "^", linestyle="", markersize=11,
                      markerfacecolor=GROUP_COLORS[g], markeredgecolor="white", label=GROUP_LABELS[g])
               for g in groups]
    if "A3" in members:
        handles.append(Line2D([], [], marker="o", linestyle="", markersize=11, markerfacecolor="white",
                              markeredgecolor=INK, markeredgewidth=2.2, label="A3: QC outlier"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)
fig.suptitle("Gene-level PCA, each comparison separately", x=0.01, ha="left",
             fontsize=15, fontweight="semibold")
fig.tight_layout()
footnote(fig, "Method (Phase 5): per comparison, median-of-ratios size factors; genes with normalised count "
              ">= 10 in >= 3 samples; log2(normalised count + 1);\nPCA of the 500 most variable genes "
              "(centred, unscaled). Circles: controls; triangles: treated. Not the VST-based PCA of Phase 6.",
         y=-0.01)
save(fig, "pca")

# ------------------------------------------------------------------ volcanoes
interp = pd.read_csv(T / "phase8_gene_interpretation.csv")
de_summary = pd.read_csv(T / "phase6_de_summary.csv")
BIOTYPE_TEXT = {"processed_pseudogene": "processed pseudogene",
                "transcribed_processed_pseudogene": "transcribed pseudogene",
                "protein_coding": "protein coding"}


def target_label(row):
    sym = row.fact_target_symbol
    if row.fact_ensembl_gene_id == "ENSG00000286192":
        kind = "novel readthrough model, no known function"
    else:
        kind = BIOTYPE_TEXT.get(row.fact_gene_biotype, row.fact_gene_biotype.replace("_", " "))
    if sym == "MATR3":
        role = "primary record" if "primary" in row.fact_record_role else "overlapping gene model"
        head = f"MATR3 {row.fact_ensembl_gene_id} ({role})"
    else:
        head = f"{sym}: {kind}"
    conf = row.conf_biological_interpretation.split(" (")[0].lower()
    return f"{head}\nInterpretation confidence: {conf}"


def volcano(comp, de_file, summary_label, title, offsets, key):
    print(f"Figure {key}")
    de = pd.read_csv(T / de_file)
    tab = de[de.padj.notna()].copy()
    tab["y"] = -np.log10(tab.padj.clip(lower=1e-300))
    rec = de_summary[de_summary.analysis == summary_label].iloc[0]
    n_up, n_down = int((tab.direction == "up").sum()), int((tab.direction == "down").sum())
    check(f"{comp} volcano counts match phase6_de_summary",
          len(tab) == rec.genes_with_padj and n_up == rec.up and n_down == rec.down,
          f"{len(tab)} with padj, {n_up} up, {n_down} down")
    targets = interp[interp.fact_comparison == comp]
    check(f"{comp} Phase 8 targets all significant in this table",
          bool(tab.set_index("gene_id").loc[targets.fact_ensembl_gene_id, "significant"].all()),
          f"{len(targets)} records")

    fig, ax = plt.subplots(figsize=(11, 8.6))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.25)
    ns = tab[tab.direction == "not_significant"]
    ax.scatter(ns.log2FoldChange_MLE, ns.y, s=7, color=COL_NS, alpha=0.6, linewidth=0, zorder=1)
    for d, col, mk in [("up", COL_UP, "o"), ("down", COL_DOWN, "o")]:
        sub = tab[tab.direction == d]
        ax.scatter(sub.log2FoldChange_MLE, sub.y, s=26, color=col, alpha=0.85, marker=mk,
                   edgecolor="white", linewidth=0.4, zorder=2)
    ax.axhline(-np.log10(0.05), color=MUTED, ls="--", lw=0.8)
    for x in (-1, 1):
        ax.axvline(x, color=AXIS, ls=":", lw=0.9)
    ax.set_ylim(-0.04 * tab.y.max(), tab.y.max() * 1.08)
    t_idx = tab.set_index("gene_id")
    for _, row in targets.iterrows():
        g = t_idx.loc[row.fact_ensembl_gene_id]
        ax.scatter(g.log2FoldChange_MLE, g.y, s=110, facecolor="none", edgecolor=INK, linewidth=1.6, zorder=4)
        fx, fy = offsets[row.fact_ensembl_gene_id]      # label position, axes fraction
        ax.annotate(target_label(row), (g.log2FoldChange_MLE, g.y), xytext=(fx, fy),
                    textcoords="axes fraction", fontsize=10.5, color=INK,
                    ha="right" if fx > 0.5 else "left", va="center",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=AXIS, linewidth=0.8),
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkB=6), zorder=5)
    ax.set_xlabel("log2 fold change, treated vs control (unshrunken MLE)")
    ax.set_ylabel("−log10 adjusted p-value")
    ax.set_title(title)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", color=COL_UP, markersize=8, label=f"Higher in treated (padj < 0.05), n = {n_up}"),
        Line2D([], [], marker="o", ls="", color=COL_DOWN, markersize=8, label=f"Lower in treated (padj < 0.05), n = {n_down}"),
        Line2D([], [], marker="o", ls="", color=COL_NS, markersize=8, label="Not significant"),
        Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=INK, markeredgewidth=1.6,
               markersize=10, label="Phase 8 target record (labelled)")],
        loc="lower center", bbox_to_anchor=(0.5, 0.085), ncol=2)
    fig.text(0.01, 0.01, f"{len(tab):,} genes with adjusted p-values (Wald test, Benjamini–Hochberg). "
                         "Dashed line: padj = 0.05. Dotted lines: |log2FC| = 1, a visual guide only.\n"
                         "Labels: Ensembl 97 annotation and Phase 8 interpretation confidence. "
                         "An RNA change is not evidence of a protein or functional change.",
             ha="left", va="bottom", fontsize=10, color=MUTED)
    save(fig, key)
    return tab


volcano("C1", "C1_differential_expression.csv", "C1 (all samples, primary)",
        "C1 H9 hESC: 5-OXO-ETE vs ethanol (3 vs 3, A3 included)",
        {"ENSG00000237506": (0.99, 0.95), "ENSG00000127589": (0.99, 0.80),
         "ENSG00000139269": (0.01, 0.92), "ENSG00000286192": (0.01, 0.62)}, "c1_volcano")
volcano("C2", "C2_differential_expression.csv", "C2 (all samples)",
        "C2 T21 hiPSC: Calcitriol vs DMSO (3 vs 3)",
        {"ENSG00000213064": (0.99, 0.84), "ENSG00000280987": (0.99, 0.68), "ENSG00000157214": (0.99, 0.54),
         # EOGT label above its point; overlapping-MATR3 label below the EOGT leader line so lines do not cross labels
         "ENSG00000163378": (0.01, 0.55), "ENSG00000015479": (0.01, 0.21)}, "c2_volcano")

# ------------------------------------------------------------------ 25b sensitivity
print("Figure 25b: A3 sensitivity")
sens = pd.read_csv(T / "C1_A3_sensitivity_comparison.csv")
summ = pd.read_csv(T / "C1_A3_sensitivity_summary.csv").set_index("metric")["value"]
both = sens[sens.tested_in == "both"]
cat_all = sens.significance_category.value_counts()
cat_plot = both.significance_category.value_counts()
r = both.lfc_mle_with.corr(both.lfc_mle_without)
check("Denominator: genes in both result tables", len(both) == summ["genes tested in both"] == 19594, f"{len(both)}")
check("All plotted genes have an MLE LFC in both analyses",
      bool(both[["lfc_mle_with", "lfc_mle_without"]].notna().all().all()))
check("Pearson r matches summary", round(r, 3) == summ["Pearson r of MLE LFC (genes tested in both)"], f"{r:.4f}")
for k, m in [("significant_in_both", "significant in both"), ("significant_with_A3_only", "significant with A3 only"),
             ("significant_without_A3_only", "significant without A3 only")]:
    check(f"{m} matches summary", cat_all.get(k, 0) == summ[m], f"{cat_all.get(k, 0)}")
n_padj_both = int(both[["padj_with", "padj_without"]].notna().all(axis=1).sum())

styles = {"significant_in_both": ("#6250d6", "o", "Significant in both"),
          "significant_with_A3_only": ("#eb6834", "s", "Significant with A3 only"),
          "significant_without_A3_only": ("#1baf7a", "^", "Significant without A3 only")}
fig, ax = plt.subplots(figsize=(9.5, 10.6))
fig.subplots_adjust(left=0.11, right=0.97, top=0.95, bottom=0.255)
lim =float(np.ceil(max(both.lfc_mle_with.abs().max(), both.lfc_mle_without.abs().max()) + 0.5))
ax.plot([-lim, lim], [-lim, lim], color=MUTED, ls="--", lw=0.9, zorder=1)
ns = both[both.significance_category == "not_significant"]
ax.scatter(ns.lfc_mle_with, ns.lfc_mle_without, s=6, color=COL_NS, alpha=0.5, linewidth=0, zorder=1)
handles = [Line2D([], [], marker="o", ls="", color=COL_NS, markersize=7,
                  label=f"Not significant in either, n = {cat_plot.get('not_significant', 0):,}")]
for k, (col, mk, lab) in styles.items():
    sub = both[both.significance_category == k]
    ax.scatter(sub.lfc_mle_with, sub.lfc_mle_without, s=40, color=col, marker=mk,
               edgecolor="white", linewidth=0.5, zorder=3 if k != "significant_without_A3_only" else 2)
    total = cat_all.get(k, 0)
    extra = "" if len(sub) == total else f" ({total - len(sub)} more not plotted*)"
    handles.append(Line2D([], [], marker=mk, ls="", color=col, markersize=9,
                          label=f"{lab}, n = {len(sub)} plotted{extra}"))
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.set_xlabel("MLE log2 fold change, A3 included (primary)")
ax.set_ylabel("MLE log2 fold change, A3 excluded (exploratory)")
ax.set_title("C1: log2 fold changes with vs without sample A3")
ax.text(0.02, 0.98, f"Plotted: {len(both):,} genes retained after\nprefiltering in both analyses\n"
                    f"Pearson r (MLE log2FC) = {r:.3f}",
        transform=ax.transAxes, va="top", fontsize=11.5, color=INK,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=AXIS))
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.11, 0.085), ncol=1)
fig.text(0.01, 0.01, f"Significant = padj < 0.05. Totals: {int(summ['significant with A3'])} significant with A3, "
              f"{int(summ['significant without A3'])} without A3.\n*Genes retained after prefiltering only "
              f"without A3 ({(sens.tested_in == 'without_A3_only').sum():,}) have no with-A3 estimate and "
              f"cannot be plotted. {n_padj_both:,} of the\n{len(both):,} plotted genes have an adjusted p-value "
              "in both analyses. Dashed line: y = x. This correlation is for gene-level\nfold changes, not the "
              "separate Phase 7 GSEA sensitivity comparison.",
         ha="left", va="bottom", fontsize=10, color=MUTED)
save(fig, "sens")

# ------------------------------------------------------------------ GSEA theme summaries
HISTONE_RE = re.compile(r"^(H1-\d|H1-10|H2AC\d|H2BC\d|H2AZ\d|H2AX|H2AJ|H2AW|H2BW|H3C\d|H3-\d|H3-3|H4C\d|H4-\d|MACROH2A)")
HISTONE_THRESHOLD = 0.5                      # fraction of leading-edge genes; see sensitivity table
SENSITIVITY_THRESHOLDS = (0.4, 0.5, 0.6)
INFECTION = "Infection, immune and disease-labelled sets"
# RNA species as whole words, case-sensitive: RNA, mRNA, pre-mRNA, rRNA, snRNA, microRNAs, piRNAs, ...
# "(?-i:...)" switches off case-insensitive matching so "RNA" cannot match inside "external"; RNA
# polymerase (transcription) is excluded.
RNA_WORD = r"(?-i:(?<![A-Za-z])(?:m|r|t|sn|sno|mi|micro|pi|si|lnc)?RNAs?\b)(?! [Pp]olymerase)"
THEME_RULES = [
    # Pathogen names first, so e.g. "Dengue Virus Genome ... Replication" is not read as DNA replication.
    # Viral-protein / host-pathogen context terms (deliberately not a general "host" trigger):
    #   "Rev"          the HIV-1 Rev protein (case-sensitive whole word, so "REV-ERBA" is not matched);
    #   "host nucleus" viral material entering the host nucleus (Reactome influenza pathway).
    # In the Phase 7 source tables both terms occur only in viral pathways. They group "Nuclear import of
    # Rev protein", "Interactions of Rev with host cellular proteins" and "Transport of Ribonucleoproteins
    # into the Host Nucleus" with their explicitly viral sibling pathways, which share the same host
    # nuclear-pore leading-edge genes.
    (INFECTION, r"virus|viral|HIV|HCMV|SARS|Dengue|DENV|NS1|NEP/NS2|influenza|parasite|Leishmania|"
                r"(?-i:\bRev\b)|\bhost nucleus\b"),
    # "helicase activity" and "nucleic acid conformation isomerase activity" are not matched here: their
    # leading edges mix RNA helicases (DDX/DHX) with DNA replication/repair helicases, so they stay in "Other".
    ("DNA replication and repair",
     r"repair|replication|recombination|HDR|HRR|NHEJ|double.strand|excision|glycosylase|depurination|"
     r"depyrimidination|BRCA|PALB2|RAD51|D-loop|Holliday|strand exchange|telomer|synthesis of DNA|ORC complex|"
     r"pre-replicative|acting on DNA|DNA metabolic|origins"),
    # Septins and tubulin folding are cytoskeletal, not specifically cell-cycle terms, so they stay in "Other".
    ("Cell cycle, mitosis and chromosome organisation",
     r"cell cycle|mitotic|M Phase|G2/M|S Phase|chromosom|chromatid|kinetochore|spindle|centrosom|APC/C|PLK1|"
     r"AURKA|nuclear division|organelle fission|metaphase|anaphase|prophase|prometaphase|cohesion|NuMA|"
     r"microtubule organizing|Meio"),
    ("Chromatin and epigenetic regulation",
     r"chromatin|nucleosome|histone|epigenetic|methylat|PRC2|MLL|WDR5|CENP|silencing|retroelement|SIRT1|NoRC|"
     r"B-WICH|ERCC6|deacetylat|acetylat|demethylat|maternal to zygotic|HOX"),
    # (HIV Rev pathways are matched by the pathogen rule above; "REV-ERBA" matches neither rule.)
    ("RNA processing, splicing and transport",
     RNA_WORD + r"|splic|ribonucleoprotein|snRNP|Sm-like|ribosome biogenesis|nuclear pore|nucleocytoplasmic|"
     r"nucleic acid transport|polyadenylation|3'-end|intronless|SLBP|\bU2\b"),
    ("Extracellular matrix, collagen and adhesion",
     r"extracellular matrix|ECM|collagen|external encapsulating|proteoglycan|glycosaminoglycan|GAG|B3GAT3|"
     r"B3GALT6|B4GALT7|EXT1|EXT2|cadherin|CDH1|adhesion|L1CAM|NCAM|IgSF CAM|\bL1$|exostoses"),
    # Plasma-membrane and cell-surface terms only: not mitochondrial membranes, nuclear receptors, or
    # ligand / regulator activities ("receptor ligand activity", "signaling receptor regulator activity").
    ("Membrane channels, transporters and receptors",
     r"channel|transporter|transport across|SLC|plasma membrane|synaptic membrane|postsynaptic|membrane raft|"
     r"membrane microdomain|side of membrane|basal part of cell|receptor complex|receptor (protein )?(tyrosine )?"
     r"kinase|receptor activity|coreceptor|synap|neurotransmitter|potassium|detection of|stimulus|Sensory"),
    (INFECTION,
     r"lupus|Alcoholism|phagocyt|interleukin|\bT cell|PD-1|PD-L1|immun|FCGR|Fcgamma|pyroptosis|necroptosis|"
     r"disease|defective|cancer|carcinoma|amyloid"),
]
THEME_ORDER = (["Histone-gene dominated (shared histone genes)"]
               + list(dict.fromkeys(t for t, _ in THEME_RULES[1:])) + [INFECTION, "Other"])
THEME_ORDER = list(dict.fromkeys(THEME_ORDER))


def histone_fraction(core):
    genes = [g for g in str(core).split("/") if g]
    return sum(bool(HISTONE_RE.match(g)) for g in genes) / len(genes) if genes else 0.0


def assign_theme(desc, core):
    frac = histone_fraction(core)
    if frac >= HISTONE_THRESHOLD:
        return THEME_ORDER[0], frac
    for theme, pattern in THEME_RULES:
        if re.search(pattern, desc, flags=re.IGNORECASE):
            return theme, frac
    return "Other", frac


assign_rows, summary_rows, sensitivity_rows = [], [], []
for comp, key, title in [("C1", "c1_gsea", "C1 GSEA themes: H9 hESC, 5-OXO-ETE vs ethanol (primary, A3 included)"),
                         ("C2", "c2_gsea", "C2 GSEA themes: T21 hiPSC, Calcitriol vs DMSO")]:
    print(f"Figure {key}")
    g = pd.read_csv(T / f"phase7_{comp}_GSEA_results.csv")
    sig = g[g.significant == True].copy()        # noqa: E712 (column is boolean)
    check(f"{comp} significant sets all have padj < 0.05", bool((sig["p.adjust"] < 0.05).all()), f"{len(sig)} sets")
    check(f"{comp} no duplicate gene sets", not sig.duplicated(["database", "ID"]).any())
    fractions = sig.core_enrichment.map(histone_fraction)
    for thr in SENSITIVITY_THRESHOLDS:
        hit = fractions >= thr
        sensitivity_rows.append({"comparison": comp, "threshold": thr, "histone_dominated_sets": int(hit.sum()),
                                 "higher_in_treatment": int((hit & (sig.direction == "higher in treatment")).sum()),
                                 "lower_in_treatment": int((hit & (sig.direction == "lower in treatment")).sum())})
    sig[["theme", "histone_fraction_leading_edge"]] = [assign_theme(d, c) for d, c in zip(sig.Description, sig.core_enrichment)]
    sig["comparison"] = comp
    assign_rows.append(sig[["comparison", "database", "ID", "Description", "direction", "NES", "p.adjust",
                            "setSize", "histone_fraction_leading_edge", "theme"]])
    tally = sig.groupby(["theme", "direction"]).size().unstack(fill_value=0).reindex(THEME_ORDER, fill_value=0)
    for col in ("higher in treatment", "lower in treatment"):
        if col not in tally:
            tally[col] = 0
    check(f"{comp} theme counts sum to all significant sets", int(tally.values.sum()) == len(sig))
    for theme in THEME_ORDER:
        summary_rows.append({"comparison": comp, "theme": theme,
                             "sets_higher_in_treatment": int(tally.at[theme, "higher in treatment"]),
                             "sets_lower_in_treatment": int(tally.at[theme, "lower in treatment"])})

    fig, ax = plt.subplots(figsize=(11, 6.2))
    y = np.arange(len(THEME_ORDER))[::-1]
    up = tally["higher in treatment"].values
    dn = tally["lower in treatment"].values
    bars_up = ax.barh(y, up, color=COL_UP, height=0.62, edgecolor="white", linewidth=1)
    bars_dn = ax.barh(y, -dn, color=COL_DOWN, height=0.62, edgecolor="white", linewidth=1)
    check(f"{comp} figure bars equal the summary counts",
          [int(round(b.get_width())) for b in bars_up] == [int(v) for v in up]
          and [int(round(-b.get_width())) for b in bars_dn] == [int(v) for v in dn])
    xmax = max(up.max(), dn.max()) * 1.25
    for yi, u, d in zip(y, up, dn):
        if u:
            ax.text(u + xmax * 0.015, yi, str(u), va="center", ha="left", fontsize=11, color=INK)
        if d:
            ax.text(-d - xmax * 0.015, yi, str(d), va="center", ha="right", fontsize=11, color=INK)
    ax.axvline(0, color=INK_2, lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(THEME_ORDER, fontsize=11.5)
    ax.set_xlim(-xmax, xmax)
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(abs(t))) for t in ticks])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Number of significant gene sets: lower in treated (left), higher in treated (right)")
    fig.suptitle(title, x=0.01, ha="left", fontsize=15, fontweight="semibold")
    ax.set_title(f"All {len(sig)} significant gene sets (GO, KEGG, Reactome; BH-adjusted p < 0.05): "
                 f"{int(up.sum())} higher and {int(dn.sum())} lower in treated",
                 fontsize=11.5, fontweight="normal", color=INK_2, loc="left")
    footnote(fig, "Themes are approximate groupings, not biological classifications. Gene sets overlap heavily, "
                  "so counts are not independent signals.\nHistone-gene dominated: >= 50% of a set's "
                  "leading-edge genes are histone genes (a chosen threshold). Other sets: first matching keyword\n"
                  "rule on the set name. Every assignment: tables/gsea_theme_assignments.csv. Full plot of the "
                  f"top sets: figure {'27' if comp == 'C1' else '30'}. Hypothesis-generating only.", y=-0.01)
    save(fig, key)

with open(OUTPUTS["theme_rows"], "x", newline="", encoding="utf-8") as fh:
    pd.concat(assign_rows).to_csv(fh, index=False)
print("saved", OUTPUTS["theme_rows"].name)
with open(OUTPUTS["theme_summary"], "x", newline="", encoding="utf-8") as fh:
    pd.DataFrame(summary_rows).to_csv(fh, index=False)
print("saved", OUTPUTS["theme_summary"].name)
with open(OUTPUTS["histone_sensitivity"], "x", newline="", encoding="utf-8") as fh:
    pd.DataFrame(sensitivity_rows).to_csv(fh, index=False)
print("saved", OUTPUTS["histone_sensitivity"].name)
