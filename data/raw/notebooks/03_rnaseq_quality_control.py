"""
Phase 3: RNA-seq quality control and exploratory data analysis
Dataset: GSE213410_transcripts.readcount.txt.gz (transcript-level raw read counts)

Steps
  1. Total counts (library size) per sample
  2. Library size bar chart
  3. Detected transcripts per sample (plus a depth-matched comparison)
  4. Sample-to-sample correlation matrix and heatmap
  5. Principal Component Analysis (PCA)
  6. PCA scatter plot (+ scree plot)
  7. Focused investigation of sample A3 as a potential outlier
  8. Figures  -> results/figures/
  9. Tables   -> results/tables/

Constraints honoured
  - The original compressed dataset is only read, never written.
  - A3 is kept in every analysis (it is investigated, not removed).
  - No differential expression analysis is performed.

Run from anywhere:
    python data/raw/notebooks/03_rnaseq_quality_control.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # render to files only, no GUI window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sklearn.decomposition import PCA

# ------------------------------------------------------------------
# Paths and parameters
# ------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_FILE = REPO_ROOT / "data/raw/GSE213410_transcripts.readcount.txt.gz"
FIG_DIR = REPO_ROOT / "results/figures"
TABLE_DIR = REPO_ROOT / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

ID_COLUMN = "transcriptID"
FOCUS_SAMPLE = "A3"       # sample under investigation
MIN_CPM = 1.0             # expression filter: CPM >= 1 ...
MIN_SAMPLES = 3           # ... in at least this many samples
N_TOP_VARIABLE = 500      # transcripts used for PCA (common convention)
RANDOM_SEED = 42          # makes the downsampling reproducible

# ------------------------------------------------------------------
# Plot style: clean scientific look, recessive axes, consistent colors.
# Blue marks ordinary samples; orange is reserved for the focus sample A3.
# ------------------------------------------------------------------
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"
BLUE_LIGHT = "#86b6ef"
ORANGE = "#eb6834"
SEQ_BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "axes.titlesize": 12,
    "axes.titleweight": "semibold",
    "axes.titlelocation": "left",
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "legend.frameon": False,
})


def section(title):
    print("\n" + "=" * 64)
    print(title)
    print("=" * 64)


def save_fig(fig, name):
    path = FIG_DIR / name
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print(f"  saved figure: {path.relative_to(REPO_ROOT)}")


def save_table(df, name, index=True):
    path = TABLE_DIR / name
    df.to_csv(path, index=index)
    print(f"  saved table:  {path.relative_to(REPO_ROOT)}")


def sample_colors(samples):
    return [ORANGE if s == FOCUS_SAMPLE else BLUE for s in samples]


def bold_focus_ticklabel(ax, axis="x"):
    labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for lbl in labels:
        if lbl.get_text() == FOCUS_SAMPLE:
            lbl.set_fontweight("bold")
            lbl.set_color(INK)


# ------------------------------------------------------------------
# Load data (read-only). Samples are reordered numerically (A1..A12),
# because the file stores them in text order (A1, A10, A11, ...).
# ------------------------------------------------------------------
section("Loading data")
raw = pd.read_csv(DATA_FILE, sep="\t", compression="gzip", index_col=ID_COLUMN)
samples = sorted(raw.columns, key=lambda s: int(s[1:]))
counts = raw[samples]
print(f"Matrix: {counts.shape[0]:,} transcripts x {counts.shape[1]} samples")

# ==================================================================
# Step 1: Total counts per sample (library size)
# Library size = sum of all read counts in a sample. Large differences
# mean raw counts are not directly comparable without normalisation.
# ==================================================================
section("Step 1: Library size per sample")
lib_size = counts.sum()
median_lib = lib_size.median()
lib_table = pd.DataFrame({
    "total_counts": lib_size,
    "million_reads": (lib_size / 1e6).round(2),
    "fraction_of_median": (lib_size / median_lib).round(3),
})
print(lib_table.to_string())
print(f"Median library size: {median_lib / 1e6:.1f} M; "
      f"range {lib_size.min() / 1e6:.1f}-{lib_size.max() / 1e6:.1f} M "
      f"({lib_size.max() / lib_size.min():.1f}-fold)")

# ==================================================================
# Step 2: Library size bar chart
# ==================================================================
section("Step 2: Library size bar chart")
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.bar(samples, lib_size / 1e6, color=sample_colors(samples), width=0.7,
       edgecolor="white", linewidth=1)
ax.axhline(median_lib / 1e6, color=INK_2, linestyle="--", linewidth=1)
ax.text(len(samples) - 0.5, median_lib / 1e6, f"median {median_lib / 1e6:.1f} M",
        ha="right", va="bottom", color=INK_2, fontsize=9)
a3_idx = samples.index(FOCUS_SAMPLE)
ax.annotate(f"{lib_size[FOCUS_SAMPLE] / 1e6:.1f} M\n"
            f"({lib_size[FOCUS_SAMPLE] / median_lib:.0%} of median)",
            xy=(a3_idx, lib_size[FOCUS_SAMPLE] / 1e6), xytext=(0, 6),
            textcoords="offset points", ha="center", va="bottom", fontsize=9, color=INK)
ax.set_title("Library size per sample")
ax.set_xlabel("Sample")
ax.set_ylabel("Total read counts (millions)")
ax.grid(axis="x", visible=False)
ax.legend(handles=[Patch(color=BLUE, label="Sample"),
                   Patch(color=ORANGE, label=f"{FOCUS_SAMPLE} (under investigation)")],
          loc="upper right", bbox_to_anchor=(1, 1.12), ncol=2, fontsize=9)
bold_focus_ticklabel(ax)
save_fig(fig, "01_library_sizes.png")

# ==================================================================
# Step 3: Detected transcripts per sample
# A transcript is "detected" if it has at least 1 read. Detection rises
# with sequencing depth, so a shallow library will look poorer even if it
# is technically fine. To separate depth from quality, every sample is
# randomly downsampled (binomial thinning) to A3's library size and
# detection is recounted at that equal depth.
# ==================================================================
section("Step 3: Detected transcripts per sample")
detected = (counts > 0).sum()

rng = np.random.default_rng(RANDOM_SEED)
target_depth = lib_size[FOCUS_SAMPLE]
detected_matched = {}
for s in samples:
    keep_prob = min(1.0, target_depth / lib_size[s])
    thinned = rng.binomial(counts[s].to_numpy(dtype=np.int64), keep_prob)
    detected_matched[s] = int((thinned > 0).sum())
detected_matched = pd.Series(detected_matched)

detect_table = pd.DataFrame({
    "detected_transcripts": detected,
    "pct_of_all_transcripts": (detected / len(counts) * 100).round(2),
    f"detected_at_{FOCUS_SAMPLE}_depth": detected_matched,
})
print(detect_table.to_string())

fig, ax = plt.subplots(figsize=(8.5, 4.2))
x = np.arange(len(samples))
w = 0.38
ax.bar(x - w / 2, detected / 1e3, width=w, color=BLUE, edgecolor="white",
       linewidth=1, label="All reads")
ax.bar(x + w / 2, detected_matched / 1e3, width=w, color=BLUE_LIGHT, edgecolor="white",
       linewidth=1, label=f"Downsampled to {FOCUS_SAMPLE} depth "
                          f"({target_depth / 1e6:.1f} M reads)")
ax.set_xticks(x, samples)
ax.set_title("Detected transcripts per sample (count ≥ 1)")
ax.set_xlabel("Sample")
ax.set_ylabel("Detected transcripts (thousands)")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=9)
bold_focus_ticklabel(ax)
save_fig(fig, "02_detected_transcripts.png")

# ==================================================================
# Normalisation for correlation and PCA
# CPM (counts per million) removes library-size differences; log2(CPM + 1)
# tames the extreme right skew of count data. Very low-expressed
# transcripts are filtered out because they are dominated by noise.
# ==================================================================
section("Normalisation and filtering")
cpm = counts / lib_size * 1e6
keep = (cpm >= MIN_CPM).sum(axis=1) >= MIN_SAMPLES
log_cpm = np.log2(cpm[keep] + 1)
print(f"Filter: CPM >= {MIN_CPM} in >= {MIN_SAMPLES} samples -> "
      f"{int(keep.sum()):,} of {len(counts):,} transcripts kept")

# ==================================================================
# Step 4: Sample-to-sample correlation matrix and heatmap
# Pearson correlation of log2 CPM between every pair of samples.
# Spearman (rank-based, robust to outliers) is saved as a cross-check.
# ==================================================================
section("Step 4: Sample correlation")
corr = log_cpm.corr(method="pearson")
corr_spearman = log_cpm.corr(method="spearman")
save_table(corr.round(4), "correlation_matrix_pearson_log2cpm.csv")
save_table(corr_spearman.round(4), "correlation_matrix_spearman_log2cpm.csv")

off_diag = corr.where(~np.eye(len(samples), dtype=bool))
mean_corr = off_diag.mean()
print("Mean Pearson correlation with the other samples:")
print(mean_corr.round(4).sort_values().to_string())

fig, ax = plt.subplots(figsize=(7.4, 6.2))
cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUES)
vmin = np.floor(off_diag.min().min() * 100) / 100
im = ax.imshow(corr.values, cmap=cmap, vmin=vmin, vmax=1.0)
ax.set_xticks(range(len(samples)), samples)
ax.set_yticks(range(len(samples)), samples)
ax.tick_params(length=0)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(False)
threshold = vmin + (1 - vmin) * 0.55   # switch to white text on dark cells
for i in range(len(samples)):
    for j in range(len(samples)):
        val = corr.iat[i, j]
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5,
                color="white" if val >= threshold else INK)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Pearson r (log2 CPM)", color=INK_2)
cbar.outline.set_visible(False)
ax.set_title(f"Sample-to-sample correlation ({int(keep.sum()):,} filtered transcripts)")
bold_focus_ticklabel(ax, "x")
bold_focus_ticklabel(ax, "y")
save_fig(fig, "03_correlation_heatmap.png")

# ==================================================================
# Step 5: Principal Component Analysis
# PCA on the N most variable transcripts (log2 CPM). Each transcript is
# mean-centred (sklearn does this); not scaled, so highly variable
# transcripts carry more weight, the usual choice for RNA-seq QC.
# ==================================================================
section("Step 5: PCA")
top_var = log_cpm.var(axis=1).sort_values(ascending=False).index[:N_TOP_VARIABLE]
X = log_cpm.loc[top_var].T            # rows = samples, columns = transcripts
n_pcs = min(len(samples), X.shape[1])
pca = PCA(n_components=n_pcs)
scores = pca.fit_transform(X)
pc_names = [f"PC{i + 1}" for i in range(n_pcs)]
pca_scores = pd.DataFrame(scores, index=samples, columns=pc_names)
var_ratio = pd.Series(pca.explained_variance_ratio_, index=pc_names)

variance_table = pd.DataFrame({
    "variance_explained_pct": (var_ratio * 100).round(2),
    "cumulative_pct": (var_ratio.cumsum() * 100).round(2),
})
print(variance_table.head(6).to_string())
save_table(pca_scores.round(4), "pca_sample_scores.csv")
save_table(variance_table, "pca_variance_explained.csv")

# Transcripts that drive PC1 (largest absolute loadings) - useful to see
# whether PC1 is biology or a technical effect.
loadings = pd.DataFrame(pca.components_[:2].T, index=top_var, columns=["PC1", "PC2"])
top_pc1 = loadings.reindex(loadings["PC1"].abs().sort_values(ascending=False).index).head(20)
save_table(top_pc1.round(5), "pca_top20_pc1_loadings.csv")

# ==================================================================
# Step 6: PCA scatter plot (PC1 vs PC2) and scree plot
# ==================================================================
section("Step 6: PCA plots")
fig, ax = plt.subplots(figsize=(6.8, 5.6))
ax.scatter(pca_scores["PC1"], pca_scores["PC2"], s=70, c=sample_colors(samples),
           edgecolor="white", linewidth=1.5, zorder=3)
# Default label offset is up-right; nudge labels for points that sit close
# together (checked visually) so no label overlaps another.
label_offsets = {"A7": (-8, -12), "A9": (8, -6), "A8": (8, 2)}
for s in samples:
    dx, dy = label_offsets.get(s, (6, 5))
    ax.annotate(s, (pca_scores.at[s, "PC1"], pca_scores.at[s, "PC2"]),
                xytext=(dx, dy), textcoords="offset points", fontsize=9,
                ha="right" if dx < 0 else "left",
                color=INK if s == FOCUS_SAMPLE else INK_2,
                fontweight="bold" if s == FOCUS_SAMPLE else "normal")
ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
ax.set_xlabel(f"PC1 ({var_ratio['PC1']:.1%} of variance)")
ax.set_ylabel(f"PC2 ({var_ratio['PC2']:.1%} of variance)")
ax.set_title(f"PCA of samples (top {N_TOP_VARIABLE} variable transcripts, log2 CPM)")
ax.legend(handles=[
    Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=BLUE,
           markeredgecolor="white", label="Sample"),
    Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=ORANGE,
           markeredgecolor="white", label=f"{FOCUS_SAMPLE} (under investigation)"),
], loc="best", fontsize=9)
ax.margins(0.12)
save_fig(fig, "04_pca_scatter.png")

fig, ax = plt.subplots(figsize=(7, 3.8))
shown = pc_names[:10]
ax.bar(shown, var_ratio[shown] * 100, color=BLUE, width=0.7, edgecolor="white")
ax.plot(shown, var_ratio[shown].cumsum() * 100, color=INK_2, marker="o",
        markersize=4, linewidth=1.5, label="Cumulative")
ax.set_ylim(0, 105)
ax.set_title("PCA variance explained")
ax.set_xlabel("Principal component")
ax.set_ylabel("Variance explained (%)")
ax.grid(axis="x", visible=False)
ax.legend(handles=[Patch(color=BLUE, label="Per component"),
                   Line2D([], [], color=INK_2, marker="o", markersize=4, label="Cumulative")],
          loc="center right", fontsize=9)
save_fig(fig, "05_pca_scree.png")

# ==================================================================
# Step 7: Investigating A3
# Metrics computed for every sample so A3 can be judged against peers:
#   - library size and its robust z-score (median/MAD on log10 scale)
#   - detected transcripts, raw and at equal depth
#   - read concentration: share of reads in the top 10 / top 100 transcripts
#     (high concentration = low library complexity)
#   - share of reads assigned to novel (TCONS) transcripts
#   - mean correlation with other samples
#   - distance from the PCA centroid (PC1-PC2)
# A robust z-score beyond |3.5| is a common outlier flag.
# ==================================================================
section(f"Step 7: Investigating {FOCUS_SAMPLE}")


def robust_z(series):
    med = series.median()
    mad = (series - med).abs().median()
    return 0.6745 * (series - med) / mad if mad > 0 else series * 0


sorted_props = np.sort((counts / lib_size).to_numpy(), axis=0)[::-1]
top10_share = pd.Series(sorted_props[:10].sum(axis=0), index=samples)
top100_share = pd.Series(sorted_props[:100].sum(axis=0), index=samples)
is_novel = counts.index.str.startswith("TCONS")
novel_share = counts[is_novel].sum() / lib_size
pc_dist = np.sqrt(pca_scores["PC1"] ** 2 + pca_scores["PC2"] ** 2)  # centroid is 0

qc = pd.DataFrame({
    "million_reads": (lib_size / 1e6).round(2),
    "robust_z_log_lib_size": robust_z(np.log10(lib_size)).round(2),
    "detected_transcripts": detected,
    f"detected_at_{FOCUS_SAMPLE}_depth": detected_matched,
    "robust_z_detected_matched": robust_z(detected_matched).round(2),
    "pct_reads_top10_transcripts": (top10_share * 100).round(2),
    "pct_reads_top100_transcripts": (top100_share * 100).round(2),
    "pct_reads_novel_TCONS": (novel_share * 100).round(2),
    "mean_pearson_r_to_others": mean_corr.round(4),
    "robust_z_mean_r": robust_z(mean_corr).round(2),
    "pca_distance_from_centroid": pc_dist.round(2),
})
print(qc.to_string())
save_table(qc, "sample_qc_metrics.csv")

# Which transcripts dominate A3? Compare each transcript's share of A3's
# reads to its median share in the other samples.
others = [s for s in samples if s != FOCUS_SAMPLE]
share = counts / lib_size
top_a3 = share[FOCUS_SAMPLE].sort_values(ascending=False).head(15).index
a3_top = pd.DataFrame({
    f"{FOCUS_SAMPLE}_reads": counts.loc[top_a3, FOCUS_SAMPLE],
    f"{FOCUS_SAMPLE}_pct_of_library": (share.loc[top_a3, FOCUS_SAMPLE] * 100).round(3),
    "others_median_pct_of_library": (share.loc[top_a3, others].median(axis=1) * 100).round(3),
})
a3_top["fold_vs_others_median"] = (
    a3_top[f"{FOCUS_SAMPLE}_pct_of_library"]
    / a3_top["others_median_pct_of_library"].replace(0, np.nan)
).round(1)
print(f"\nTop 15 transcripts by share of {FOCUS_SAMPLE}'s library:")
print(a3_top.to_string())
save_table(a3_top, f"{FOCUS_SAMPLE}_top15_transcripts.csv")

# Library complexity curve: cumulative share of reads captured by the top-N
# transcripts. A curve that rises faster means reads are concentrated in
# fewer transcripts.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ranks = np.arange(1, sorted_props.shape[0] + 1)
cum = np.cumsum(sorted_props, axis=0)
for j, s in enumerate(samples):
    if s != FOCUS_SAMPLE:
        ax.plot(ranks, cum[:, j] * 100, color=BLUE_LIGHT, linewidth=1.2, zorder=2)
ax.plot(ranks, cum[:, a3_idx] * 100, color=ORANGE, linewidth=2, zorder=3)
ax.set_xscale("log")
ax.set_xlim(1, 2e5)
ax.set_ylim(0, 100)
ax.set_title("Library complexity: cumulative share of reads in top-N transcripts")
ax.set_xlabel("Top N transcripts (ranked by count, log scale)")
ax.set_ylabel("Cumulative share of reads (%)")
ax.legend(handles=[Line2D([], [], color=BLUE_LIGHT, linewidth=1.5, label="Other samples"),
                   Line2D([], [], color=ORANGE, linewidth=2, label=FOCUS_SAMPLE)],
          loc="lower right", fontsize=9)
save_fig(fig, "06_library_complexity.png")

# Combined per-sample table (library size + detection) for convenience
save_table(lib_table.join(detect_table), "library_size_and_detection.csv")

# ------------------------------------------------------------------
# Plain-language verdict for A3 based on the metrics above
# ------------------------------------------------------------------
section(f"Summary for {FOCUS_SAMPLE}")
a3 = qc.loc[FOCUS_SAMPLE]
print(f"Library size: {a3['million_reads']} M reads "
      f"(robust z = {a3['robust_z_log_lib_size']}; median {median_lib / 1e6:.1f} M)")
print(f"Detected at equal depth: {int(a3[f'detected_at_{FOCUS_SAMPLE}_depth']):,} "
      f"vs. peer median {int(detected_matched[others].median()):,} "
      f"(robust z = {a3['robust_z_detected_matched']})")
print(f"Reads in top 10 transcripts: {a3['pct_reads_top10_transcripts']}% "
      f"vs. peer median {top10_share[others].median() * 100:.2f}%")
print(f"Mean correlation to others: {a3['mean_pearson_r_to_others']} "
      f"vs. peer median {mean_corr[others].median():.4f} "
      f"(robust z = {a3['robust_z_mean_r']})")
print(f"PCA: PC1 = {pca_scores.at[FOCUS_SAMPLE, 'PC1']:.2f}, "
      f"PC2 = {pca_scores.at[FOCUS_SAMPLE, 'PC2']:.2f}")
print("\nQuality control completed. A3 retained in all analyses.")
