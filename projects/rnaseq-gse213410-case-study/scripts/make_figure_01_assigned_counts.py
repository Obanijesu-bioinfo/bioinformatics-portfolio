"""
Case-study figure 01: assigned transcript counts per sample.

Redraws projects/rnaseq-gse213410-case-study/figures/01_library_sizes.png with
clearer labels. The plotted values are the same as in the Phase 3 figure
(results/figures/01_library_sizes.png, made by
data/raw/notebooks/03_rnaseq_quality_control.py): the per-sample sum of the GEO
transcript count matrix. That sum is the number of counts assigned to
transcripts, not sequencing depth, so this version says so on the axis and
shows A3's sequenced read pairs for contrast.

Inputs (read only):
  data/raw/GSE213410_transcripts.readcount.txt.gz   transcript count matrix
  results/tables/count_provenance_check.csv         SRA read pairs ("spots")
Output:
  projects/rnaseq-gse213410-case-study/figures/01_library_sizes.png

The Phase 3 script, its figure and all Phase 3-7 outputs are not modified.
The script stops if the totals it computes differ from the recorded
transcript_level_counts column.

Run from anywhere:
    python projects/rnaseq-gse213410-case-study/scripts/make_figure_01_assigned_counts.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_FILE = REPO_ROOT / "data/raw/GSE213410_transcripts.readcount.txt.gz"
PROVENANCE_FILE = REPO_ROOT / "results/tables/count_provenance_check.csv"
OUT_FILE = REPO_ROOT / "projects/rnaseq-gse213410-case-study/figures/01_library_sizes.png"
FOCUS_SAMPLE = "A3"

# Same style as the Phase 3 figure
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 10, "text.color": INK, "axes.labelcolor": INK_2,
    "axes.titlesize": 12, "axes.titleweight": "semibold", "axes.titlelocation": "left",
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2, "legend.frameon": False,
})

# Assigned counts: per-sample sum of the transcript count matrix (as in Phase 3)
raw = pd.read_csv(DATA_FILE, sep="\t", compression="gzip", index_col="transcriptID")
samples = sorted(raw.columns, key=lambda s: int(s[1:]))
assigned = raw[samples].sum()

prov = pd.read_csv(PROVENANCE_FILE, index_col="sample_id")
mismatch = (assigned - prov.loc[samples, "transcript_level_counts"]).abs().max()
if mismatch != 0:
    raise SystemExit(f"Assigned totals differ from {PROVENANCE_FILE.name} (max diff {mismatch})")
read_pairs = prov.loc[samples, "spots"]

median_assigned = assigned.median()
a3_idx = samples.index(FOCUS_SAMPLE)

fig, ax = plt.subplots(figsize=(8, 4.6))
ax.bar(samples, assigned / 1e6, color=[ORANGE if s == FOCUS_SAMPLE else BLUE for s in samples],
       width=0.7, edgecolor="white", linewidth=1)
ax.axhline(median_assigned / 1e6, color=INK_2, linestyle="--", linewidth=1)
ax.text(len(samples) - 0.5, median_assigned / 1e6, f"median {median_assigned / 1e6:.1f} M",
        ha="right", va="bottom", color=INK_2, fontsize=9)
ax.annotate(f"{assigned[FOCUS_SAMPLE] / 1e6:.1f} M assigned\n"
            f"({assigned[FOCUS_SAMPLE] / median_assigned:.0%} of median)",
            xy=(a3_idx, assigned[FOCUS_SAMPLE] / 1e6), xytext=(0, 6),
            textcoords="offset points", ha="center", va="bottom", fontsize=9, color=INK,
            bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="none"))
ax.set_title("Assigned transcript counts per sample")
ax.set_xlabel("Sample")
ax.set_ylabel("Total assigned transcript counts (millions)")
ax.grid(axis="x", visible=False)
ax.legend(handles=[Patch(color=BLUE, label="Sample"),
                   Patch(color=ORANGE, label=f"{FOCUS_SAMPLE} (QC outlier)")],
          loc="upper right", bbox_to_anchor=(1, 1.12), ncol=2, fontsize=9)
for lbl in ax.get_xticklabels():
    if lbl.get_text() == FOCUS_SAMPLE:
        lbl.set_fontweight("bold")
        lbl.set_color(INK)
fig.text(0.0, -0.04,
         "Bars: sum of the GEO transcript count matrix per sample (counts assigned to transcripts), "
         f"not sequencing depth.\nSequenced read pairs (SRA): {FOCUS_SAMPLE} "
         f"{read_pairs[FOCUS_SAMPLE] / 1e6:.1f} M; all samples "
         f"{read_pairs.min() / 1e6:.1f}–{read_pairs.max() / 1e6:.1f} M. "
         f"{FOCUS_SAMPLE}'s low total reflects low assignment, not shallow sequencing.",
         ha="left", va="top", fontsize=8, color=MUTED)

fig.savefig(OUT_FILE, facecolor="white")
plt.close(fig)
print(f"saved {OUT_FILE.relative_to(REPO_ROOT)}")
print(f"A3: {assigned[FOCUS_SAMPLE]:,} assigned counts; {read_pairs[FOCUS_SAMPLE]:,} read pairs")
