"""
Phase 5: Transcript-to-gene mapping, gene-level aggregation and gene-level QC
Dataset: NCBI GEO GSE213410 (GSE213410_transcripts.readcount.txt.gz)

Annotation
    Ensembl release 97, GRCh38.p12 (data/reference/Homo_sapiens.GRCh38.97.gtf.gz)
    Evidence: GEO states GRCh38/hg38, and the count file's 226,788 ENST IDs are
    exactly the release-97 transcript set (see data/reference/ANNOTATION_SOURCE.md).
    This script checks that exact match again on every run.

Steps
    0. SHA-256 of the original count file (checked again at the end)
    1. Verify and parse the Ensembl 97 GTF (transcript -> gene)
    2. Load counts, check samples against the metadata, map ENST -> ENSG
    3. Summarise TCONS transcripts (not annotated; see report)
    4. Sum raw integer counts per gene (mapped ENST only)
    5. Mapping statistics
    6. Gene-level QC (library size, detection, distribution, correlation, PCA)
    7. A3 sensitivity check (C1 with vs. without A3; A3 is not removed)
    8. Validation report and PHASE5_REPORT.md

Exploratory normalisation (QC only)
    - Size factors: DESeq2-style median-of-ratios (genes with counts > 0 in
      every sample of the set being analysed).
    - Transform: log2(normalised count + 1).
    - Expression filter: normalised count >= 10 in >= 3 samples (3 = group size).
    - PCA: top 500 most variable filtered genes, centred, unscaled.
    - PCA and the sensitivity check are done separately per comparison
      (C1, C2). The two comparisons are never pooled for those analyses.
    - Correlation heatmap and distribution plots show all 12 samples side by
      side for descriptive QC only (global size factors).
    - Transcript- vs gene-level comparison uses the Phase 3 method exactly
      (log2 CPM + 1, CPM >= 1 in >= 3 samples) so the two are comparable.

No differential expression is performed. The original count file is read-only.

Run from the repository root:
    python data/raw/notebooks/05_transcript_to_gene_aggregation.py
"""

import gzip
import hashlib
import platform
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from sklearn.decomposition import PCA

# ------------------------------------------------------------------
# Paths, annotation identity and parameters
# ------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[3]
COUNTS_FILE = REPO_ROOT / "data/raw/GSE213410_transcripts.readcount.txt.gz"
METADATA_FILE = REPO_ROOT / "data/metadata/sample_metadata.csv"
GTF_FILE = REPO_ROOT / "data/reference/Homo_sapiens.GRCh38.97.gtf.gz"
GTF_URL = ("https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/"
           "Homo_sapiens.GRCh38.97.gtf.gz")
GTF_SHA256 = "8d2b55f734f4e351b37d24eddba7d955034c98390564260318341ee464590a19"
ANNOTATION_LABEL = "Ensembl release 97 (GRCh38.p12)"
FIG_DIR = REPO_ROOT / "results/figures"
TABLE_DIR = REPO_ROOT / "results/tables"
REPORT_FILE = REPO_ROOT / "results/PHASE5_REPORT.md"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_SAMPLES = [f"A{i}" for i in range(1, 13)]
OUTLIER_SAMPLE = "A3"
MIN_NORM_COUNT = 10        # expression filter for QC: normalised count >= 10 ...
MIN_SAMPLES = 3            # ... in >= 3 samples (smallest group size)
N_TOP_VARIABLE = 500       # genes used for PCA
CPM_MIN = 1.0              # Phase 3 filter, reused for the like-for-like comparison
RANDOM_SEED = 42           # downsampling reproducibility

# ------------------------------------------------------------------
# Plot style (same scientific style as Phase 3). Colours identify the
# experimental group; A3 is marked with a dark ring / hatch, never by colour.
# Palette validated for colour-vision deficiency (adjacent pairs).
# ------------------------------------------------------------------
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, EXCLUDED_GREY = "#e1e0d9", "#c3c2b7", "#c3c2b7"
GROUP_COLORS = {
    "hESC_Ethanol": "#2a78d6",
    "hESC_5-OXO-ETE": "#eb6834",
    "T21-hiPSC_DMSO": "#1baf7a",
    "T21-hiPSC_Calcitriol": "#eda100",
}
GROUP_LABELS = {
    "hESC_Ethanol": "C1 hESC: Ethanol (control)",
    "hESC_5-OXO-ETE": "C1 hESC: 5-OXO-ETE",
    "T21-hiPSC_DMSO": "C2 T21 hiPSC: DMSO (control)",
    "T21-hiPSC_Calcitriol": "C2 T21 hiPSC: Calcitriol",
}
SEQ_BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

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

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
validation_checks = []


def section(title):
    print("\n" + "=" * 68)
    print(title)
    print("=" * 68)


def record_check(name, passed, detail=""):
    """Record a validation check; failures are collected and stop the run at the end."""
    status = "PASS" if passed else "FAIL"
    validation_checks.append({"check": name, "result": status, "detail": str(detail)})
    print(f"  [{status}] {name}" + (f" - {detail}" if detail != "" else ""))
    return passed


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_table(df, name, index=False):
    path = TABLE_DIR / name
    df.to_csv(path, index=index)
    print(f"  saved table:  {path.relative_to(REPO_ROOT)}")


def save_fig(fig, name):
    path = FIG_DIR / name
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    print(f"  saved figure: {path.relative_to(REPO_ROOT)}")


def strip_version(identifier):
    """Remove an Ensembl version suffix: ENST00000456328.2 -> ENST00000456328."""
    return re.sub(r"\.\d+$", "", identifier)


def median_of_ratios_size_factors(count_df):
    """DESeq2-style size factors from genes with counts > 0 in every sample."""
    positive = count_df[(count_df > 0).all(axis=1)]
    if positive.empty:
        raise ValueError("No genes are expressed in all samples; size factors undefined.")
    log_counts = np.log(positive.astype(float))
    log_ratios = log_counts.sub(log_counts.mean(axis=1), axis=0)
    return np.exp(log_ratios.median())


def log_normalise_and_filter(count_df):
    """Size-factor normalise, filter low expression, return (log2 matrix, size factors)."""
    size_factors = median_of_ratios_size_factors(count_df)
    normalised = count_df / size_factors
    keep = (normalised >= MIN_NORM_COUNT).sum(axis=1) >= min(MIN_SAMPLES, count_df.shape[1])
    return np.log2(normalised[keep] + 1), size_factors


def run_pca(log_matrix):
    """PCA on the most variable genes. Returns scores, variance ratios and gene count used."""
    top_genes = log_matrix.var(axis=1).sort_values(ascending=False).index[:N_TOP_VARIABLE]
    sample_by_gene = log_matrix.loc[top_genes].T
    n_components = min(sample_by_gene.shape)
    pca = PCA(n_components=n_components)
    scores = pd.DataFrame(pca.fit_transform(sample_by_gene), index=sample_by_gene.index,
                          columns=[f"PC{i + 1}" for i in range(n_components)])
    variance = pd.Series(pca.explained_variance_ratio_, index=scores.columns)
    return scores, variance, len(top_genes)


def mean_offdiagonal_r(corr, members=None):
    """Mean correlation of each sample with the other samples in `members`."""
    members = list(corr.index) if members is None else members
    sub = corr.loc[members, members]
    return sub.where(~np.eye(len(members), dtype=bool)).mean()


def pairwise_mean(corr, set_a, set_b):
    values = [corr.at[a, b] for a in set_a for b in set_b if a != b]
    return float(np.mean(values)) if values else np.nan


def robust_z(series):
    median = series.median()
    mad = (series - median).abs().median()
    return 0.6745 * (series - median) / mad if mad > 0 else series * 0.0


def log2_cpm_corr(count_df):
    """Phase 3 method: Pearson r of log2(CPM + 1), CPM >= 1 in >= 3 samples."""
    cpm = count_df / count_df.sum() * 1e6
    keep = (cpm >= CPM_MIN).sum(axis=1) >= MIN_SAMPLES
    return np.log2(cpm[keep] + 1).corr(), int(keep.sum())


# ==================================================================
# Step 0: Checksum of the original dataset (before)
# ==================================================================
section("Step 0: Original dataset checksum")
if not COUNTS_FILE.exists():
    sys.exit(f"ERROR: count file not found: {COUNTS_FILE}")
counts_sha256_before = sha256_of(COUNTS_FILE)
print(f"SHA-256 (before): {counts_sha256_before}")

# ==================================================================
# Step 1: Annotation - download if absent, verify checksum, parse
# ==================================================================
section("Step 1: Ensembl 97 annotation")
if not GTF_FILE.exists():
    print(f"Annotation not found locally; downloading {GTF_URL}")
    GTF_FILE.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(GTF_URL, GTF_FILE)
gtf_sha256 = sha256_of(GTF_FILE)
if not record_check("Annotation GTF SHA-256 matches recorded value", gtf_sha256 == GTF_SHA256,
                    gtf_sha256):
    sys.exit("ERROR: annotation file checksum mismatch - refusing to continue.")

genome_build = None
transcript_rows = []
attribute_pattern = re.compile(r'(\S+) "([^"]*)"')
with gzip.open(GTF_FILE, "rt", encoding="utf-8") as gtf:
    for line in gtf:
        if line.startswith("#!genome-build "):
            genome_build = line.split(" ", 1)[1].strip()
            continue
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        attributes = dict(attribute_pattern.findall(fields[8]))
        transcript_rows.append({
            "transcript_id": attributes.get("transcript_id"),
            "transcript_version": attributes.get("transcript_version", ""),
            "gene_id": attributes.get("gene_id"),
            "gene_version": attributes.get("gene_version", ""),
            "gene_name": attributes.get("gene_name", ""),
            "gene_biotype": attributes.get("gene_biotype", ""),
            "transcript_biotype": attributes.get("transcript_biotype", ""),
            "chromosome": fields[0],
        })
gtf_transcripts = pd.DataFrame(transcript_rows)
print(f"Genome build in GTF header: {genome_build}")
print(f"Transcript records in GTF: {len(gtf_transcripts):,}")
record_check("GTF genome build is GRCh38", bool(genome_build) and genome_build.startswith("GRCh38"),
             genome_build)
record_check("Every GTF transcript has a transcript_id and gene_id",
             gtf_transcripts[["transcript_id", "gene_id"]].notna().all().all())

# Ambiguity: a transcript ID (version-stripped) linked to more than one gene.
gtf_transcripts["transcript_id_unversioned"] = gtf_transcripts["transcript_id"].map(strip_version)
genes_per_transcript = gtf_transcripts.groupby("transcript_id_unversioned")["gene_id"].nunique()
ambiguous_in_annotation = set(genes_per_transcript[genes_per_transcript > 1].index)
print(f"Transcript IDs linked to >1 gene in annotation: {len(ambiguous_in_annotation)}")
annotation_lookup = (gtf_transcripts
                     .drop_duplicates("transcript_id_unversioned")
                     .set_index("transcript_id_unversioned"))

# ==================================================================
# Step 2: Load counts, check samples, map ENST -> ENSG
# ==================================================================
section("Step 2: Load counts and map transcripts to genes")
raw_counts = pd.read_csv(COUNTS_FILE, sep="\t", compression="gzip", index_col="transcriptID")
metadata = pd.read_csv(METADATA_FILE)
print(f"Count matrix: {raw_counts.shape[0]:,} transcripts x {raw_counts.shape[1]} samples")

record_check("All 12 samples present in count file",
             sorted(raw_counts.columns) == sorted(EXPECTED_SAMPLES), list(raw_counts.columns))
record_check("Count-file samples match sample_metadata.csv",
             set(raw_counts.columns) == set(metadata["sample_id"]),
             f"metadata: {metadata['sample_id'].tolist()}")
record_check("Transcript IDs unique in count file", raw_counts.index.is_unique)
record_check("Original counts are non-negative integers",
             all(pd.api.types.is_integer_dtype(t) for t in raw_counts.dtypes)
             and bool((raw_counts >= 0).all().all()))

transcript_counts = raw_counts[EXPECTED_SAMPLES]      # enforce A1..A12 order
metadata = metadata.set_index("sample_id").loc[EXPECTED_SAMPLES]

original_ids = transcript_counts.index.to_series()
id_class = np.select(
    [original_ids.str.startswith("ENST"), original_ids.str.startswith("TCONS")],
    ["ENST", "TCONS"], default="other")
id_class = pd.Series(id_class, index=original_ids.index)
unversioned_ids = original_ids.map(strip_version)
n_with_version_suffix = int((unversioned_ids != original_ids).sum())
collisions = unversioned_ids[id_class == "ENST"].duplicated(keep=False)
record_check("Stripping version suffixes creates no ID collisions", not collisions.any(),
             f"{n_with_version_suffix} IDs carried a version suffix")

# Release-match check: file ENST set must equal the release-97 transcript set.
file_enst = set(unversioned_ids[id_class == "ENST"])
annotation_enst = set(annotation_lookup.index)
only_in_file = file_enst - annotation_enst
only_in_annotation = annotation_enst - file_enst
record_check("Count-file ENST IDs exactly match Ensembl 97 transcript set",
             not only_in_file and not only_in_annotation,
             f"only in file: {len(only_in_file)}, only in annotation: {len(only_in_annotation)}")

# Assign a mapping status to every transcript. Nothing is mapped by guesswork:
# a transcript is 'mapped' only if its ID is in the annotation with one gene.
mapping = pd.DataFrame({
    "original_transcript_id": original_ids.values,
    "transcript_id_unversioned": unversioned_ids.values,
    "id_type": id_class.values,
}, index=original_ids.index)
in_annotation = mapping["transcript_id_unversioned"].isin(annotation_lookup.index)
is_ambiguous = mapping["transcript_id_unversioned"].isin(ambiguous_in_annotation)
mapping["mapping_status"] = np.select(
    [
        (mapping["id_type"] == "ENST") & in_annotation & ~is_ambiguous,
        (mapping["id_type"] == "ENST") & is_ambiguous,
        (mapping["id_type"] == "ENST") & ~in_annotation,
        mapping["id_type"] == "TCONS",
    ],
    ["mapped", "ambiguous_multiple_genes", "unmapped_not_in_annotation",
     "unmapped_novel_TCONS_no_assembly_annotation"],
    default="unmapped_unrecognised_prefix")
annotation_columns = ["gene_id", "gene_name", "gene_biotype", "transcript_biotype",
                      "transcript_version", "chromosome"]
mapping = mapping.join(annotation_lookup[annotation_columns], on="transcript_id_unversioned")
mapping.loc[mapping["mapping_status"] != "mapped", annotation_columns] = np.nan
mapping["total_count_all_samples"] = transcript_counts.sum(axis=1)
mapping["n_samples_detected"] = (transcript_counts > 0).sum(axis=1)
print(mapping["mapping_status"].value_counts().to_string())
save_table(mapping, "transcript_gene_mapping.csv")

# ==================================================================
# Step 3: TCONS transcripts (novel StringTie assemblies)
# GEO provides no merged assembly GTF or sequences, so these cannot be
# placed against Ensembl genes. They are reported, not mapped.
# ==================================================================
section("Step 3: TCONS transcripts")
is_tcons = mapping["id_type"] == "TCONS"
tcons_counts = transcript_counts[is_tcons.values]
library_totals = transcript_counts.sum()
tcons_share = tcons_counts.sum() / library_totals * 100
print(f"TCONS transcripts: {int(is_tcons.sum()):,}; "
      f"detected in >=1 sample: {int((tcons_counts.sum(axis=1) > 0).sum()):,}; "
      f"zero in all samples: {int((tcons_counts.sum(axis=1) == 0).sum()):,}")
print(f"Share of counts from TCONS: {tcons_counts.values.sum() / library_totals.sum():.2%} overall")
print(tcons_share.round(2).rename("pct_counts_TCONS").to_string())

unmapped = mapping[mapping["mapping_status"] != "mapped"][
    ["original_transcript_id", "id_type", "mapping_status",
     "total_count_all_samples", "n_samples_detected"]].copy()
unmapped["reason"] = unmapped["mapping_status"].map({
    "unmapped_novel_TCONS_no_assembly_annotation":
        "Novel StringTie transcript; no assembly GTF/sequences published on GEO",
    "unmapped_not_in_annotation": "ENST ID absent from Ensembl 97",
    "ambiguous_multiple_genes": "Transcript ID linked to more than one gene",
    "unmapped_unrecognised_prefix": "Unrecognised identifier prefix",
})
unmapped = unmapped.join(transcript_counts)
save_table(unmapped, "unmapped_transcripts.csv")

# ==================================================================
# Step 4: Gene-level counts (sum of raw integer counts per gene)
# ==================================================================
section("Step 4: Aggregate to gene level")
mapped_mask = (mapping["mapping_status"] == "mapped").values
mapped_counts = transcript_counts[mapped_mask]
gene_ids_for_mapped = mapping.loc[mapped_mask, "gene_id"].values
gene_counts = mapped_counts.groupby(gene_ids_for_mapped).sum().sort_index()
gene_counts.index.name = "gene_id"
gene_counts = gene_counts[EXPECTED_SAMPLES].astype(np.int64)
print(f"Gene-level matrix: {gene_counts.shape[0]:,} genes x {gene_counts.shape[1]} samples")
save_table(gene_counts, "gene_level_counts.csv", index=True)

gene_annotation = (mapping[mapped_mask]
                   .groupby("gene_id")
                   .agg(gene_name=("gene_name", "first"), gene_biotype=("gene_biotype", "first"),
                        chromosome=("chromosome", "first"),
                        n_transcripts=("original_transcript_id", "size"))
                   .loc[gene_counts.index])
save_table(gene_annotation, "gene_annotation.csv", index=True)

record_check("Gene-level counts are non-negative integers",
             all(pd.api.types.is_integer_dtype(t) for t in gene_counts.dtypes)
             and bool((gene_counts >= 0).all().all()))
record_check("No duplicate gene IDs", gene_counts.index.is_unique)
record_check("No missing values in gene-level matrix", not gene_counts.isna().any().any())
record_check("Gene matrix columns are A1..A12 in order",
             list(gene_counts.columns) == EXPECTED_SAMPLES)
record_check("Gene IDs are Ensembl ENSG identifiers",
             bool(gene_counts.index.str.fullmatch(r"ENSG\d{11}").all()))

# ==================================================================
# Step 5: Mapping statistics and reconciliation
# ==================================================================
section("Step 5: Mapping statistics")
status_counts = mapping["mapping_status"].value_counts()
n_total = len(mapping)
n_enst = int((mapping["id_type"] == "ENST").sum())
n_tcons = int(is_tcons.sum())
n_other = int((mapping["id_type"] == "other").sum())
n_mapped = int(status_counts.get("mapped", 0))
n_ambiguous = int(status_counts.get("ambiguous_multiple_genes", 0))
n_enst_unmapped = int(status_counts.get("unmapped_not_in_annotation", 0))
total_original_counts = int(transcript_counts.values.sum())
total_gene_counts = int(gene_counts.values.sum())
total_tcons_counts = int(tcons_counts.values.sum())
pct_counts_in_genes = total_gene_counts / total_original_counts * 100

summary_rows = [
    ("annotation", ANNOTATION_LABEL),
    ("annotation_file", GTF_FILE.name),
    ("annotation_sha256", gtf_sha256),
    ("genome_build", genome_build),
    ("total_transcripts", n_total),
    ("ENST_transcripts", n_enst),
    ("TCONS_transcripts", n_tcons),
    ("other_prefix_transcripts", n_other),
    ("ENST_with_version_suffix", n_with_version_suffix),
    ("ENST_mapped_to_one_gene", n_mapped),
    ("ENST_unmapped", n_enst_unmapped),
    ("ENST_ambiguous_multiple_genes", n_ambiguous),
    ("TCONS_annotated", 0),
    ("TCONS_unannotated", n_tcons),
    ("unique_genes_in_matrix", gene_counts.shape[0]),
    ("genes_with_nonzero_counts", int((gene_counts.sum(axis=1) > 0).sum())),
    ("total_counts_original", total_original_counts),
    ("total_counts_in_gene_matrix", total_gene_counts),
    ("total_counts_TCONS_excluded", total_tcons_counts),
    ("total_counts_unmapped_or_ambiguous_ENST_excluded",
     total_original_counts - total_gene_counts - total_tcons_counts),
    ("pct_counts_in_gene_matrix", round(pct_counts_in_genes, 3)),
    ("pct_counts_TCONS", round(total_tcons_counts / total_original_counts * 100, 3)),
]
mapping_summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
print(mapping_summary.to_string(index=False))
save_table(mapping_summary, "transcript_mapping_summary.csv")

# Reconciliation: every count and every transcript must be accounted for.
record_check("ENST: mapped + unmapped + ambiguous = total ENST",
             n_mapped + n_enst_unmapped + n_ambiguous == n_enst,
             f"{n_mapped} + {n_enst_unmapped} + {n_ambiguous} = {n_enst}")
record_check("ENST + TCONS + other = total transcripts",
             n_enst + n_tcons + n_other == n_total, f"{n_enst} + {n_tcons} + {n_other} = {n_total}")
record_check("Mapping table has one row per original transcript",
             len(mapping) == len(transcript_counts) and mapping["original_transcript_id"].is_unique)
record_check("Gene count = unique gene IDs among mapped transcripts",
             gene_counts.shape[0] == len(set(gene_ids_for_mapped)))
per_sample_gene_totals = gene_counts.sum()
record_check("Per-sample gene totals equal mapped-transcript totals",
             per_sample_gene_totals.equals(mapped_counts.sum()))
excluded_totals = transcript_counts[~mapped_mask].sum()
record_check("Per-sample gene totals + excluded totals = original totals",
             (per_sample_gene_totals + excluded_totals).equals(library_totals))

# ==================================================================
# Step 6: Gene-level QC
# ==================================================================
section("Step 6: Gene-level QC")
sample_groups = metadata["experimental_group"]
sample_comparison = metadata["comparison"]
comparison_samples = {c: list(s.index) for c, s in sample_comparison.groupby(sample_comparison, sort=False)}
C1, C2 = list(comparison_samples)          # C1 = hESC, C2 = T21 hiPSC (metadata order)
colors = [GROUP_COLORS[sample_groups[s]] for s in EXPECTED_SAMPLES]
outlier_index = EXPECTED_SAMPLES.index(OUTLIER_SAMPLE)


def group_legend_handles(groups, kind="patch"):
    handles = []
    for g in groups:
        if kind == "patch":
            handles.append(Patch(color=GROUP_COLORS[g], label=GROUP_LABELS[g]))
        else:
            handles.append(Line2D([], [], marker="o", linestyle="", markersize=8,
                                  markerfacecolor=GROUP_COLORS[g], markeredgecolor="white",
                                  label=GROUP_LABELS[g]))
    return handles


def emphasise_outlier_ticklabels(ax, axes=("x",)):
    for axis in axes:
        labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
        for label in labels:
            if label.get_text() == OUTLIER_SAMPLE:
                label.set_fontweight("bold")
                label.set_color(INK)


# --- Library size and detection ---------------------------------------------
gene_library_size = gene_counts.sum()
detected_genes = (gene_counts > 0).sum()
rng = np.random.default_rng(RANDOM_SEED)
target_depth = gene_library_size[OUTLIER_SAMPLE]
detected_at_outlier_depth = {}
for sample in EXPECTED_SAMPLES:
    keep_probability = min(1.0, target_depth / gene_library_size[sample])
    thinned = rng.binomial(gene_counts[sample].to_numpy(), keep_probability)
    detected_at_outlier_depth[sample] = int((thinned > 0).sum())
detected_at_outlier_depth = pd.Series(detected_at_outlier_depth)
sorted_gene_props = np.sort((gene_counts / gene_library_size).to_numpy(), axis=0)[::-1]
top10_gene_share = pd.Series(sorted_gene_props[:10].sum(axis=0), index=EXPECTED_SAMPLES)

# --- Figure 07: gene-level library size (stacked with excluded counts) -------
fig, ax = plt.subplots(figsize=(8.8, 4.6))
x = np.arange(len(EXPECTED_SAMPLES))
ax.bar(x, gene_library_size / 1e6, color=colors, width=0.7, edgecolor="white", linewidth=1)
ax.bar(x, excluded_totals / 1e6, bottom=gene_library_size / 1e6, color=EXCLUDED_GREY,
       width=0.7, edgecolor="white", linewidth=1)
ax.patches[outlier_index].set_hatch("///")
ax.patches[outlier_index].set_edgecolor(INK)
median_gene_lib = gene_library_size.median()
ax.axhline(median_gene_lib / 1e6, color=INK_2, linestyle="--", linewidth=1)
ax.annotate(f"{gene_library_size[OUTLIER_SAMPLE] / 1e6:.1f} M\n"
            f"({gene_library_size[OUTLIER_SAMPLE] / median_gene_lib:.0%} of median)",
            xy=(outlier_index, library_totals[OUTLIER_SAMPLE] / 1e6), xytext=(0, 6),
            textcoords="offset points", ha="center", va="bottom", fontsize=9)
ax.set_xticks(x, EXPECTED_SAMPLES)
ax.set_title("Gene-level library size (Ensembl 97 genes) and excluded transcript counts")
ax.set_xlabel("Sample")
ax.set_ylabel("Read counts (millions)")
ax.grid(axis="x", visible=False)
ax.legend(handles=group_legend_handles(GROUP_COLORS) + [
    Patch(color=EXCLUDED_GREY, label="Excluded (TCONS, not in gene matrix)"),
    Patch(facecolor="white", edgecolor=INK, hatch="///", label=f"{OUTLIER_SAMPLE}: QC outlier"),
    Line2D([], [], color=INK_2, linestyle="--", linewidth=1,
           label=f"Gene-level median ({median_gene_lib / 1e6:.1f} M)")],
    loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=8.5)
emphasise_outlier_ticklabels(ax)
save_fig(fig, "07_gene_library_sizes.png")

# --- Figure 08: detected genes ----------------------------------------------
fig, ax = plt.subplots(figsize=(8.8, 4.6))
bars = ax.bar(x, detected_genes / 1e3, color=colors, width=0.7, edgecolor="white", linewidth=1)
bars[outlier_index].set_hatch("///")
bars[outlier_index].set_edgecolor(INK)
ax.scatter(x, detected_at_outlier_depth / 1e3, marker="_", s=420, linewidths=2.2,
           color=INK, zorder=3)
ax.set_xticks(x, EXPECTED_SAMPLES)
ax.set_title("Detected genes per sample (count ≥ 1)")
ax.set_xlabel("Sample")
ax.set_ylabel("Detected genes (thousands)")
ax.grid(axis="x", visible=False)
ax.legend(handles=group_legend_handles(GROUP_COLORS) + [
    Line2D([], [], marker="_", linestyle="", markersize=16, markeredgewidth=2.2, color=INK,
           label=f"Detected after downsampling to {OUTLIER_SAMPLE} depth "
                 f"({target_depth / 1e6:.1f} M)"),
    Patch(facecolor="white", edgecolor=INK, hatch="///", label=f"{OUTLIER_SAMPLE}: QC outlier")],
    loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=8.5)
emphasise_outlier_ticklabels(ax)
save_fig(fig, "08_detected_genes.png")

# --- Global normalisation for descriptive all-sample plots ------------------
global_log, global_size_factors = log_normalise_and_filter(gene_counts)
print(f"Global filter (norm >= {MIN_NORM_COUNT} in >= {MIN_SAMPLES}): "
      f"{len(global_log):,} genes kept")

# --- Figure 09: expression distribution -------------------------------------
fig, (ax_box, ax_density) = plt.subplots(1, 2, figsize=(12, 4.6),
                                         gridspec_kw={"width_ratios": [1.25, 1]})
box = ax_box.boxplot([global_log[s] for s in EXPECTED_SAMPLES], patch_artist=True,
                     showfliers=False, widths=0.6,
                     medianprops={"color": INK, "linewidth": 1.4},
                     whiskerprops={"color": MUTED}, capprops={"color": MUTED})
for patch, sample, color in zip(box["boxes"], EXPECTED_SAMPLES, colors):
    patch.set_facecolor(color)
    patch.set_edgecolor(INK if sample == OUTLIER_SAMPLE else "white")
    if sample == OUTLIER_SAMPLE:
        patch.set_hatch("///")
ax_box.set_xticks(range(1, len(EXPECTED_SAMPLES) + 1), EXPECTED_SAMPLES)
ax_box.set_title(f"a  Filtered genes (n = {len(global_log):,})")
ax_box.set_xlabel("Sample")
ax_box.set_ylabel("log2(normalised count + 1)")
ax_box.grid(axis="x", visible=False)
emphasise_outlier_ticklabels(ax_box)

# Empirical cumulative distribution (ECDF): avoids the spiky artefacts that
# histogram densities show for low, discrete counts.
all_detected_log = np.log2(gene_counts / global_size_factors + 1)
for sample, color in zip(EXPECTED_SAMPLES, colors):
    values = np.sort(all_detected_log.loc[gene_counts[sample] > 0, sample].to_numpy())
    cumulative = np.arange(1, len(values) + 1) / len(values)
    is_outlier = sample == OUTLIER_SAMPLE
    ax_density.plot(values, cumulative, color=INK if is_outlier else color,
                    linewidth=2 if is_outlier else 1.2, linestyle="--" if is_outlier else "-",
                    zorder=3 if is_outlier else 2)
ax_density.set_ylim(0, 1.02)
ax_density.set_title("b  All detected genes (count ≥ 1): cumulative distribution")
ax_density.set_xlabel("log2(normalised count + 1)")
ax_density.set_ylabel("Cumulative fraction of detected genes")
ax_density.legend(handles=[Line2D([], [], color=GROUP_COLORS[g], linewidth=1.5, label=GROUP_LABELS[g])
                           for g in GROUP_COLORS]
                  + [Line2D([], [], color=INK, linewidth=2, linestyle="--", label=OUTLIER_SAMPLE)],
                  loc="lower right", fontsize=8.5)
fig.suptitle("Gene-level expression distribution (median-of-ratios normalised)",
             x=0.01, ha="left", fontsize=12, fontweight="semibold")
fig.tight_layout()
save_fig(fig, "09_gene_expression_distribution.png")

# --- Correlation (descriptive, all 12) + Figure 10 --------------------------
gene_corr = global_log.corr(method="pearson")
save_table(gene_corr.round(4), "gene_level_correlation.csv", index=True)

fig, ax = plt.subplots(figsize=(8, 7))
offdiag = gene_corr.where(~np.eye(12, dtype=bool))
vmin = np.floor(offdiag.min().min() * 100) / 100
image = ax.imshow(gene_corr.values, cmap=LinearSegmentedColormap.from_list("blues", SEQ_BLUES),
                  vmin=vmin, vmax=1.0)
threshold = vmin + (1 - vmin) * 0.55
for i in range(12):
    for j in range(12):
        value = gene_corr.iat[i, j]
        ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7.5,
                color="white" if value >= threshold else INK)
ax.set_xticks(range(12), EXPECTED_SAMPLES)
ax.set_yticks(range(12), EXPECTED_SAMPLES)
ax.tick_params(length=0)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.axhline(5.5, color="white", linewidth=3)
ax.axvline(5.5, color="white", linewidth=3)
for j, sample in enumerate(EXPECTED_SAMPLES):   # group colour strip above the matrix
    ax.add_patch(Rectangle((j - 0.5, -1.05), 1, 0.35, color=GROUP_COLORS[sample_groups[sample]],
                           clip_on=False))
ax.set_ylim(11.5, -1.15)
colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
colorbar.set_label("Pearson r, log2(normalised count + 1)", color=INK_2)
colorbar.outline.set_visible(False)
ax.set_title(f"Gene-level sample correlation ({len(global_log):,} filtered genes)", pad=12)
ax.legend(handles=group_legend_handles(GROUP_COLORS), loc="upper center",
          bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=8.5)
emphasise_outlier_ticklabels(ax, ("x", "y"))
save_fig(fig, "10_gene_correlation_heatmap.png")

# --- Per-comparison normalisation and PCA -----------------------------------
per_comparison = {}
for comparison, members in comparison_samples.items():
    log_matrix, size_factors = log_normalise_and_filter(gene_counts[members])
    scores, variance, n_pca_genes = run_pca(log_matrix)
    per_comparison[comparison] = {
        "members": members, "log": log_matrix, "size_factors": size_factors,
        "scores": scores, "variance": variance, "corr": log_matrix.corr(),
    }
    print(f"{comparison}: {len(log_matrix):,} genes after filter; "
          f"PC1 {variance['PC1']:.1%}, PC2 {variance['PC2']:.1%}")


def draw_pca(ax, scores, variance, title, label_offsets=None):
    label_offsets = label_offsets or {}
    for sample in scores.index:
        is_outlier = sample == OUTLIER_SAMPLE
        ax.scatter(scores.at[sample, "PC1"], scores.at[sample, "PC2"], s=80,
                   color=GROUP_COLORS[sample_groups[sample]],
                   edgecolor=INK if is_outlier else "white",
                   linewidth=2 if is_outlier else 1.5, zorder=3)
        dx, dy = label_offsets.get(sample, (7, 5))
        ax.annotate(sample, (scores.at[sample, "PC1"], scores.at[sample, "PC2"]),
                    xytext=(dx, dy), textcoords="offset points", fontsize=9,
                    ha="right" if dx < 0 else "left",
                    color=INK if is_outlier else INK_2,
                    fontweight="bold" if is_outlier else "normal")
    ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
    ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
    ax.set_xlabel(f"PC1 ({variance['PC1']:.1%} of variance)")
    ax.set_ylabel(f"PC2 ({variance['PC2']:.1%} of variance)")
    ax.set_title(title)
    ax.margins(0.2)
    groups_present = list(dict.fromkeys(sample_groups[scores.index]))
    handles = group_legend_handles(groups_present, kind="marker")
    if OUTLIER_SAMPLE in scores.index:
        handles.append(Line2D([], [], marker="o", linestyle="", markersize=8,
                              markerfacecolor="white", markeredgecolor=INK, markeredgewidth=2,
                              label=f"{OUTLIER_SAMPLE}: QC outlier"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=1,
              fontsize=8.5)


# --- Figure 11: gene-level PCA, each comparison separately ------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
for ax, (comparison, panel) in zip(axes, per_comparison.items()):
    letter = "a" if comparison == C1 else "b"
    draw_pca(ax, panel["scores"], panel["variance"],
             f"{letter}  {comparison.split('_')[0]}: "
             f"{'H9 hESC' if comparison == C1 else 'T21 hiPSC'} (n = {len(panel['members'])})")
fig.suptitle(f"Gene-level PCA per comparison (top {N_TOP_VARIABLE} variable genes, "
             "log2 median-of-ratios normalised)", x=0.01, ha="left",
             fontsize=12, fontweight="semibold")
fig.tight_layout()
save_fig(fig, "11_gene_level_pca.png")

# ==================================================================
# Step 7: A3 sensitivity check (C1 with vs. without A3)
# A3 is removed only inside this temporary comparison; the primary
# gene-level matrix keeps all 12 samples.
# ==================================================================
section(f"Step 7: {OUTLIER_SAMPLE} sensitivity (C1)")
c1_members = comparison_samples[C1]
c1_without = [s for s in c1_members if s != OUTLIER_SAMPLE]
c1_log_wo, c1_sf_wo = log_normalise_and_filter(gene_counts[c1_without])
c1_scores_wo, c1_var_wo, _ = run_pca(c1_log_wo)
c1_corr_wo = c1_log_wo.corr()


def design_metrics(label, members, log_matrix, corr, scores, variance, size_factors):
    """Descriptive structure metrics for one sample set (no statistical testing)."""
    groups = sample_groups[members]
    group_names = list(dict.fromkeys(groups))
    first = [s for s in members if groups[s] == group_names[0]]
    second = [s for s in members if groups[s] == group_names[1]]
    pc1 = scores["PC1"]
    separates = (pc1[first].max() < pc1[second].min()) or (pc1[second].max() < pc1[first].min())
    return {
        "analysis": label,
        "n_samples": len(members),
        "samples": ", ".join(members),
        "genes_after_filter": len(log_matrix),
        "PC1_pct_variance": round(variance["PC1"] * 100, 2),
        "PC2_pct_variance": round(variance["PC2"] * 100, 2),
        "PC1_separates_groups": bool(separates),
        f"mean_r_within_{group_names[0]}": round(pairwise_mean(corr, first, first), 4),
        f"mean_r_within_{group_names[1]}": round(pairwise_mean(corr, second, second), 4),
        "mean_r_between_groups": round(pairwise_mean(corr, first, second), 4),
        "size_factor_range": f"{size_factors.min():.3f}-{size_factors.max():.3f}",
    }


c1_with = per_comparison[C1]
c2_panel = per_comparison[C2]
sensitivity = pd.DataFrame([
    design_metrics("C1 with A3 (primary)", c1_members, c1_with["log"], c1_with["corr"],
                   c1_with["scores"], c1_with["variance"], c1_with["size_factors"]),
    design_metrics("C1 without A3 (sensitivity)", c1_without, c1_log_wo, c1_corr_wo,
                   c1_scores_wo, c1_var_wo, c1_sf_wo),
    design_metrics("C2 (reference, unaffected by A3)", c2_panel["members"], c2_panel["log"],
                   c2_panel["corr"], c2_panel["scores"], c2_panel["variance"],
                   c2_panel["size_factors"]),
])
# Size factors are only defined up to a common scale, so both sets are
# rescaled to a geometric mean of 1 over the five shared samples before
# comparing; otherwise a uniform rescale would look like a real change.
def rescale_to_unit_geomean(size_factors):
    return size_factors / np.exp(np.log(size_factors).mean())


sf_shift = (rescale_to_unit_geomean(c1_sf_wo)
            / rescale_to_unit_geomean(c1_with["size_factors"][c1_without]) - 1) * 100
print(sensitivity.T.to_string(header=False))
print("Relative size-factor change for remaining C1 samples when A3 is excluded (%):")
print(sf_shift.round(2).to_string())
save_table(sensitivity, "gene_level_A3_sensitivity.csv")

# --- Figure 12: C1 PCA with vs. without A3 ----------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
draw_pca(axes[0], c1_with["scores"], c1_with["variance"], "a  C1 with A3 (n = 6, primary)")
draw_pca(axes[1], c1_scores_wo, c1_var_wo, "b  C1 without A3 (n = 5, sensitivity)")
fig.suptitle("A3 sensitivity: C1 gene-level PCA (each panel normalised and filtered "
             "on its own samples)", x=0.01, ha="left", fontsize=12, fontweight="semibold")
fig.tight_layout()
save_fig(fig, "12_gene_pca_without_A3.png")

# --- Transcript-level vs gene-level (identical Phase 3 method) ---------------
transcript_corr, n_tx_kept = log2_cpm_corr(transcript_counts)
gene_cpm_corr, n_gene_kept = log2_cpm_corr(gene_counts)
peers = [s for s in EXPECTED_SAMPLES if s != OUTLIER_SAMPLE]
level_rows = []
for level, corr, n_kept in [("transcript", transcript_corr, n_tx_kept),
                            ("gene", gene_cpm_corr, n_gene_kept)]:
    r_all = mean_offdiagonal_r(corr)
    within = pd.concat([mean_offdiagonal_r(corr, m) for m in comparison_samples.values()])
    level_rows.append({
        "level": level, "features_after_filter": n_kept,
        "median_r_all_non_A3": round(r_all[peers].median(), 4),
        "A3_r_all": round(r_all[OUTLIER_SAMPLE], 4),
        "median_r_within_comparison_non_A3": round(within[peers].median(), 4),
        "A3_r_within_comparison": round(within[OUTLIER_SAMPLE], 4),
        "min_r_within_comparison_non_A3": round(within[peers].min(), 4),
    })
level_comparison = pd.DataFrame(level_rows)
print("\nTranscript- vs gene-level correlation (log2 CPM, CPM >= 1 in >= 3):")
print(level_comparison.to_string(index=False))
save_table(level_comparison, "transcript_vs_gene_qc_comparison.csv")

# --- Sample-level summary ----------------------------------------------------
within_comparison_r = pd.concat(
    [mean_offdiagonal_r(p["corr"]) for p in per_comparison.values()])[EXPECTED_SAMPLES]
within_sf = pd.concat([p["size_factors"] for p in per_comparison.values()])[EXPECTED_SAMPLES]
pc_scores = pd.concat([p["scores"][["PC1", "PC2"]] for p in per_comparison.values()]).loc[EXPECTED_SAMPLES]
sample_summary = pd.DataFrame({
    "sample_id": EXPECTED_SAMPLES,
    "comparison": sample_comparison.values,
    "experimental_group": sample_groups.values,
    "metadata_qc_flag": metadata["qc_flag"].values,
    "transcript_level_total_counts": library_totals.values,
    "gene_level_total_counts": gene_library_size.values,
    "pct_counts_in_gene_matrix": (gene_library_size / library_totals * 100).round(2).values,
    "pct_counts_TCONS": tcons_share.round(2).values,
    "detected_genes": detected_genes.values,
    f"detected_genes_at_{OUTLIER_SAMPLE}_depth": detected_at_outlier_depth.values,
    "pct_counts_top10_genes": (top10_gene_share * 100).round(2).values,
    "size_factor_global": global_size_factors.round(4).values,
    "size_factor_within_comparison": within_sf.round(4).values,
    "mean_r_within_comparison": within_comparison_r.round(4).values,
    "PC1_within_comparison": pc_scores["PC1"].round(3).values,
    "PC2_within_comparison": pc_scores["PC2"].round(3).values,
})
sample_summary["robust_z_log_gene_library"] = robust_z(np.log10(gene_library_size)).round(2).values
sample_summary["robust_z_detected_at_equal_depth"] = robust_z(detected_at_outlier_depth).round(2).values
sample_summary["robust_z_mean_r_within_comparison"] = robust_z(within_comparison_r).round(2).values
print("\n" + sample_summary.to_string(index=False))
save_table(sample_summary, "gene_level_sample_summary.csv")

# ==================================================================
# Step 8: Final validation and report
# ==================================================================
section("Step 8: Validation")
record_check("Sample summary covers all 12 samples in order",
             sample_summary["sample_id"].tolist() == EXPECTED_SAMPLES)
record_check(f"{OUTLIER_SAMPLE} retained in primary gene-level matrix",
             OUTLIER_SAMPLE in gene_counts.columns)
record_check("Correlation matrix is 12 x 12 with unit diagonal",
             gene_corr.shape == (12, 12) and np.allclose(np.diag(gene_corr), 1.0))
expected_figures = ["07_gene_library_sizes.png", "08_detected_genes.png",
                    "09_gene_expression_distribution.png", "10_gene_correlation_heatmap.png",
                    "11_gene_level_pca.png", "12_gene_pca_without_A3.png"]
missing_figures = [f for f in expected_figures if not (FIG_DIR / f).exists()]
record_check("All Phase 5 figures written", not missing_figures, missing_figures or "6 figures")
counts_sha256_after = sha256_of(COUNTS_FILE)
record_check("Original count file unchanged (SHA-256 before = after)",
             counts_sha256_after == counts_sha256_before, counts_sha256_after)

validation = pd.DataFrame(validation_checks)
save_table(validation, "phase5_validation_report.csv")
n_failed = int((validation["result"] == "FAIL").sum())
print(f"\n{len(validation) - n_failed}/{len(validation)} validation checks passed")

# --- Markdown report (all numbers computed above) ----------------------------
a3_row = sample_summary.set_index("sample_id").loc[OUTLIER_SAMPLE]
non_outlier = sample_summary[sample_summary["sample_id"] != OUTLIER_SAMPLE]
tx_level, gene_level = level_comparison.iloc[0], level_comparison.iloc[1]
sens = sensitivity.set_index("analysis")
c1_with_row = sens.loc["C1 with A3 (primary)"]
c1_wo_row = sens.loc["C1 without A3 (sensitivity)"]
c2_row = sens.loc["C2 (reference, unaffected by A3)"]


def yes_no(flag):
    return "yes" if flag else "no"


report = f"""# Phase 5: Transcript-to-gene aggregation and gene-level QC

*Generated by `data/raw/notebooks/05_transcript_to_gene_aggregation.py` on
{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}.
All numbers below are computed by the script.*

**Validation: {len(validation) - n_failed}/{len(validation)} checks passed**
(`results/tables/phase5_validation_report.csv`).

## Annotation

- **{ANNOTATION_LABEL}**, `{GTF_FILE.name}` (SHA-256 `{gtf_sha256[:16]}…`).
- GEO states the assembly is GRCh38/hg38 (HISAT2 2.0.4 alignment, StringTie
  1.3.3b quantification). No annotation release is given on GEO.
- The {n_enst:,} ENST IDs in the count file are **exactly** the Ensembl 97
  transcript set: {len(only_in_file)} IDs are only in the file and
  {len(only_in_annotation)} only in the annotation. Releases 96, 98, 99 and
  101 do not match (see `data/reference/ANNOTATION_SOURCE.md`).

## Mapping

| Metric | Value |
|---|---|
| Total transcripts | {n_total:,} |
| ENST transcripts | {n_enst:,} |
| TCONS transcripts | {n_tcons:,} |
| ENST mapped to exactly one gene | {n_mapped:,} |
| ENST unmapped | {n_enst_unmapped:,} |
| ENST ambiguous (more than one gene) | {n_ambiguous:,} |
| Unique genes in matrix | {gene_counts.shape[0]:,} |
| Genes with any counts | {int((gene_counts.sum(axis=1) > 0).sum()):,} |
| TCONS annotated | 0 |
| Original counts in gene-level matrix | {pct_counts_in_genes:.2f}% |

IDs had no version suffixes ({n_with_version_suffix} found). Matching
nevertheless strips `.N` suffixes from both sides and checks that this causes
no collisions.

**TCONS.** The {n_tcons:,} novel StringTie transcripts
({int((tcons_counts.sum(axis=1) > 0).sum()):,} detected in at least one
sample) carry {total_tcons_counts / total_original_counts:.2%} of all counts.
That is {tcons_share.drop(OUTLIER_SAMPLE).min():.1f}–{tcons_share.drop(OUTLIER_SAMPLE).max():.1f}%
per sample, and {tcons_share[OUTLIER_SAMPLE]:.1f}% in {OUTLIER_SAMPLE}. GEO
publishes no assembly GTF or transcript sequences, so they cannot be placed
against Ensembl genes. They are excluded from the gene matrix and listed in
`unmapped_transcripts.csv`.

## Gene-level QC

Normalisation (QC only): median-of-ratios size factors, log2(normalised + 1),
genes kept if normalised count ≥ {MIN_NORM_COUNT} in ≥ {MIN_SAMPLES} samples,
and PCA on the top {N_TOP_VARIABLE} variable genes. PCA and normalisation are
done **per comparison**. The 12-sample heatmap and distribution plots are
descriptive only.

**Transcript-level vs. gene-level**, same method as Phase 3
(log2 CPM, CPM ≥ 1 in ≥ 3 samples):

| Level | Features | Median r, non-A3 (all samples) | Median r, non-A3 (within comparison) | A3 r (within C1) |
|---|---|---|---|---|
| Transcript | {int(tx_level['features_after_filter']):,} | {tx_level['median_r_all_non_A3']:.3f} | {tx_level['median_r_within_comparison_non_A3']:.3f} | {tx_level['A3_r_within_comparison']:.3f} |
| Gene | {int(gene_level['features_after_filter']):,} | {gene_level['median_r_all_non_A3']:.3f} | {gene_level['median_r_within_comparison_non_A3']:.3f} | {gene_level['A3_r_within_comparison']:.3f} |

**{OUTLIER_SAMPLE}**, gene level:
- {a3_row['gene_level_total_counts'] / 1e6:.1f} M counts; the non-A3 median is
  {non_outlier['gene_level_total_counts'].median() / 1e6:.1f} M.
- {int(a3_row[f'detected_genes_at_{OUTLIER_SAMPLE}_depth']):,} genes detected at its own depth.
  Other samples downsampled to the same depth detect
  {int(non_outlier[f'detected_genes_at_{OUTLIER_SAMPLE}_depth'].min()):,}–{int(non_outlier[f'detected_genes_at_{OUTLIER_SAMPLE}_depth'].max()):,}.
- {a3_row['pct_counts_top10_genes']:.1f}% of counts fall in its top 10 genes
  (others: {non_outlier['pct_counts_top10_genes'].min():.1f}–{non_outlier['pct_counts_top10_genes'].max():.1f}%).
- Mean r within C1 is {a3_row['mean_r_within_comparison']:.3f}; the other samples
  range {non_outlier['mean_r_within_comparison'].min():.3f}–{non_outlier['mean_r_within_comparison'].max():.3f}.

## A3 sensitivity (C1)

| | C1 with A3 | C1 without A3 | C2 (reference) |
|---|---|---|---|
| Genes after filter | {c1_with_row['genes_after_filter']:,} | {c1_wo_row['genes_after_filter']:,} | {c2_row['genes_after_filter']:,} |
| PC1 % variance | {c1_with_row['PC1_pct_variance']} | {c1_wo_row['PC1_pct_variance']} | {c2_row['PC1_pct_variance']} |
| PC1 separates the two groups | {yes_no(c1_with_row['PC1_separates_groups'])} | {yes_no(c1_wo_row['PC1_separates_groups'])} | {yes_no(c2_row['PC1_separates_groups'])} |
| Mean r between groups | {c1_with_row['mean_r_between_groups']} | {c1_wo_row['mean_r_between_groups']} | {c2_row['mean_r_between_groups']} |

Within-group mean r for C1 control (Ethanol): {c1_with_row['mean_r_within_hESC_Ethanol']} with A3,
{c1_wo_row['mean_r_within_hESC_Ethanol']} without. Relative size factors of the
remaining C1 samples change by at most {sf_shift.abs().max():.2f}% when A3 is
excluded (compared after rescaling both sets to a geometric mean of 1).
Full table: `gene_level_A3_sensitivity.csv`.
These are descriptive QC metrics, not tests of treatment effects.

## Limitations

- TCONS transcripts ({total_tcons_counts / total_original_counts:.1%} of counts)
  cannot be annotated without the submitters' StringTie merged GTF. Reads from
  unannotated loci, or from novel isoforms of known genes, are missing from the
  gene matrix.
- Summing StringTie transcript counts to genes is a plain sum. StringTie
  transcript counts are estimates for reads that fit several isoforms, so
  gene totals inherit that estimation, not unique read counts.
- QC normalisation here is exploratory and is not the model to be used for
  differential expression.
- No batch, RNA-integrity or library-preparation metadata are available.

## Status

{"All validation checks passed. The gene-level matrix (`results/tables/gene_level_counts.csv`) is ready for Phase 6, to be analysed **per comparison** with A3 retained in the primary C1 analysis and a without-A3 sensitivity analysis." if n_failed == 0 else f"**{n_failed} validation check(s) FAILED - do not proceed to Phase 6 until resolved.**"}

Software: Python {platform.python_version()}, pandas {pd.__version__}, numpy {np.__version__},
scikit-learn {sklearn.__version__}, matplotlib {matplotlib.__version__}. Random seed {RANDOM_SEED}.
"""
REPORT_FILE.write_text(report, encoding="utf-8")
print(f"  saved report: {REPORT_FILE.relative_to(REPO_ROOT)}")

if n_failed:
    print(validation[validation["result"] == "FAIL"].to_string(index=False))
    sys.exit(f"ERROR: {n_failed} validation check(s) failed.")
print("\nPhase 5 completed successfully. No differential expression performed.")
