"""
Phase 4: Sample metadata and experimental design
Source: NCBI GEO accession GSE213410
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213410
        (sample records GSM6585843-GSM6585854)

Design
  Comparison 1 (C1) - H9 hESCs, normal karyotype, 24 h
      A1-A3  Ethanol (vehicle control)
      A4-A6  5-OXO-ETE
  Comparison 2 (C2) - P hiPSCs, trisomy 21, 24 h
      A7-A9   DMSO (vehicle control)
      A10-A12 Calcitriol

What this script does
  1. Loads the curated metadata (data/metadata/sample_metadata.csv)
  2. Validates it: required columns, no gaps, 12 unique samples, A1..A12 order,
     exact match to the count-matrix columns, 3 replicates per group,
     one control + one treatment group per comparison, A3 flagged
  3. Joins Phase 3 QC metrics to the metadata
  4. Checks whether the Phase 3 sample structure agrees with the design
  5. Writes summary tables to results/tables/

No differential expression is performed. The count dataset is only read
(header row only), never modified.

Run from anywhere:
    python data/raw/notebooks/04_sample_metadata.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
METADATA_FILE = REPO_ROOT / "data/metadata/sample_metadata.csv"
COUNTS_FILE = REPO_ROOT / "data/raw/GSE213410_transcripts.readcount.txt.gz"
TABLE_DIR = REPO_ROOT / "results/tables"
QC_METRICS_FILE = TABLE_DIR / "sample_qc_metrics.csv"
CORR_FILE = TABLE_DIR / "correlation_matrix_pearson_log2cpm.csv"
PCA_FILE = TABLE_DIR / "pca_sample_scores.csv"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_SAMPLES = [f"A{i}" for i in range(1, 13)]
REQUIRED_COLUMNS = [
    "sample_id", "geo_accession", "cell_type", "genotype", "treatment",
    "experimental_group", "biological_replicate", "comparison",
]
EXPECTED_DESIGN = {   # verified assignments (NCBI GEO GSE213410)
    "A1": "Ethanol", "A2": "Ethanol", "A3": "Ethanol",
    "A4": "5-OXO-ETE", "A5": "5-OXO-ETE", "A6": "5-OXO-ETE",
    "A7": "DMSO", "A8": "DMSO", "A9": "DMSO",
    "A10": "Calcitriol", "A11": "Calcitriol", "A12": "Calcitriol",
}
OUTLIER_SAMPLE = "A3"

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 40)


def section(title):
    print("\n" + "=" * 64)
    print(title)
    print("=" * 64)


def save_table(df, name, index=False):
    path = TABLE_DIR / name
    df.to_csv(path, index=index)
    print(f"  saved table: {path.relative_to(REPO_ROOT)}")


# ------------------------------------------------------------------
# Step 1: Load metadata
# qc_note is free text and empty for most samples, so keep empty strings
# instead of NaN to avoid tripping the missing-value check.
# ------------------------------------------------------------------
section("Step 1: Load metadata")
meta = pd.read_csv(METADATA_FILE, dtype={"biological_replicate": int},
                   keep_default_na=False, na_values=[""])
meta["qc_note"] = meta["qc_note"].fillna("")
print(f"Loaded {len(meta)} rows x {meta.shape[1]} columns from "
      f"{METADATA_FILE.relative_to(REPO_ROOT)}")
print(meta[["sample_id", "geo_accession", "cell_type", "genotype", "treatment",
            "experimental_group", "biological_replicate", "comparison",
            "qc_flag"]].to_string(index=False))

# ------------------------------------------------------------------
# Step 2: Validation. Every check is recorded; the script stops at the
# end if any check fails, so problems can't slip through silently.
# ------------------------------------------------------------------
section("Step 2: Validation")
checks = []


def check(name, passed, detail=""):
    checks.append({"check": name, "result": "PASS" if passed else "FAIL",
                   "detail": detail})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


missing_cols = [c for c in REQUIRED_COLUMNS if c not in meta.columns]
check("All required columns present", not missing_cols,
      f"missing: {missing_cols}" if missing_cols else f"{len(REQUIRED_COLUMNS)} required")

required_na = int(meta[REQUIRED_COLUMNS].isna().sum().sum())
check("No missing values in required columns", required_na == 0,
      f"{required_na} missing")

ids = meta["sample_id"].tolist()
dupes = sorted({s for s in ids if ids.count(s) > 1})
check("Each sample appears exactly once", not dupes and len(ids) == 12,
      f"{len(ids)} rows, duplicates: {dupes or 'none'}")

absent = [s for s in EXPECTED_SAMPLES if s not in ids]
extra = [s for s in ids if s not in EXPECTED_SAMPLES]
check("All 12 expected samples (A1-A12) present, no extras", not absent and not extra,
      f"absent: {absent or 'none'}, unexpected: {extra or 'none'}")

check("Samples ordered A1 through A12", ids == EXPECTED_SAMPLES,
      "order: " + ", ".join(ids))

gsm_ok = meta["geo_accession"].str.fullmatch(r"GSM\d+").all() and meta["geo_accession"].is_unique
check("GEO accessions well-formed and unique", bool(gsm_ok),
      f"{meta['geo_accession'].iloc[0]} .. {meta['geo_accession'].iloc[-1]}")

# Count matrix: read only the header row to get the sample columns.
count_samples = pd.read_csv(COUNTS_FILE, sep="\t", compression="gzip", nrows=0).columns[1:]
only_counts = sorted(set(count_samples) - set(ids))
only_meta = sorted(set(ids) - set(count_samples))
check("Metadata samples match count-matrix columns", not only_counts and not only_meta,
      f"count matrix has {len(count_samples)} samples; "
      f"only in counts: {only_counts or 'none'}, only in metadata: {only_meta or 'none'}")

wrong_treat = [s for s, t in EXPECTED_DESIGN.items()
               if s in ids and meta.set_index("sample_id").at[s, "treatment"] != t]
check("Treatments match the verified GEO assignments", not wrong_treat,
      f"mismatches: {wrong_treat or 'none'}")

reps = meta.groupby("experimental_group")["biological_replicate"].apply(sorted)
bad_reps = [g for g, r in reps.items() if r != [1, 2, 3]]
check("Each group has replicates 1, 2, 3", not bad_reps,
      f"{len(reps)} groups; problems: {bad_reps or 'none'}")

roles = meta.groupby("comparison")["treatment_role"].apply(lambda r: sorted(set(r)))
bad_roles = [c for c, r in roles.items() if r != ["control", "treatment"]]
check("Each comparison has one control and one treatment group", not bad_roles,
      f"{len(roles)} comparisons; problems: {bad_roles or 'none'}")

# Within a comparison, cell type and genotype must be constant; otherwise
# the treatment effect would be confounded with cell background.
uniform = meta.groupby("comparison")[["cell_type", "genotype", "cell_line"]].nunique()
check("Cell type, line and genotype constant within each comparison",
      bool((uniform == 1).all().all()))

flagged = meta.loc[meta["qc_flag"] == "outlier", "sample_id"].tolist()
check(f"{OUTLIER_SAMPLE} preserved and flagged as QC outlier",
      OUTLIER_SAMPLE in ids and flagged == [OUTLIER_SAMPLE],
      f"flagged samples: {flagged or 'none'}")

validation = pd.DataFrame(checks)
save_table(validation, "metadata_validation_report.csv")
n_fail = int((validation["result"] == "FAIL").sum())
print(f"\n{len(validation) - n_fail}/{len(validation)} checks passed")
if n_fail:
    raise SystemExit("Metadata validation failed - fix the issues above before continuing.")

# ------------------------------------------------------------------
# Step 3: Join Phase 3 QC metrics so each sample carries its QC profile
# ------------------------------------------------------------------
section("Step 3: Metadata with QC metrics")
annotated = meta.copy()
if QC_METRICS_FILE.exists():
    qc = pd.read_csv(QC_METRICS_FILE, index_col=0)
    qc_cols = ["million_reads", "detected_transcripts",
               "pct_reads_top10_transcripts", "mean_pearson_r_to_others"]
    annotated = annotated.join(qc[qc_cols], on="sample_id")
    print(annotated[["sample_id", "experimental_group", "biological_replicate",
                     "qc_flag"] + qc_cols].to_string(index=False))
else:
    print("  Phase 3 QC metrics not found - run 03_rnaseq_quality_control.py first.")
save_table(annotated, "sample_metadata_annotated.csv")

# ------------------------------------------------------------------
# Step 4: Experimental design summary (one row per experimental group)
# "usable" = replicates that passed QC; shows the effect of the A3 flag.
# ------------------------------------------------------------------
section("Step 4: Experimental design summary")
agg = {
    "cell_type": ("cell_type", "first"),
    "cell_line": ("cell_line", "first"),
    "genotype": ("genotype", "first"),
    "treatment": ("treatment", "first"),
    "role": ("treatment_role", "first"),
    "n_replicates": ("sample_id", "size"),
    "n_passing_qc": ("qc_flag", lambda f: int((f == "pass").sum())),
    "samples": ("sample_id", ", ".join),
    "geo_accessions": ("geo_accession", ", ".join),
}
if "million_reads" in annotated:
    agg["mean_million_reads"] = ("million_reads", lambda x: round(x.mean(), 2))
design = (annotated.groupby(["comparison", "experimental_group"], sort=False)
          .agg(**agg).reset_index())
design["qc_flagged"] = [
    ", ".join(annotated.loc[(annotated["experimental_group"] == g)
                            & (annotated["qc_flag"] != "pass"), "sample_id"]) or "-"
    for g in design["experimental_group"]
]
print(design.drop(columns=["geo_accessions"]).to_string(index=False))
save_table(design, "experimental_design_summary.csv")

overview = pd.DataFrame({
    "item": ["GEO series", "Organism", "Samples", "Comparisons", "Groups",
             "Replicates per group", "Treatment duration", "QC-flagged samples"],
    "value": ["GSE213410", "Homo sapiens", str(len(meta)),
              str(meta["comparison"].nunique()),
              str(meta["experimental_group"].nunique()),
              "3", "24 h", ", ".join(flagged)],
})
print("\n" + overview.to_string(index=False))
save_table(overview, "metadata_overview.csv")

# ------------------------------------------------------------------
# Step 5: Does the Phase 3 sample structure agree with the design?
# Compare average correlation within vs. between comparisons and groups,
# and PC1 position per comparison. Exploratory only - not a DE test.
# ------------------------------------------------------------------
section("Step 5: Consistency of Phase 3 QC with the design")
if CORR_FILE.exists() and PCA_FILE.exists():
    corr = pd.read_csv(CORR_FILE, index_col=0)
    pca = pd.read_csv(PCA_FILE, index_col=0)
    m = meta.set_index("sample_id")

    rows = []
    for i, a in enumerate(EXPECTED_SAMPLES):
        for b in EXPECTED_SAMPLES[i + 1:]:
            same_comp = m.at[a, "comparison"] == m.at[b, "comparison"]
            same_group = m.at[a, "experimental_group"] == m.at[b, "experimental_group"]
            relation = ("same group" if same_group
                        else "same comparison, different treatment" if same_comp
                        else "different comparison (cell type)")
            rows.append({"relation": relation, "r": corr.at[a, b],
                         "involves_A3": OUTLIER_SAMPLE in (a, b)})
    pairs = pd.DataFrame(rows)
    consistency = (pairs.groupby(["involves_A3", "relation"])["r"]
                   .agg(n_pairs="size", mean_r="mean", min_r="min", max_r="max")
                   .round(4).reset_index())
    print("Pairwise Pearson r (log2 CPM) by sample relationship:")
    print(consistency.to_string(index=False))
    save_table(consistency, "design_correlation_consistency.csv")

    pc_by_comp = (pca[["PC1", "PC2"]].join(m["comparison"])
                  .groupby("comparison").agg(["mean", "min", "max"]).round(2))
    print("\nPC1/PC2 scores by comparison:")
    print(pc_by_comp.to_string())
else:
    print("  Phase 3 correlation/PCA tables not found - skipped.")

print("\nMetadata phase completed. No differential expression performed.")
