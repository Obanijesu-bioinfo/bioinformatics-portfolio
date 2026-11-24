"""
Phase 2: RNA-seq dataset exploration
Dataset: GSE213410_transcripts.readcount.txt.gz

Loads the transcript-level read count matrix and reports its structure,
sample columns, data quality checks, and summary statistics.
The original compressed file is only read, never modified.

Run from the repository root:
    python data/raw/notebooks/02_dataset_exploration.py
"""

from pathlib import Path

import pandas as pd

DATA_FILE = Path("data/raw/GSE213410_transcripts.readcount.txt.gz")
ID_COLUMN = "transcriptID"

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)


def section(title):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)


# ------------------------------------------------------------------
# Step 1: Load the complete dataset
# pandas decompresses .gz files on the fly, so the original file is
# read-only here. The file is tab-separated. The ID column is kept as a
# regular column (not the index) so duplicates can be checked later.
# ------------------------------------------------------------------
section("Step 1: Load dataset")
counts = pd.read_csv(DATA_FILE, sep="\t", compression="gzip")
print(f"Loaded: {DATA_FILE}")

# ------------------------------------------------------------------
# Step 2: Dimensions
# ------------------------------------------------------------------
section("Step 2: Dataset dimensions")
n_rows, n_cols = counts.shape
print(f"Rows (transcripts): {n_rows:,}")
print(f"Columns (ID + samples): {n_cols}")

# ------------------------------------------------------------------
# Step 3: First 10 rows
# ------------------------------------------------------------------
section("Step 3: First 10 rows")
print(counts.head(10).to_markdown(index=False))

# ------------------------------------------------------------------
# Step 4: Sample column names
# Every column except the identifier column is a sample.
# ------------------------------------------------------------------
section("Step 4: Sample columns")
sample_cols = [c for c in counts.columns if c != ID_COLUMN]
print(f"Identifier column: {ID_COLUMN}")
print(f"Number of samples: {len(sample_cols)}")
print("Samples:", ", ".join(sample_cols))

# ------------------------------------------------------------------
# Step 5: Missing values and duplicate identifiers
# ------------------------------------------------------------------
section("Step 5: Data quality checks")
missing = counts.isna().sum()
print(f"Total missing values: {int(missing.sum())}")
if missing.sum() > 0:
    print(missing[missing > 0].to_string())

dup_mask = counts[ID_COLUMN].duplicated(keep=False)
n_dup_ids = counts.loc[dup_mask, ID_COLUMN].nunique()
print(f"Duplicate identifiers: {n_dup_ids} IDs across {int(dup_mask.sum())} rows")
if n_dup_ids > 0:
    print(counts.loc[dup_mask].sort_values(ID_COLUMN).head(10).to_markdown(index=False))

non_numeric = [c for c in sample_cols if not pd.api.types.is_numeric_dtype(counts[c])]
print(f"Non-numeric sample columns: {non_numeric if non_numeric else 'none'}")

# Identifier source: Ensembl (ENST) vs. novel assembled transcripts (TCONS)
id_prefix = counts[ID_COLUMN].str.extract(r"^([A-Za-z]+)", expand=False)
print("\nIdentifier prefixes:")
print(id_prefix.value_counts().to_string())

# ------------------------------------------------------------------
# Step 6: Basic dataset summary
# ------------------------------------------------------------------
section("Step 6: Dataset summary")
sample_data = counts[sample_cols]

per_sample = pd.DataFrame({
    "total_reads": sample_data.sum(),
    "detected_transcripts(>0)": (sample_data > 0).sum(),
    "mean_count": sample_data.mean().round(2),
    "median_count": sample_data.median(),
    "max_count": sample_data.max(),
})
print("Per-sample statistics:")
print(per_sample.to_markdown(floatfmt=("", ",.0f", ",.0f", ",.2f", ",.0f", ",.0f")))

row_totals = sample_data.sum(axis=1)
total_cells = sample_data.size
zero_cells = int((sample_data == 0).sum().sum())
print("\nOverall statistics:")
print(f"Total reads across all samples: {int(row_totals.sum()):,}")
print(f"Zero entries: {zero_cells:,} of {total_cells:,} ({zero_cells / total_cells:.1%})")
print(f"Transcripts with zero counts in all samples: {int((row_totals == 0).sum()):,}")
print(f"Transcripts detected in all samples: {int((sample_data > 0).all(axis=1).sum()):,}")
print(f"Memory usage: {counts.memory_usage(deep=True).sum() / 1e6:.1f} MB")

print("\nDataset exploration completed.")
