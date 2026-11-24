
import gzip
from pathlib import Path

# Locate the original RNA-seq dataset
data_file = Path("data/raw/GSE213410_transcripts.readcount.txt.gz")

print("RNA-seq Dataset Inspection")
print("=" * 40)

# Read the compressed file
with gzip.open(data_file, "rt", encoding="utf-8") as file:
    for line_number, line in enumerate(file):
        if line_number < 10:
            print(line.rstrip())

print("\nDataset inspection completed.")

