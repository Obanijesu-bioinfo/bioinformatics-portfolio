# Reference annotation: provenance

| Field | Value |
|---|---|
| File | `Homo_sapiens.GRCh38.97.gtf.gz` |
| Source | Ensembl, release 97 (July 2019) |
| Genome assembly | GRCh38.p12 (`#!genome-build GRCh38.p12`, `NCBI:GCA_000001405.27`) |
| Gene build last updated | 2019-03 (GTF header) |
| Download URL | https://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz |
| Downloaded | 2026-11-24 (UTC) |
| Size | 46,425,997 bytes |
| SHA-256 | `8d2b55f734f4e351b37d24eddba7d955034c98390564260318341ee464590a19` |
| Ensembl CHECKSUMS (BSD `sum`) | `46698 45338` (matches the value published in the release-97 `CHECKSUMS` file) |
| Equivalent GENCODE release | GENCODE 31 (same gene build; GENCODE uses versioned IDs and duplicates chrY PAR features) |

## Why this annotation

1. **Assembly.** The GEO sample records for GSE213410 (GSM6585843 to GSM6585854)
   state `Assembly: Homo Sapiens(GRCh38/hg38)`. Reads were aligned with
   HISAT2 2.0.4 and quantified with StringTie 1.3.3b. No annotation release is
   stated on GEO.
2. **Release.** The count file's 226,788 Ensembl transcript IDs carry no version
   suffix. StringTie reports every reference transcript (74,288 rows are zero in
   all samples), so the ID set should equal the annotation's transcript set. I
   compared the IDs against Ensembl GRCh38 GTFs:

   | Ensembl release | Transcripts in GTF | IDs only in count file | IDs only in GTF |
   |---|---|---|---|
   | 96 | 208,527 | 18,684 | 423 |
   | **97** | **226,788** | **0** | **0** |
   | 98 | 227,368 | 58 | 638 |
   | 99 | 227,818 | 138 | 1,168 |
   | 101 | 229,487 | 705 | 3,404 |

   Only release 97 matches exactly, so it is the annotation the submitters used.
   This comparison was run once during setup. `05_transcript_to_gene_aggregation.py`
   checks the exact release-97 match again on every run.

## TCONS transcripts

The 27,832 `TCONS_*` IDs are novel transcripts from StringTie's merge step. GEO
provides only the count table (`GSE213410_transcripts.readcount.txt.gz`), not
the merged assembly GTF or transcript sequences. Without coordinates or
sequences, TCONS IDs cannot be placed against Ensembl genes, so they are left
unannotated. Raw reads are linked through BioProject PRJNA880672, but
re-assembling them would not reproduce the same TCONS IDs.
