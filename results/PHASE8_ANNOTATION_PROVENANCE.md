# Phase 8, Step 2: annotation source and provenance

**Purpose:** official Ensembl release 97 annotation for the eight target genes (nine Ensembl gene IDs, because MATR3 has two).
**Script:** `data/raw/notebooks/08_gene_annotation_retrieval.py`
**Retrieved:** 2026-11-26 (UTC). Exact per-file times are in `results/tables/phase8_raw/ensembl97_source_files.csv`.

## Source

| Item | Value |
|---|---|
| Release | Ensembl 97 (July 2019), GRCh38.p12, the same release as the local GTF |
| Database | `homo_sapiens_core_97_38`, the release-97 human core database |
| Location | https://ftp.ensembl.org/pub/release-97/mysql/homo_sapiens_core_97_38/ |
| Files used | schema `homo_sapiens_core_97_38.sql.gz` and the tables `gene`, `transcript`, `translation`, `seq_region`, `attrib_type`, `gene_attrib`, `transcript_attrib`, `external_db`, `xref`, `object_xref`, `external_synonym` |
| File dates on server | 3 Sep 2019 (all files) |
| Integrity | Every file matches the BSD `sum` value in Ensembl's published `CHECKSUMS` file for that directory. SHA-256 values are also recorded. |
| Current Ensembl records used | None |

### Why the database files and not the archive website

The requested source was the July 2019 archive website, `jul2019.archive.ensembl.org`. That site has been retired. On 2026-11-26, every URL on it was redirected to https://www.ensembl.org/help/articles/archives. That included the home page, gene pages and the BioMart service at `/biomart/martservice`.

The release-97 core database is the data the archive website displayed. Ensembl still publishes it on its FTP site, with the original 2019 file dates and checksums. It is therefore the same release from the same provider.

## Reproducing

```
python data/raw/notebooks/08_gene_annotation_retrieval.py [cache_dir]
```

- **Cache:** the full tables (117.0 MB, 116,981,509 bytes) are downloaded to a cache outside the repository. By default this is `%LOCALAPPDATA%\bioinformatics-portfolio-cache\ensembl97_core_db`. They are not committed. The script re-downloads any missing file and stops if any checksum fails.
- **GTF cross-check:** the script checks each gene's ID, version, biotype, chromosome, start, end and strand against `data/reference/Homo_sapiens.GRCh38.97.gtf.gz`. It stops on any difference. All 9 genes matched.
- **No overwriting:** the script refuses to run if any of its output files already exist.
- **Validation:** `python data/raw/notebooks/08_validate_gene_annotation.py [cache_dir]` repeats the checks and writes `results/tables/phase8_annotation_validation.csv`. On 2026-11-26, 35/35 checks passed.

## Outputs

| File | Content |
|---|---|
| `results/tables/phase8_raw/ensembl97_gene_records.csv` | The 9 gene rows exactly as stored in the database, plus chromosome name and display symbol |
| `results/tables/phase8_raw/ensembl97_transcripts.csv` | All 71 transcripts of those genes, with biotype, coordinates, canonical flag and protein ID |
| `results/tables/phase8_raw/ensembl97_attributes.csv` | All 660 gene and transcript attributes, such as readthrough flag, Havana curation terms, TSL and MANE/Ensembl Select |
| `results/tables/phase8_raw/ensembl97_xrefs.csv` | All 1,183 external references on the genes, transcripts and proteins, with synonyms |
| `results/tables/phase8_raw/ensembl97_source_files.csv` | Source URL, server date, retrieval time, size, published and observed checksum, and SHA-256 for each file |
| `results/tables/phase8_target_gene_annotation.csv` | Curated table: one row per gene ID with the main official fields. It contains annotation facts only and no biological interpretation. |
| `results/tables/phase8_annotation_validation.csv` | 35 validation checks, each with expected value, actual value and result, written by `data/raw/notebooks/08_validate_gene_annotation.py`. Covers checksums of the 13 downloaded database files, comparison with Ensembl's current published `CHECKSUMS`, the retired archive website, the GTF checksum, the 9 gene records and 71 transcripts against the GTF, row counts, the curated table and Phase 5 consistency, and a check for local paths. |

The raw files are unedited extracts. The curated table is derived from them by the same script. Biological interpretation will be kept in a separate file.
