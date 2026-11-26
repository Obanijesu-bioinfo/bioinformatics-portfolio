# Phase 8, Step 3: functional annotation source and provenance

**Purpose:** UniProtKB functional annotation for the nine Ensembl 97 gene IDs. The eight targets include two MATR3 records, which are kept separate.
**Script:** `data/raw/notebooks/08b_uniprot_functional_annotation.py`
**Retrieved:** 2026-11-26 (UTC). Every request's URL, time, HTTP status, release header and SHA-256 is in `results/tables/phase8_raw/uniprot_2026_04/request_log.csv`.
**Scope:** annotation facts only. Biological interpretation is reported separately in Step 4 (`results/PHASE8_BIOLOGICAL_INTERPRETATION.md`).

## Sources

| Source | Release | Access | Check |
|---|---|---|---|
| UniProtKB, human records only (taxon 9606) | **2026_04**, released 4 Nov 2026 | REST API, https://rest.uniprot.org | Every response (39 in total) carried `X-UniProt-Release: 2026_04`. The script stops if any response comes from another release. |
| UniParc | 2026_04 | REST API | Same release header check |
| Ensembl 97 cross-references and transcripts | Release 97 | Step 2 extracts in `results/tables/phase8_raw/` | See `PHASE8_ANNOTATION_PROVENANCE.md` |
| Ensembl 97 protein sequences | Release 97 (server file date 25 May 2019) | `https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz` | Matches the BSD `sum` value in Ensembl's published `CHECKSUMS` file. The FASTA contains all 51 target translations, with the same versions as the Step 2 transcript table. |

**Rerunning after 2026_04 is replaced:** the UniProt API serves only the current release. Once 2026_04 is replaced, a rerun of the script will stop at the release check rather than mix releases. The saved raw responses are the record of this retrieval.

## Mapping rules

No mapping was made by gene symbol.

1. **Candidates come only from explicit cross-references.** A UniProt accession is a candidate for an Ensembl 97 gene only through one of these:
   - (a) an Ensembl 97 cross-reference from one of the gene's proteins (`Uniprot/SWISSPROT` or `Uniprot/SPTREMBL`);
   - (b) the UniProt 2026_04 entry's own Ensembl cross-reference naming that gene ID.
2. **Name-based links are excluded.** Ensembl's gene-level `Uniprot_gn` links are assigned by gene name, so they were not used.
3. **Reviewed records are preferred.**
   - If a gene has a reviewed (Swiss-Prot) candidate, only reviewed records are used.
   - Unreviewed (TrEMBL) records are used only for genes with no reviewed candidate, and are labelled as unreviewed.
4. **Every pair is checked by sequence.** For each gene and accession, the Ensembl 97 protein sequences are compared for exact identity with the UniProt canonical and isoform sequences.
5. **Links to another target's transcripts are not transferred.** If an entry cross-references a transcript that belonged to a different target gene in Ensembl 97, that link is recorded as ambiguous and not used.

## Results by gene

| Gene ID | Record used | Review status | Notes |
|---|---|---|---|
| RPSAP15, ENSG00000237506.3 | none | — | Processed pseudogene. No UniProt cross-reference in Ensembl 97 or UniProt 2026_04, so no protein annotation. |
| TUBBP1, ENSG00000127589.4 | none | — | Transcribed processed pseudogene. No cross-reference in either source, so no protein annotation. |
| INHBE, ENSG00000139269.3 | P58166 | reviewed | Sequence identical to canonical ENSP00000266646.2 |
| AC069288.1, ENSG00000286192.1 | A0A3B3ITW8 | **unreviewed** | See below |
| SFT2D2, ENSG00000213064.10 | O95562 | reviewed | Sequence identical to canonical ENSP00000271375.3 |
| MATR3, ENSG00000280987.4 (primary) | A0A0R4J2E8, plus 7 further unreviewed records | **unreviewed** | See below |
| MATR3, ENSG00000015479.18 | P43243 | reviewed | Identical to ENSP00000378284.3 and ENSP00000482895.1. ENSP00000423587.1 is identical to isoform P43243-2. |
| STEAP2, ENSG00000157214.14 | Q8NFT2 | reviewed | Sequence identical to canonical ENSP00000378119.2 and other translations |
| EOGT, ENSG00000163378.14 | Q5NDL2 | reviewed | Identical to canonical ENSP00000373206.3 and ENSP00000443780.2. Two other translations are identical to isoform Q5NDL2-3. |

### AC069288.1 (ENSG00000286192)

- **Ensembl 97 made no UniProtKB link** for this gene, only a UniParc sequence identifier (UPI0000168663).
- **UniProt 2026_04 links one unreviewed record to the gene ID: A0A3B3ITW8.**
  - Its sequence (155 amino acids) is identical to the gene's only Ensembl 97 protein, ENSP00000498895.1.
  - UniProt gives it the gene name **MRM2** and the *submitted* (not curated) name "Mitochondrial rRNA methyltransferase 2".
- **Annotation on that record:**
  - no function text;
  - no subcellular location;
  - one GO term, *methyltransferase activity* (GO:0008168), assigned automatically from InterPro (`IEA:InterPro`).
- **Other identical records:** UniParc lists two further unreviewed records with the identical sequence, A0ACP2BFA1 and A0ACP2EWM7. Neither is cross-referenced to this gene, so neither was used.

### MATR3: the two records are handled separately

- **ENSG00000015479 (overlapping gene model)** maps to the reviewed record **P43243** through explicit Ensembl 97 cross-references, confirmed by sequence identity.
- **ENSG00000280987 (primary record)** maps only to unreviewed records:
  - **A0A0R4J2E8** is the primary one. Ensembl 97 cross-references it to this gene's canonical protein ENSP00000354346.4, UniProt 2026_04 also cross-references it to this gene, and the sequences are identical.
  - Seven more unreviewed records are each linked by Ensembl 97 to one of this gene's proteins, with identical sequences: A8MXP9, B3KM87, D6R8Z5, D6R9N0, D6RAM9, D6RB45 and D6RCM3. UniProt 2026_04 does not cross-reference any of them to this gene ID.
  - B3KM87 is linked by Ensembl 97 to a protein of each MATR3 gene ID: ENSP00000421218 (ENSG00000280987) and ENSP00000426030 (ENSG00000015479). Its sequence is identical to both of these Ensembl 97 proteins, which share UniParc ID UPI000015F327. It is used only for ENSG00000280987, which has no reviewed record. For ENSG00000015479 it is not used, because the reviewed record P43243 is available.
- **P43243 is not used for ENSG00000280987.** The evidence is recorded as *ambiguous – not transferred* in `phase8_uniprot_mapping.csv`:
  - The canonical protein of ENSG00000280987 (ENSP00000354346.4) has exactly the same sequence as P43243. Both carry UniParc ID UPI0000000DEE, and so does A0A0R4J2E8.
  - UniProt 2026_04 cross-references P43243-1 to two transcripts that belonged to ENSG00000280987 in Ensembl 97, ENST00000361059.7 and ENST00000509990.5. It assigns both transcripts to ENSG00000015479.
  - Ensembl 97 linked these proteins to A0A0R4J2E8, not to P43243.
- **The function text on A0A0R4J2E8 is automatic annotation.** It matches the reviewed MATR3 text, but its evidence is `ECO:0000256` with source ARBA (UniProt's automatic rule-based annotation). Its three GO terms are all electronic: `IEA:UniProtKB-SubCell`, `IEA:UniProtKB-UniRule` and `IEA:InterPro`. It is reported as unreviewed and automatic.

## Evidence codes in the tables

| Code | Meaning |
|---|---|
| ECO:0000269 | Experimental evidence used in a manual assertion |
| ECO:0000250 | Sequence-similarity evidence used in a manual assertion (UniProt shows it as "By similarity") |
| ECO:0000255 | Match to a sequence model, used in a manual assertion |
| ECO:0000256 | Match to a sequence model, used in an automatic assertion |
| IEA:… (GO) | Inferred from electronic annotation, with the source named after the colon |

## Missing fields

| Gene ID | Missing |
|---|---|
| RPSAP15, TUBBP1 | All UniProt fields; no protein record is cross-referenced |
| AC069288.1 | Function, subcellular location and reviewed record. Name and gene name are submitted (unreviewed) values. |
| MATR3, ENSG00000280987 | Reviewed record. Function, location and GO are automatic annotation only. |
| INHBE, SFT2D2, STEAP2, EOGT, MATR3 ENSG00000015479 | None of the requested fields |

Two notes on the reviewed records:
- **Evidence type varies.** Some fields rest on similarity or sequence-model evidence rather than experiments:
  - SFT2D2 function: ECO:0000250 (similarity); location: ECO:0000255 (sequence model).
  - STEAP2 function: ECO:0000250 (similarity).
  - EOGT location: ECO:0000255 (sequence model).
- **Reference support:** only 7 PubMed references are cited directly as evidence for a function or location statement. SFT2D2 has none.

## Outputs

| File | Content |
|---|---|
| `results/tables/phase8_raw/uniprot_2026_04/entries/*.json`, `*_isoforms.fasta` | 14 UniProtKB entries and their isoform sequences, exactly as received |
| `results/tables/phase8_raw/uniprot_2026_04/searches/*.json` | UniProt searches for each gene ID's Ensembl cross-reference, exactly as received |
| `results/tables/phase8_raw/uniprot_2026_04/uniparc/*.json` | UniParc records for the canonical proteins of the two genes without a reviewed record |
| `results/tables/phase8_raw/uniprot_2026_04/request_log.csv` | Every request: URL, time, status, release header, size, SHA-256 |
| `results/tables/phase8_raw/uniprot_2026_04/ensembl97_pep_source.csv`, `ensembl97_pep_CHECKSUMS` | Provenance of the Ensembl 97 protein FASTA |
| `results/tables/phase8_raw/uniprot_2026_04/uniparc_evidence.csv` | Active human UniProtKB records with sequences identical to those canonical proteins |
| `results/tables/phase8_uniprot_mapping.csv` | Every gene–accession pair considered (90), with evidence, decision and reason |
| `results/tables/phase8_functional_annotation.csv` | One row per gene ID: record used, review status, name, function, location, evidence codes, GO counts, entry version, missing fields |
| `results/tables/phase8_functional_go_terms.csv` | 46 GO terms with aspect and evidence |
| `results/tables/phase8_functional_references.csv` | 7 PubMed references cited as evidence for function or location |
| `results/tables/phase8_functional_validation.csv` | 18 validation checks, all PASS |

The Ensembl 97 protein FASTA (14 MB) is cached outside the repository, in `%LOCALAPPDATA%\bioinformatics-portfolio-cache\ensembl97_pep`.
