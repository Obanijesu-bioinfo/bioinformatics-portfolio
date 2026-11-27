# Phase 8: Biological interpretation of the target genes

**Status:** Step 4 of Phase 8 (interpretation). Generated 2026-11-26; revised 2026-11-26 (see Section 7).
**Table:** `results/tables/phase8_gene_interpretation.csv`, with one row for each of the nine gene records.
**Validation:** see `results/tables/phase8_gene_interpretation_validation.csv`.
**Script:** `data/raw/notebooks/08c_biological_interpretation.py`

## 1. What this report covers

The report covers the eight target genes chosen from the differential-expression results. Four come from comparison **C1** (H9 embryonic stem cells, 5-OXO-ETE vs ethanol) and four from **C2** (trisomy 21 induced pluripotent stem cells, Calcitriol vs DMSO). MATR3 has two Ensembl gene IDs, and they are reported as separate records, which gives nine records in total.

C1 and C2 are separate experiments, with different cells, treatments and controls. They are never combined, and their fold changes are never compared with each other.

Each gene section separates four things:

- **What we know:** official annotation facts from Ensembl release 97 and UniProtKB release 2026_04. The source files are listed in Section 5.
- **What the results show:** the Phase 6 differential-expression result for that gene's own comparison.
- **What it may suggest:** a cautious hypothesis, clearly labelled.
- **What we cannot conclude.**

In the table, columns starting with `fact_` are copied or computed from saved outputs. Columns starting with `interp_` are interpretation, and columns starting with `conf_` are confidence levels.

**A change in RNA level is not evidence of a change in protein, activity or function.** Every interpretation below is a hypothesis that would need further experiments.

## 2. How to read the numbers and confidence levels

- **Fold change.** "log2FC" is DESeq2's unshrunken (MLE) log2 fold change, treated vs control. A log2FC of 1 means about 2 times higher; −1 means about 2 times lower.
- **Adjusted p-value.** The Benjamini–Hochberg adjusted p-value. All nine records are significant (adjusted p < 0.05) in their own comparison.
- **Normalised counts.** DESeq2 median-of-ratios values. They are given so you can see how consistent the replicates are.
- **95% interval.** An approximate interval for the fold change (log2FC ± 1.96 × standard error), computed from the saved Phase 6 values. It shows how exactly the size of each change is known.

Three separate confidence ratings are used. Each is assigned by written rules in the script.

**Expression-signal confidence** (`conf_expression_signal`). This rates how strongly and consistently the data support a real, reproducible change in the stated direction. It refers to the strength and consistency of the signal, **not certainty about the exact size of the fold change.** It counts these caveats:

1. Adjusted p ≥ 0.001.
2. |log2FC| < 1 (a modest change).
3. For C1 only: not significant, or a different direction, without A3.
4. The gene-level estimate may be ambiguous because the gene model is a readthrough model or overlaps another target record.
5. The gene is a pseudogene. Reads from related parent or paralogous sequences may contribute to its gene-level estimate. Sequence similarity was not assessed, so this is a possible concern, not evidence that reads were misassigned.

The rating is **High** with no caveats, **Moderate** with one, and **Low** with two or more.

**Effect-size precision** (`conf_effect_size_precision`). This rates how exactly the size of the change is known. It counts these issues:

1. Standard error ≥ 0.5 log2 units.
2. Mean normalised count below 50 in either group, which means a very low baseline.
3. The shrunken (apeglm) estimate collapses towards zero. Phase 6 defines this as "the test detects a difference, but its size is uncertain".
4. For C1 only: the estimate shifts by ≥ 0.5 log2 units when A3 is excluded.

The rating is **Low** with any issue. Without issues, it is **Good** when the standard error is < 0.25 and **Moderate** otherwise. This rating is purely statistical. It does not address whether a gene-level estimate is ambiguous, which the signal rating covers.

**Biological-interpretation confidence** (`conf_biological_interpretation`). This rates how much the expression change can tell us about biology:

- **Very low:** pseudogenes, and genes with no known function.
- **Low:** unreviewed or automatic-only annotation; function known only by similarity; or an expression signal with caveats.
- **Moderate:** a reviewed record with experimentally supported function **and** High expression-signal confidence.

No record is rated High, because an expression change on its own cannot establish what a gene does in these cells.

**Evidence codes** (UniProt):

| Code | Meaning |
|---|---|
| ECO:0000269 | Experimental |
| ECO:0000250 | By similarity to a related protein |
| ECO:0000255 | Predicted from a sequence model |
| ECO:0000256 | Automatic annotation |
| IEA (GO) | Electronic annotation |
| IDA, IPI, HDA (GO) | Experimental |
| ISS (GO) | Inferred from sequence similarity |
| IBA (GO) | Inferred from related genes (phylogenetic) |

### How the counts were produced

These details matter for genes whose models overlap.

- **StringTie (by the data submitters).** The data submitters aligned reads with HISAT2 and estimated **transcript-level** coverage with StringTie 1.3.3b. GEO provides only transcript-level counts.
- **prepDE.py.** As established in Phase 6, the integer counts follow StringTie's `prepDE.py` formula: coverage × transcript length / read length, rounded up, for each transcript. `prepDE.py` derives counts from coverage. It does not decide which gene a read belongs to.
- **Phase 5 gene sums.** Phase 5 summed transcript counts into Ensembl 97 genes. Each transcript belongs to exactly one gene (Phase 5 found 0 transcripts linked to more than one gene), so this step is a plain sum. It allocates nothing between genes.
- **Where ambiguity can arise.** When transcripts from two overlapping gene models are compatible with the same reads, any ambiguity lies in StringTie's transcript-level estimates. How StringTie handled those reads was not examined in this project. The separate gene-level estimates for overlapping models therefore **may not be independent**.

## 3. Summary

| Comparison | Gene | Ensembl gene ID | log2FC | Adjusted p | Direction | UniProt | Expression signal | Effect-size precision | Interpretation confidence |
|---|---|---|---|---|---|---|---|---|---|
| C1 | RPSAP15 | ENSG00000237506.3 | 1.39 | 1.73e-10 | higher | none (pseudogene) | Moderate | Good | Very low |
| C1 | TUBBP1 | ENSG00000127589.4 | 5.09 | 2.75e-08 | higher | none (pseudogene) | Moderate | Low | Very low |
| C1 | INHBE | ENSG00000139269.3 | −1.67 | 3.17e-08 | lower | P58166 (reviewed) | High | Good | Moderate |
| C1 | AC069288.1 | ENSG00000286192.1 | −3.48 | 8.17e-08 | lower | A0A3B3ITW8 (unreviewed; not used for function) | Moderate | Low | Very low |
| C2 | SFT2D2 | ENSG00000213064.10 | 2.61 | 1.52e-46 | higher | O95562 (reviewed) | High | Good | Low |
| C2 | MATR3 (primary record) | ENSG00000280987.4 | 1.77 | 2.64e-26 | higher | A0A0R4J2E8 (unreviewed, automatic only) | Moderate | Good | Low |
| C2 | MATR3 (overlapping model) | ENSG00000015479.18 | −0.59 | 3.26e-04 | lower | P43243 (reviewed) | Low | Good | Low |
| C2 | STEAP2 | ENSG00000157214.14 | 1.11 | 8.73e-23 | higher | Q8NFT2 (reviewed) | High | Good | Low |
| C2 | EOGT | ENSG00000163378.14 | −0.81 | 7.51e-15 | lower | Q5NDL2 (reviewed) | Moderate | Good | Low |

"Higher" and "lower" refer to the treated group compared with its own control.

### A3 sensitivity (C1 only)

Sample A3, one of the three C1 controls, was a quality-control outlier. Phase 6 therefore repeated C1 without it.

| Gene | log2FC with A3 | log2FC without A3 | Adjusted p without A3 | Still significant, same direction? | Shrunken estimate without A3 |
|---|---|---|---|---|---|
| RPSAP15 | 1.39 | 1.35 | 7.40e-30 | yes | stable |
| TUBBP1 | 5.09 | 6.40 | 4.13e-87 | yes | stable |
| INHBE | −1.67 | −1.61 | 2.24e-13 | yes | stable |
| AC069288.1 | −3.48 | −3.53 | 7.56e-08 | yes | **collapses to −0.07** |

All four C1 targets stay significant in the same direction when A3 is excluded. However, **the analysis without A3 is exploratory.** It has only two control replicates, and the decision to exclude A3 was informed by the QC results. It supports the robustness of the primary results but does not replace them.

## 4. Gene-by-gene interpretation

### 4.1 RPSAP15 (C1, ENSG00000237506.3)

**What we know.**
- Ensembl 97 annotates RPSAP15 as a **processed pseudogene**, "ribosomal protein SA pseudogene 15" (HGNC:31464), with one transcript.
- There is no UniProtKB record, so no protein function is known.

**What the results show.**
- Signal assigned to RPSAP15 was about **2.6 times higher** with 5-OXO-ETE (log2FC 1.39, adjusted p 1.73e-10; approximate 95% interval 2.1–3.4 times).
- The replicates agree closely: controls 683–769, treated 1,863–2,114 normalised counts.
- The result is unchanged without A3.

**What it may suggest.**
Something in the treated cells is producing more RNA that matches this pseudogene. It could be transcription of the pseudogene itself. It could also be reads from the related ribosomal protein SA gene (the gene its name refers to) being assigned here, because processed pseudogenes closely resemble the gene they came from.

**What we cannot conclude.**
That any protein is made, that the change has a biological effect, or where the reads truly come from. GEO provides only counts, so read-level assignment cannot be checked.

**Confidence.** Expression signal: Moderate. The change is strong and consistent, but this is a pseudogene, so reads from the related parent gene or similar sequences may contribute (sequence similarity was not assessed). Effect-size precision: Good. Biological interpretation: Very low, because this is a pseudogene with no protein function.

### 4.2 TUBBP1 (C1, ENSG00000127589.4)

**What we know.**
- Ensembl 97 annotates TUBBP1 as a **transcribed processed pseudogene**, "tubulin beta pseudogene 1" (HGNC:12414), with two transcripts. Its canonical transcript is annotated as a long non-coding RNA.
- There is no UniProtKB record, and the retrieved records document no function for the RNA.

**What the results show.**
- Signal was **much higher** with 5-OXO-ETE (log2FC 5.09, adjusted p 2.75e-08).
- Two controls are very low (A1 14, A2 15; control mean 39) and all three treated samples are high (948–1,360).
- **The direction of the change is well supported.** The adjusted p-value is strong, all three treated samples are consistently high, and the change stays significant in the same direction without A3.
- **The exact size of the change is not well determined.** The point estimate is about 34-fold, but:
  - the standard error is large (0.73 log2 units), giving an approximate 95% interval of about **13- to 92-fold**;
  - the control baseline is very low;
  - A3 is higher than the other controls (88), so excluding it moves the estimate to log2FC 6.40 (about 84-fold).

  The change is best described as "large", not as a precise 34-fold increase.

**What it may suggest.**
The locus is annotated as transcribed, so the treatment may be associated with its transcription. Some reads may instead come from related tubulin beta sequences.

**What we cannot conclude.**
That a protein is made, that the RNA has a function, that the change matters biologically, or the exact size of the increase.

**Confidence.**
- **Expression signal: Moderate.** The direction and consistency of the change are well supported. The one caveat is that this is a pseudogene, so reads from related tubulin beta sequences may contribute (sequence similarity was not assessed).
- **Effect-size precision: Low,** because of the large standard error, the low control baseline and the shift without A3.
- **Biological interpretation: Very low** (pseudogene).

### 4.3 INHBE (C1, ENSG00000139269.3)

**What we know.**
- **Ensembl 97:** protein coding, "inhibin subunit beta E" (HGNC:24029).
- **Function** (reviewed UniProt record P58166): INHBE belongs to the inhibin/activin family of secreted signalling proteins. UniProt describes activin E (made of two INHBE subunits) as secreted by the liver and involved in lipid metabolism and energy homeostasis. This is experimental evidence (PubMed:38533769, "Activin E is a transforming growth factor beta ligand that signals specifically through activin receptor-like kinase 7"). Its location, "secreted", is by similarity.
- **GO terms:**
  - cell-surface receptor serine/threonine kinase signalling (IBA);
  - cytokine activity (ISS);
  - negative regulation of adipose tissue development (ISS).
- **Ensembl 97 Reactome links:** Glycoprotein hormones; Peptide hormone biosynthesis; Peptide hormone metabolism.

**What the results show.**
- INHBE RNA was about **3.2 times lower** with 5-OXO-ETE (log2FC −1.67, adjusted p 3.17e-08; approximate 95% interval 2.3–4.4 times lower): controls 1,575–1,954, treated 504–675.
- The result is the same without A3 (log2FC −1.61).

**What it may suggest.**
If less RNA led to less INHBE protein, these cells might produce less activin E. This is a hypothesis. Protein and secretion were not measured, and the documented functions come mainly from liver, so their relevance to embryonic stem cells is unknown.

**What we cannot conclude.**
That INHBE protein or activin E signalling changed, or that INHBE mediates any effect of 5-OXO-ETE.

**Confidence.** Expression signal: High. Effect-size precision: Good. Biological interpretation: Moderate, the highest level in this report, because the function has experimental support. It is still only a hypothesis about these cells.

### 4.4 AC069288.1 (C1, ENSG00000286192.1)

**What we know.**
- **Ensembl 97:**
  - The gene's description is only "Novel protein", and its name is a placeholder based on a sequencing clone. There is no HGNC or NCBI Gene record.
  - Its single transcript (ENST00000651235.1) carries Ensembl's official **readthrough transcript** flag. Its biotype is **nonsense-mediated decay**, a class of transcripts expected to be degraded, so protein production cannot be assumed.
  - The gene spans the MAD1L1 locus and extends into MRM2 on the same strand. This is a coordinate fact noted in Step 1.
- **UniProt:**
  - The only record cross-referenced to it is the **unreviewed** A0A3B3ITW8. That record carries another gene's name, has no function or location, and has one **electronic-only** GO term.
  - That record's annotation is **not** used for AC069288.1.

**What the results show.**
- Signal assigned to this gene model was about **11 times lower** with 5-OXO-ETE (log2FC −3.48, adjusted p 8.17e-08; approximate 95% interval 5.5–23 times lower).
- The controls vary about three-fold (A1 195, A2 632, A3 376), while the treated samples are low (28–41).
- Phase 6 also found **MRM2 significantly lower** in C1 (log2FC −1.29, adjusted p 3.10e-05), while **MAD1L1 was not significant** (log2FC −0.67, adjusted p 0.91).
- Without A3, the unshrunken change is similar (−3.53), but the shrunken estimate collapses towards zero (−0.07). The size of the effect is therefore uncertain.

**What it may suggest.**
The decrease may be related to the MRM2 decrease in the same samples.
- The two gene models share genomic sequence, so their transcripts can be compatible with some of the same reads.
- Their gene-level counts are separate sums of StringTie transcript estimates (Section 2), so the two estimates **may not be independent**.
- Gene-level counts cannot tell whether a distinct readthrough RNA changed on its own.

**What we cannot conclude.**
Any function, that a protein is made, the exact size of the change, or that this RNA changed independently of MRM2.

**Confidence.**
- **Expression signal: Moderate.** The change is significant, large, and robust in direction without A3. The one caveat is that the gene-level estimate may be ambiguous, because this is a readthrough model overlapping other genes.
- **Effect-size precision: Low:** standard error 0.52, a very low treated mean (36), and a shrunken estimate that collapses without A3.
- **Biological interpretation: Very low,** because no function is known.

(Before this revision, the single expression rating was Low. The collapse of the shrunken estimate now counts as an effect-size issue rather than a signal issue; see Section 7.)

### 4.5 SFT2D2 (C2, ENSG00000213064.10)

**What we know.**
- **Ensembl 97:** protein coding, "SFT2 domain containing 2" (HGNC:25140).
- **Reviewed UniProt record O95562** ("Vesicle transport protein SFT2B"). Its only function statement is that it "may be involved in fusion of retrograde transport vesicles derived from an endocytic compartment with the Golgi complex".
- **The evidence for that function is by similarity only** (ECO:0000250), meaning it is inferred from related proteins. It has not been shown by experiments on this protein.
- Its membrane location comes from a sequence model (ECO:0000255).
- **GO terms:**
  - vesicle-mediated transport: electronic only;
  - extracellular exosomes: the only experimental-type term, a high-throughput detection (HDA, PubMed:19056867).

**What the results show.**
- SFT2D2 was the **most significant C2 gene**: about **6 times higher** with Calcitriol (log2FC 2.61, adjusted p 1.52e-46; approximate 95% interval 4.8–7.7 times).
- The replicates agree closely: controls 699–819, treated 3,529–5,381.

**What it may suggest.**
If the similarity-based function holds in these cells, higher SFT2D2 might be connected with changes in vesicle trafficking. This is speculative, because the function itself is not experimentally established.

**What we cannot conclude.**
What SFT2D2 does in these cells, whether Calcitriol affects it directly or indirectly, or whether vesicle transport changed.

**Confidence.** Expression signal: High. Effect-size precision: Good. Biological interpretation: Low, because the function is not experimentally established.

### 4.6 and 4.7 MATR3: two overlapping gene models (C2)

**What we know.**

Ensembl release 97 contains two gene records named MATR3 on chromosome 5. They share 14 identical exons, and most of the smaller model's exon sequence lies inside the larger one's.

**ENSG00000280987.4: primary record (the selected C2 result)**

- **Ensembl 97:**
  - The description "matrin 3" comes from NCBI Gene 9782; no HGNC record is attached.
  - **10 of its 11 transcripts are officially flagged as readthrough transcripts.**
  - Havana curation terms: "ncRNA host", "overlapping locus".
- **UniProt:**
  - Only unreviewed records are linked. The main one, A0A0R4J2E8, carries function text and GO terms that are **automatic only** (ECO:0000256 / IEA).
  - Its protein sequence is identical to the reviewed MATR3 record P43243. However, no Ensembl 97 cross-reference links this gene to P43243, so the reviewed annotation is **not transferred** (see `phase8_uniprot_mapping.csv`).

**ENSG00000015479.18: overlapping gene model**

- **Ensembl 97:** "matrin 3" (HGNC:6912). It has 28 transcripts, none flagged as readthrough.
- **Reviewed UniProt record P43243 (Matrin-3),** with experimentally supported functions:
  - a nuclear matrix protein that, with the SFPQ-NONO complex, helps keep defective RNAs in the nucleus (PubMed:11525732);
  - a role in innate immune responses to DNA viruses (PubMed:28712728);
  - binding of m6A-modified mRNAs, including MYC (PubMed:32245947);
  - binding of specific miRNA hairpins (PubMed:28431233).
- Experimental GO terms include post-transcriptional regulation of gene expression (IDA, PubMed:21771347).

Both IDs carry the same OMIM links in Ensembl 97, MATR3 (164015) and amyotrophic lateral sclerosis 21 (606070). These disease links are recorded as annotation only and are not used to interpret these cells.

**What the results show.**

| Record | log2FC | Adjusted p | Approx. 95% interval | Controls (A7–A9) | Calcitriol (A10–A12) |
|---|---|---|---|---|---|
| ENSG00000280987 (primary) | **+1.77** (about 3.4× higher) | 2.64e-26 | 2.8–4.2× higher | 3,860–5,538 | 14,095–17,267 |
| ENSG00000015479 (overlapping) | **−0.59** (about 1.5× lower) | 3.26e-04 | 1.3–1.8× lower | 14,418–15,488 | 8,716–11,288 |

As a descriptive check only, not a statistical test: adding the two models' normalised counts together gives 18,441–19,985 in controls and 24,685–26,918 with Calcitriol.

**What it may suggest.**
- **Each record's count is a separate sum of StringTie estimates for its own transcripts** (Section 2). The two models share exons, so their transcripts can be compatible with the same reads. The two gene-level estimates therefore **may not be independent**.
- **The opposite changes could have three explanations,** and these data cannot distinguish between them:
  - a real change in the relative abundance of different transcript structures at the MATR3 locus;
  - ambiguity in estimating overlapping transcripts;
  - both.
- **The combined total.** The descriptive sum suggests that total MATR3-locus signal was somewhat higher with Calcitriol, but this was not tested.
- **What would be needed:** a transcript-level analysis of the MATR3 locus. It was not done here, and it would require a new approved step.

**What we cannot conclude.**
- That MATR3 protein increased or decreased, or that any MATR3 function changed.
- That the two changes are independent biological effects.
- That the reviewed MATR3 functions apply to the primary record: its annotation is unreviewed and automatic.

**Confidence.**

| Record | Expression signal | Effect-size precision | Biological interpretation |
|---|---|---|---|
| ENSG00000280987 | Moderate: gene-level estimate may be ambiguous (overlapping/readthrough models) | Good | Low: unreviewed, automatic annotation only |
| ENSG00000015479 | Low: modest effect; gene-level estimate may be ambiguous | Good | Low: expression signal has caveats |

Effect-size precision is Good for both records because their standard errors are small. That rating describes only the statistical precision of each gene-level estimate. It does not resolve the possible ambiguity between the two overlapping models.

### 4.8 STEAP2 (C2, ENSG00000157214.14)

**What we know.**
- **Ensembl 97:** protein coding, "STEAP2 metalloreductase" (HGNC:17885).
- **Reviewed UniProt record Q8NFT2.** Its described enzyme function, reducing iron (Fe³⁺) chelates and copper (Cu²⁺) using NADPH, is **by similarity** (ECO:0000250).
- **Location in endosome and cell membranes** includes experimental evidence (PubMed:12429817).
- **GO terms:**
  - Experimental (PubMed:12095985): endocytosis; Golgi-to-plasma-membrane transport; regulated exocytosis; response to hormone.
  - Inferred from related genes (IBA): copper ion import; ferric-chelate reductase activity.
- **Ensembl 97 Reactome links:** Iron uptake and transport; Transferrin endocytosis and recycling.

**What the results show.**
- STEAP2 was about **2.2 times higher** with Calcitriol (log2FC 1.11, adjusted p 8.73e-23; approximate 95% interval 1.9–2.5 times).
- The replicates agree closely: controls 639–718, treated 1,392–1,523.

**What it may suggest.**
If the similarity-based enzyme function holds, higher STEAP2 RNA might be associated with more capacity to reduce metal ions during endosomal uptake. This was not tested.

**What we cannot conclude.**
That iron or copper handling changed, that STEAP2 protein or enzyme activity increased, or any mechanism.

**Confidence.** Expression signal: High. Effect-size precision: Good. Biological interpretation: Low, because the enzyme function is by similarity.

### 4.9 EOGT (C2, ENSG00000163378.14)

**What we know.**
- **Ensembl 97:** protein coding, "EGF domain specific O-linked N-acetylglucosamine transferase" (HGNC:28526).
- **Function** (reviewed UniProt record Q5NDL2), with experimental evidence (ECO:0000269, PubMed:23671640): it adds a single sugar, O-GlcNAc, to serine or threonine residues in EGF-like domains of extracellular proteins.
- **Location:** endoplasmic reticulum lumen, predicted from a sequence model.
- **OMIM link in Ensembl 97:** Adams–Oliver syndrome 4 (615297). It is recorded as annotation only and not used for interpretation.

**What the results show.**
- EOGT was modestly lower, about **1.8 times lower**, with Calcitriol (log2FC −0.81, adjusted p 7.51e-15; approximate 95% interval 1.6–2.0 times lower).
- The replicates are very consistent: controls 1,716–1,737, treated 908–1,024.

**What it may suggest.**
Less EOGT RNA might mean less O-GlcNAc modification of EGF-domain proteins, if protein levels follow. Neither EOGT protein nor the modification was measured.

**What we cannot conclude.**
That glycosylation of any protein changed, or that any downstream signalling was affected.

**Confidence.** Expression signal: Moderate, because the change is modest (|log2FC| < 1). Effect-size precision: Good. Biological interpretation: Low, because the expression signal has a caveat, despite the experimentally supported function.

## 5. Sources (all saved locally)

| Source | Files |
|---|---|
| Phase 6 differential expression | `results/tables/C1_differential_expression.csv`, `C1_differential_expression_without_A3.csv`, `C2_differential_expression.csv`, `C1_A3_sensitivity_comparison.csv`; counting workflow described in `results/PHASE6_REPORT.md` (Section 1) and `results/PHASE5_REPORT.md` |
| Sample groups | `data/metadata/sample_metadata.csv` |
| Ensembl release 97 annotation (Step 2) | `results/tables/phase8_target_gene_annotation.csv`, `results/tables/phase8_raw/ensembl97_*.csv`. Provenance: `results/PHASE8_ANNOTATION_PROVENANCE.md` |
| UniProtKB 2026_04 (Step 3) | `results/tables/phase8_functional_annotation.csv`, `phase8_functional_go_terms.csv`, `phase8_functional_references.csv`, `phase8_uniprot_mapping.csv`; raw entries in `results/tables/phase8_raw/uniprot_2026_04/entries/`. Provenance: `results/PHASE8_FUNCTIONAL_ANNOTATION_PROVENANCE.md` |

Every PubMed ID cited here is cited by the saved UniProt records, either as evidence for a function or location statement or as GO evidence. No other literature was used.

## 6. Limitations that apply to all records

- **RNA is not protein.** Differential expression shows a change in RNA signal assigned to a gene model. It does not show a change in protein, activity or function, and it does not identify a mechanism.
- **Study design.** Three replicates per group, a single treatment time point, and one cell model per comparison.
- **How counts were produced** (Section 2). Counts come from StringTie 1.3.3b transcript-level coverage estimates, converted to integers by the `prepDE.py` formula, then summed by gene in Phase 5. Neither `prepDE.py` nor Phase 5 allocates reads between genes. Where overlapping gene models share sequence, any ambiguity lies in StringTie's transcript-level estimation, which was not examined. The separate gene-level estimates may then not be independent. This matters most for AC069288.1 and the two MATR3 records.
- **Signal is not size.** A well-supported direction of change does not mean the exact fold change is known. It can still be imprecise, as for TUBBP1. Effect-size precision and the 95% intervals show this separately.
- **Annotation is time-stamped.** Gene models come from Ensembl 97 (2019). Function annotation comes from UniProtKB 2026_04. Records can change between releases, as the MATR3 cross-references show.
- **Evidence quality varies.** Several "known" functions rest on similarity, sequence models or automatic annotation rather than experiments on the protein. The table records the evidence class for every statement.
- **The C1 sensitivity analysis is exploratory.** Without A3 it has only two control replicates, and excluding A3 was informed by QC.
- **C1 and C2 are independent.** They use different cells, treatments and controls and must not be combined or compared.

## 7. Revision note (2026-11-26)

This was a targeted revision. No statistical value, gene ID, annotation record, UniProt evidence or A3 sensitivity result changed.

1. **Expression confidence was split into two ratings.**
   - The single "expression-result confidence" became **expression-signal confidence** (direction and consistency) plus **effect-size precision** (how exactly the fold change is known).
   - New table fields: `conf_effect_size_precision`, its issues, group means, and approximate 95% intervals.
   - TUBBP1 keeps High signal confidence and now has Low effect-size precision, with the reasons stated.
   - Shrinkage collapse moved from the signal rating to the precision rating, matching its Phase 6 definition. As a result, **AC069288.1's signal rating changed from Low to Moderate**, and its precision rating is Low. Its biological-interpretation rating is unchanged (Very low).
   - No other rating changed.
2. **The overlapping-gene wording was corrected.** The report no longer says that reads are allocated between overlapping gene models during quantification. It now describes the actual workflow (StringTie transcript estimates, `prepDE.py`-derived counts, Phase 5 plain sums) and states only that separate gene-level estimates for overlapping models may not be independent. This applies to AC069288.1/MRM2, both MATR3 records, and the limitations section.

## 8. Revision note (2026-11-27)

This was a targeted revision. No statistical value, gene ID, annotation record, UniProt evidence or A3 sensitivity result changed.

1. **Pseudogene status now counts as an expression-signal caveat** (Section 2, caveat 5). Reads from related parent or paralogous sequences may contribute to a pseudogene's gene-level estimate. Sequence similarity was not assessed, so this is a possible concern, not evidence of misassignment.
2. **Effect on the ratings:** RPSAP15 and TUBBP1, the only pseudogenes among the nine records, change from High to **Moderate** expression-signal confidence. Their effect-size precision (Good and Low) and biological-interpretation ratings (Very low) are unchanged. No other record's rating changed.
3. This supersedes the statement in Section 7 that TUBBP1 keeps High signal confidence.
