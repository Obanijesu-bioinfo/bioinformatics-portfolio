# Projects

This section contains the practical bioinformatics projects in this portfolio. Each project has its own folder and README, which set out the research question, data sources, methods, results and limitations, together with the scripts, tables and figures needed to check the work.

The projects are educational analyses carried out to develop research skills. Where a project uses public data, the original experiments and data generation were performed by the original study's authors, and the project documents a reanalysis of their published data.

## Featured project

### RNA-seq analysis of 5-OXO-ETE and Calcitriol treatment effects in human stem cell models (GSE213410)

**[Read the full case study →](rnaseq-gse213410-case-study/README.md)**

An independent reanalysis of a publicly available RNA-seq dataset, investigating which gene-expression changes are associated with two treatments in two human stem cell models.

| | |
|---|---|
| **Project type** | Independent educational reanalysis of public data |
| **Dataset** | [GSE213410](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213410) (NCBI GEO): 12 samples, transcript-level read-count estimates supplied by the data submitters |
| **Organism** | *Homo sapiens* |
| **Comparison C1** | H9 human embryonic stem cells: 5-OXO-ETE vs ethanol control (3 replicates per group) |
| **Comparison C2** | Trisomy 21 induced pluripotent stem cells: Calcitriol vs DMSO control (3 replicates per group) |
| **Main tools** | R (DESeq2, clusterProfiler, ReactomePA, fgsea); Python (pandas, Matplotlib, scikit-learn); Git |
| **Reference resources** | Ensembl release 97 gene models; UniProtKB release 2026_04 protein annotation |

**Research focus.** The case study asks two questions:

- What gene-expression changes are associated with 5-OXO-ETE treatment in H9 human embryonic stem cells?
- What gene-expression changes are associated with Calcitriol treatment in trisomy 21 induced pluripotent stem cells?

The two comparisons use different cell models, treatments and controls, so they were analysed independently and their results are never combined or compared.

**Analytical stages.** The case study was completed in eight phases:

1. Project setup and documentation.
2. Dataset exploration and structural validation.
3. RNA-seq quality control and exploratory data analysis, including investigation of an outlier control sample (A3).
4. Sample metadata validation and experimental-design assessment.
5. Transcript-to-gene aggregation and gene-level quality control.
6. Differential expression analysis with DESeq2, with an exploratory sensitivity analysis excluding sample A3.
7. Functional enrichment analysis using GO, KEGG and Reactome, alongside gene set enrichment analysis (GSEA).
8. Gene annotation, functional annotation and cautious biological interpretation of eight selected genes.

**Methodological emphasis.**

- **Scope stated clearly.** The analysis starts from the count estimates provided with the dataset. No raw sequencing reads were aligned or quantified in this project.
- **Quality control and transparency.** Sample A3 was a quality-control outlier whose cause could not be determined. It was kept in the primary analysis, and its influence was assessed in a separate, explicitly exploratory analysis.
- **Provenance and versioning.** Annotation sources are version-fixed (Ensembl release 97, UniProtKB 2026_04), with retrieval and validation records.
- **Reproducibility.** The scripts for every phase are in the repository, with result tables from Phase 3 onwards and validation reports from Phase 4 onwards. A validation script checks most numerical claims in the case-study README against the machine-readable results.
- **Cautious interpretation.** Results are reported as statistical associations and treated as hypothesis-generating. RNA-level changes are not presented as evidence of protein or functional changes, and the study's limitations are documented.

**Where to start.**

- [Key findings](rnaseq-gse213410-case-study/README.md#6-key-findings)
- [Challenges and limitations](rnaseq-gse213410-case-study/README.md#7-challenges-and-limitations)
- [Reproducibility and evidence map](rnaseq-gse213410-case-study/README.md#9-reproducibility-and-evidence-map)
- [Visual results](rnaseq-gse213410-case-study/README.md#10-visual-results)

## Future projects

Further bioinformatics projects will be added to this section as they are completed. Each will follow the same structure: a dedicated folder with a README describing the research question, data sources, methods, results and limitations, and the code and outputs needed to reproduce the analysis.

| Project | Focus | Status |
|---|---|---|
| [GSE213410 RNA-seq case study](rnaseq-gse213410-case-study/README.md) | Bulk RNA-seq reanalysis: differential expression, functional enrichment and gene annotation | Complete (eight phases; final scientific review 2026-11-27) |
