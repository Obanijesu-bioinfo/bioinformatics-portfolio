"""
Phase 8, Step 4: biological interpretation table for the nine target gene records.

Inputs (read only)
  Phase 6:  results/tables/C1_differential_expression.csv
            results/tables/C1_differential_expression_without_A3.csv
            results/tables/C2_differential_expression.csv
            results/tables/C1_A3_sensitivity_comparison.csv
  Step 2:   results/tables/phase8_target_gene_annotation.csv, results/tables/phase8_raw/ensembl97_*.csv
  Step 3:   results/tables/phase8_functional_annotation.csv, phase8_functional_go_terms.csv,
            phase8_functional_references.csv, phase8_uniprot_mapping.csv

Output
  results/tables/phase8_gene_interpretation.csv           (one row per gene record)
  results/tables/phase8_gene_interpretation_validation.csv

Column groups
  fact_*    copied or computed from the saved Phase 6 / Ensembl 97 / UniProt 2026_04 outputs
  conf_*    confidence levels assigned by the explicit rules in assign_confidence()
  interp_*  written interpretation, hypotheses and limitations (INTERPRETATION below)

Revision 2026-11-26 (targeted scientific revision)
  1. Confidence is split into (a) expression-signal confidence: how strongly and
     consistently the data support a real, reproducible change in the stated direction;
     and (b) effect-size precision: how exactly the size of the fold change is known.
     Shrinkage collapse moved from (a) to (b), because Phase 6 defines it as "the test
     detects a difference, but its size is uncertain".
  2. Wording about overlapping gene models no longer says that reads are divided
     between genes. It now describes the actual workflow (StringTie transcript
     estimates -> prepDE.py-derived counts -> Phase 5 plain sum per gene) and the
     resulting potential ambiguity. The report's text is also checked here.

Rules followed in the interpretation text: differential expression is not treated as
proof of function or mechanism; pseudogene expression is not treated as protein
production; MRM2 annotation is not assigned to AC069288.1; the two MATR3 gene models are
kept separate; C1 and C2 are never combined or compared.

The script refuses to overwrite existing outputs. Run from anywhere:
    python data/raw/notebooks/08c_biological_interpretation.py
"""

import csv
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
T = REPO_ROOT / "results/tables"
OUT = T / "phase8_gene_interpretation.csv"
OUT_VAL = T / "phase8_gene_interpretation_validation.csv"

COMPARISONS = {
    "C1": ("H9 human embryonic stem cells: 5-OXO-ETE vs ethanol control, 3 vs 3 replicates "
           "(primary analysis, A3 included)", "C1_differential_expression.csv"),
    "C2": ("Trisomy 21 induced pluripotent stem cells: Calcitriol vs DMSO control, 3 vs 3 replicates",
           "C2_differential_expression.csv"),
}
C1_EXPLORATORY = ("The C1 analysis without A3 is exploratory: it has only two control replicates "
                  "(A1, A2), and excluding A3 was a decision informed by the QC results.")
PSEUDO_NOTE = ("Processed pseudogenes resemble the gene they derive from, so some reads assigned here may "
               "come from similar sequences elsewhere; GEO provides counts only, so read assignment cannot be checked.")
WORKFLOW_NOTE = ("Counting workflow: the data submitters estimated transcript-level coverage with StringTie 1.3.3b; "
                 "the integer counts follow prepDE.py's formula (coverage x transcript length / read length, rounded up; "
                 "Phase 6); Phase 5 then summed transcript counts to Ensembl 97 genes, each transcript belonging to exactly "
                 "one gene. Neither prepDE.py nor Phase 5 allocates reads between genes. Where transcripts of overlapping "
                 "gene models are compatible with the same reads, any ambiguity lies in StringTie's transcript-level "
                 "estimates, which were not examined here, so the separate gene-level estimates may not be independent.")
OVERLAP_BANNED = re.compile(r"divided between|redistribut|split between|shared out|are divided", re.I)

# ---------------------------------------------------------------------------
# Written interpretation, keyed by Ensembl gene ID. Facts are not repeated here
# except where needed to explain the reasoning.
# ---------------------------------------------------------------------------
INTERPRETATION = {
    "ENSG00000237506": {  # RPSAP15
        "interp_known_function": "No protein function. RPSAP15 is annotated as a processed pseudogene "
            "(\"ribosomal protein SA pseudogene 15\") with no UniProtKB record.",
        "interp_process_pathway": "None supported by the retrieved annotation.",
        "interp_expression_change": "More RNA signal was assigned to this pseudogene in 5-OXO-ETE-treated cells "
            "than in ethanol controls, consistently across all six samples. The difference is a reproducible "
            "change in reads assigned to this annotation; it is not evidence of a protein or of a function.",
        "interp_hypotheses": "Two explanations are possible and cannot be separated with these data: (1) "
            "increased transcription of the pseudogene locus itself; (2) reads from the related ribosomal "
            "protein SA gene (the gene its name refers to) assigned to the pseudogene.",
        "interp_cannot_conclude": "That a protein is made; that the change has a biological effect; "
            "where the reads truly originate.",
        "interp_limitations": PSEUDO_NOTE + " " + C1_EXPLORATORY,
    },
    "ENSG00000127589": {  # TUBBP1
        "interp_known_function": "No protein function. TUBBP1 is annotated as a transcribed processed pseudogene "
            "(\"tubulin beta pseudogene 1\") with no UniProtKB record. Its canonical transcript is annotated as "
            "a long non-coding RNA (lncRNA), which has no documented function in the retrieved records.",
        "interp_process_pathway": "None supported by the retrieved annotation.",
        "interp_expression_change": "Signal at this locus is very low in two controls (A1, A2) and much "
            "higher in all treated samples. The direction of the change is well supported: a strong adjusted "
            "p-value, consistently high values in all three treated samples, and the same direction (still "
            "significant) without A3. The exact size of the change is not well determined: the control baseline is "
            "very low, the standard error is large, and the estimate moves substantially when A3 is excluded (see "
            "fact_fold_change_95pct_wald_interval). The fold change should be read as 'large' rather than as a "
            "precise number.",
        "interp_hypotheses": "Treatment may be associated with transcription of this locus, which is annotated "
            "as transcribed. Alternatively, some reads may derive from related tubulin beta sequences. Neither "
            "can be tested with gene-level counts.",
        "interp_cannot_conclude": "That a protein is made; that the lncRNA transcript has a function; that the "
            "change has a biological effect.",
        "interp_limitations": PSEUDO_NOTE + " Effect-size precision is low: A3 has a higher value than the "
            "other two controls, so excluding A3 increases the estimated fold change; High expression-signal "
            "confidence refers to the direction and consistency of the change, not its exact size. " + C1_EXPLORATORY,
    },
    "ENSG00000139269": {  # INHBE
        "interp_known_function": "Reviewed UniProt record P58166: INHBE is a subunit of the inhibin/activin family "
            "of secreted signalling proteins. UniProt describes activin E (an INHBE homodimer) as secreted by the "
            "liver and involved in lipid metabolism and energy homeostasis (experimental evidence, PubMed:38533769). "
            "Location: secreted (by similarity).",
        "interp_process_pathway": "GO: cell surface receptor protein serine/threonine kinase signalling pathway "
            "(IBA); cytokine activity (ISS); negative regulation of adipose tissue development (ISS). Ensembl 97 "
            "Reactome links: Glycoprotein hormones; Peptide hormone biosynthesis; Peptide hormone metabolism.",
        "interp_expression_change": "INHBE mRNA was lower in 5-OXO-ETE-treated cells, consistently across "
            "replicates and with or without A3.",
        "interp_hypotheses": "If lower mRNA led to less INHBE protein, these cells could produce less activin E "
            "ligand. This is a hypothesis only: protein levels and secretion were not measured, and the documented "
            "functions come from other tissues (mainly liver), so their relevance to embryonic stem cells is unknown.",
        "interp_cannot_conclude": "That INHBE protein or activin E signalling changed; that INHBE mediates any "
            "effect of 5-OXO-ETE.",
        "interp_limitations": "Expression change only; single time point; three replicates per group. "
            + C1_EXPLORATORY,
    },
    "ENSG00000286192": {  # AC069288.1
        "interp_known_function": "No known function. Ensembl 97 describes it only as \"Novel protein\" with a "
            "clone-based placeholder name. Its only transcript is officially flagged as a readthrough transcript and "
            "has the biotype nonsense-mediated decay, so protein production cannot be assumed. The only UniProt "
            "record cross-referenced to it (A0A3B3ITW8) is unreviewed, has no function or location, and carries "
            "another gene's name; that record's annotation is not assigned to this gene.",
        "interp_process_pathway": "None supported. The record's single GO term is electronic only and is not used "
            "for interpretation.",
        "interp_expression_change": "Signal assigned to this gene model was lower in treated cells. The model spans "
            "the MAD1L1 locus and extends into MRM2 on the same strand. MRM2 itself was also significantly lower in "
            "C1 (Phase 6), while MAD1L1 was not.",
        "interp_hypotheses": "The decrease may be related to the MRM2 decrease in the same samples. The two gene "
            "models share genomic sequence, so their transcripts can be compatible with some of the same reads, and "
            "the separate gene-level estimates may not be independent. Whether a distinct readthrough RNA changed "
            "cannot be separated from changes in MRM2 using gene-level counts.",
        "interp_cannot_conclude": "Any function; that a protein is made; that this readthrough RNA changed "
            "independently of MRM2.",
        "interp_limitations": "Gene model is a readthrough annotation overlapping two other genes. " + WORKFLOW_NOTE
            + " Control values vary about three-fold (A1 vs A2). Without A3 the MLE fold change is similar, but the "
            "shrunken estimate collapses towards zero, so the size of the effect is uncertain (low effect-size "
            "precision). " + C1_EXPLORATORY,
    },
    "ENSG00000213064": {  # SFT2D2
        "interp_known_function": "Reviewed UniProt record O95562 (\"Vesicle transport protein SFT2B\"). The only "
            "function statement, \"may be involved in fusion of retrograde transport vesicles derived from an "
            "endocytic compartment with the Golgi complex\", is supported by sequence similarity (ECO:0000250), not "
            "by experiments on this protein. Membrane location is from a sequence model (ECO:0000255).",
        "interp_process_pathway": "GO: vesicle-mediated transport (electronic, IEA:InterPro). The only "
            "experimentally derived GO term is a high-throughput detection in extracellular exosomes "
            "(HDA, PubMed:19056867).",
        "interp_expression_change": "SFT2D2 was the most significant C2 gene: about six-fold higher with "
            "Calcitriol, with tight agreement between replicates.",
        "interp_hypotheses": "If the similarity-based function applies in these cells, higher SFT2D2 could relate "
            "to changes in vesicle trafficking. This is speculative because the function itself has not been "
            "established experimentally in the retrieved record.",
        "interp_cannot_conclude": "SFT2D2's function in these cells; whether Calcitriol affects SFT2D2 directly "
            "or indirectly; any effect on vesicle transport.",
        "interp_limitations": "Function evidence is by similarity only. Expression change only; single time "
            "point; three replicates per group.",
    },
    "ENSG00000280987": {  # MATR3, primary record
        "interp_known_function": "No reviewed annotation for this gene model. Its unreviewed UniProt record "
            "A0A0R4J2E8 carries function text and GO terms that are automatic (ECO:0000256 / IEA) only. Its "
            "protein sequence is identical to the reviewed MATR3 record P43243, but that reviewed annotation is not "
            "transferred because no Ensembl 97 cross-reference links them (see phase8_uniprot_mapping.csv).",
        "interp_process_pathway": "None from reviewed evidence. Automatic GO terms only (RNA binding, zinc ion "
            "binding, nuclear matrix).",
        "interp_expression_change": "Signal assigned to this gene model was about 3.4-fold higher with Calcitriol. "
            "Ten of its eleven transcripts are officially flagged as readthrough transcripts, and the model overlaps "
            "the other MATR3 model (ENSG00000015479), which moved in the opposite direction.",
        "interp_hypotheses": "Each MATR3 record's count is a separate sum of StringTie estimates for its own "
            "transcripts. Because the two models share exons, their transcripts can be compatible with the same "
            "reads, so the two gene-level estimates may not be independent. The opposite changes could reflect a "
            "real change in the relative abundance of different transcript structures at the MATR3 locus, ambiguity "
            "in estimating overlapping transcripts, or both; these data cannot distinguish between them. A "
            "transcript-level look at the MATR3 locus (not done here) would be needed to explore this.",
        "interp_cannot_conclude": "That MATR3 protein increased; any function; that the increase is independent of "
            "the decrease in ENSG00000015479.",
        "interp_limitations": "Overlapping and readthrough gene models. " + WORKFLOW_NOTE + " Annotation automatic "
            "only; expression change only; single time point; three replicates per group.",
    },
    "ENSG00000015479": {  # MATR3, overlapping model
        "interp_known_function": "Reviewed UniProt record P43243 (Matrin-3). Experimentally supported function "
            "(ECO:0000269): a nuclear matrix protein involved, with the SFPQ-NONO heteromer, in nuclear retention of "
            "defective RNAs (PubMed:11525732); a role in DNA virus-mediated innate immune responses through the "
            "HDP-RNP complex (PubMed:28712728); binding of m6A-containing mRNAs including MYC (PubMed:32245947); "
            "binding of specific miRNA hairpins (PubMed:28431233).",
        "interp_process_pathway": "GO (experimental): activation of innate immune response (IDA, PubMed:28712728); "
            "post-transcriptional regulation of gene expression (IDA, PubMed:21771347); miRNA binding (IDA, "
            "PubMed:28431233).",
        "interp_expression_change": "Signal assigned to this gene model was modestly lower with Calcitriol "
            "(about 1.5-fold). It overlaps the primary MATR3 model, which increased.",
        "interp_hypotheses": "As for ENSG00000280987: the two records' gene-level estimates may not be independent "
            "because the models share exons, so the opposite directions cannot be read as two separate biological "
            "changes. The reviewed functions describe the MATR3 protein in general; they do not indicate what the "
            "change in this gene model means.",
        "interp_cannot_conclude": "That MATR3 protein decreased; that any MATR3 function changed; that this change "
            "is independent of the increase in ENSG00000280987.",
        "interp_limitations": "Overlapping gene models. " + WORKFLOW_NOTE + " Modest effect size (|log2FC| < 1); "
            "expression change only; single time point; three replicates per group.",
    },
    "ENSG00000157214": {  # STEAP2
        "interp_known_function": "Reviewed UniProt record Q8NFT2 (Metalloreductase STEAP2). Its described function, "
            "NADPH-dependent reduction of ferric (Fe3+) chelates and of Cu2+, is by similarity (ECO:0000250). "
            "Its location in endosome and cell membranes includes experimental evidence (PubMed:12429817).",
        "interp_process_pathway": "GO (experimental, PubMed:12095985): endocytosis; Golgi to plasma membrane transport; "
            "regulated exocytosis; response to hormone. GO (IBA): copper ion import; ferric-chelate reductase "
            "activity. Ensembl 97 Reactome links: Iron uptake and transport; Transferrin endocytosis and recycling.",
        "interp_expression_change": "STEAP2 was about two-fold higher with Calcitriol, consistently across replicates.",
        "interp_hypotheses": "Higher STEAP2 mRNA could be associated with greater capacity for metal-ion reduction "
            "linked to endosomal uptake, if the similarity-based enzymatic function applies. Not tested here.",
        "interp_cannot_conclude": "That iron or copper handling changed; that STEAP2 protein or enzyme activity "
            "increased; any mechanism.",
        "interp_limitations": "Enzymatic function is by similarity; expression change only; single time point; "
            "three replicates per group.",
    },
    "ENSG00000163378": {  # EOGT
        "interp_known_function": "Reviewed UniProt record Q5NDL2. Experimentally supported function (ECO:0000269, "
            "PubMed:23671640): transfers a single N-acetylglucosamine (O-GlcNAc) to serine or threonine residues in "
            "EGF-like domains of extracellular proteins. Location: endoplasmic reticulum lumen (sequence model).",
        "interp_process_pathway": "GO: protein O-linked glycosylation and protein O-acetylglucosaminyltransferase "
            "activity (ISS); endoplasmic reticulum lumen (IBA).",
        "interp_expression_change": "EOGT was modestly lower with Calcitriol (about 1.8-fold), consistently across "
            "replicates.",
        "interp_hypotheses": "Lower EOGT mRNA might reduce O-GlcNAc modification of EGF-domain proteins, if protein "
            "levels follow. Neither EOGT protein nor glycosylation was measured.",
        "interp_cannot_conclude": "That glycosylation of any protein changed; any downstream signalling effect.",
        "interp_limitations": "Modest effect size (|log2FC| < 1); expression change only; single time point; "
            "three replicates per group.",
    },
}


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt_p(x):
    return f"{float(x):.3g}"


def evidence_class(codes):
    """Summarise UniProt ECO codes for a statement."""
    s = {c.strip() for c in codes.split(";") if c.strip()}
    if not s:
        return "none"
    if "ECO:0000269" in s:
        return "experimental (with other evidence)" if len(s) > 1 else "experimental"
    if s <= {"ECO:0000256", "ECO:0000313", "ECO:0007829"}:
        return "automatic only"
    if s == {"ECO:0000250"}:
        return "by similarity only"
    if s == {"ECO:0000255"}:
        return "sequence model only"
    return "non-experimental (" + "; ".join(sorted(s)) + ")"


SE_IMPRECISE = 0.5        # log2 units
SE_GOOD = 0.25            # log2 units
LOW_BASELINE = 50         # mean normalised count in either group
A3_SHIFT = 0.5            # |log2FC with A3 - log2FC without A3|


def effect_size_precision(r):
    """How exactly the size of the fold change is known (separate from whether the change is real)."""
    issues = []
    se = float(r["fact_lfcSE_MLE"])
    if se >= SE_IMPRECISE:
        issues.append(f"standard error {se:.2f} >= {SE_IMPRECISE} log2 units")
    low = min(float(r["fact_control_mean_normalised"]), float(r["fact_treated_mean_normalised"]))
    if low < LOW_BASELINE:
        issues.append(f"low expression in one group (mean {low:.0f} < {LOW_BASELINE} normalised counts)")
    if r["fact_shrinkage_collapsed"] == "TRUE":
        issues.append("shrunken estimate collapses")
    if r["fact_comparison"] == "C1":
        if r["fact_A3_excluded_shrinkage_collapsed"] == "TRUE":
            issues.append("shrunken estimate collapses without A3")
        shift = abs(float(r["fact_log2FC_MLE"]) - float(r["fact_A3_excluded_log2FC_MLE"]))
        if shift >= A3_SHIFT:
            issues.append(f"estimate shifts by {shift:.2f} log2 units when A3 is excluded")
    if issues:
        level = "Low"
    else:
        level = "Good" if se < SE_GOOD else "Moderate"
    return level, "; ".join(issues) or "none"


def assign_confidence(r):
    """Explicit rules; see PHASE8_BIOLOGICAL_INTERPRETATION.md section 2.
    Expression-signal confidence = how strongly and consistently the data support a real, reproducible
    change in the stated direction. It does not express certainty about the exact fold change."""
    issues = []
    if float(r["fact_padj"]) >= 1e-3:
        issues.append("adjusted p-value >= 0.001")
    if abs(float(r["fact_log2FC_MLE"])) < 1:
        issues.append("|log2FC| < 1 (modest effect)")
    if r["fact_comparison"] == "C1":
        if r["fact_A3_excluded_significant"] != "TRUE" or r["fact_A3_direction_concordant"] != "TRUE":
            issues.append("not robust to excluding A3")
    if r["fact_gene_model_note"] != "none":
        issues.append("gene-level estimate may be ambiguous (overlapping/readthrough gene models)")
    de = "High" if not issues else ("Moderate" if len(issues) == 1 else "Low")

    fn = r["fact_function_evidence_class"]
    if "pseudogene" in r["fact_gene_biotype"]:
        bio = "Very low (pseudogene; no protein function)"
    elif not r["fact_uniprot_accession"] or fn == "none":
        bio = "Very low (no known function)"
    elif r["fact_uniprot_review_status"].startswith("UniProtKB unreviewed") or fn == "automatic only":
        bio = "Low (unreviewed / automatic annotation only)"
    elif fn.startswith("experimental"):
        bio = "Moderate" if de == "High" else "Low (expression signal has caveats)"
    else:
        bio = "Low (function not experimentally established)"
    return de, "; ".join(issues) or "none", bio


def main():
    for p in (OUT, OUT_VAL):
        if p.exists():
            raise SystemExit(f"Refusing to overwrite {p}")
    ann = {r["ensembl_gene_id"]: r for r in read(T / "phase8_target_gene_annotation.csv")}
    fun = {r["ensembl_gene_id"]: r for r in read(T / "phase8_functional_annotation.csv")}
    go = read(T / "phase8_functional_go_terms.csv")
    refs = read(T / "phase8_functional_references.csv")
    mapping = read(T / "phase8_uniprot_mapping.csv")
    attrs = read(T / "phase8_raw/ensembl97_attributes.csv")
    xrefs = read(T / "phase8_raw/ensembl97_xrefs.csv")
    de = {c: {r["gene_id"]: r for r in read(T / f)} for c, (_, f) in COMPARISONS.items()}
    de_noA3 = {r["gene_id"]: r for r in read(T / "C1_differential_expression_without_A3.csv")}
    sens = {r["gene_id"]: r for r in read(T / "C1_A3_sensitivity_comparison.csv")}
    groups = {c: {"control": [], "treatment": []} for c in COMPARISONS}
    for s in read(REPO_ROOT / "data/metadata/sample_metadata.csv"):
        groups[s["comparison"][:2]][s["treatment_role"]].append(s["sample_id"])

    order = list(ann)
    rows = []
    for g in order:
        a, f, d = ann[g], fun[g], de[ann[g]["comparison"]][g]
        comp = a["comparison"]
        lfc = float(d["log2FoldChange_MLE"])
        hv = sorted({x["value"] for x in attrs if x["gene_stable_id"] == g and x["object_type"] == "Gene"
                     and x["attrib_code"] == "havana_cv"})
        n_rt = len({x["object_stable_id"] for x in attrs if x["gene_stable_id"] == g and x["attrib_code"] == "readthrough_tra"})
        note = []
        if n_rt:
            note.append(f"{n_rt} of {a['n_transcripts']} transcripts flagged readthrough (Ensembl 97)")
        overlaps = [o for o in order if o != g and ann[o]["chromosome"] == a["chromosome"] and ann[o]["strand"] == a["strand"]
                    and int(ann[o]["start"]) <= int(a["end"]) and int(ann[o]["end"]) >= int(a["start"])]
        if overlaps:
            note.append("overlaps target gene model " + ", ".join(overlaps))
        r_hv = "; ".join(hv) or "none"
        gene_go = [x for x in go if x["ensembl_gene_id"] == g]
        auto = bool(f.get("uniprot_accession")) and all(x["go_evidence"].startswith("IEA") for x in gene_go) \
            and evidence_class(f.get("function_evidence_codes", "")) in ("automatic only", "none")
        # Reactome xrefs: 'description' holds the pathway name, 'primary_id' the stable ID
        reactome = sorted({x["description"] + " (" + x["primary_id"] + ")" for x in xrefs if x["gene_stable_id"] == g
                           and x["external_db"] == "Reactome_gene"})
        omim = sorted({x["display_label"].split(";")[0] + " (" + x["external_db"] + " " + x["primary_id"] + ")"
                       for x in xrefs if x["gene_stable_id"] == g and x["external_db"] in ("MIM_GENE", "MIM_MORBID")})
        pmids = sorted({x["pubmed_id"] for x in refs if x["ensembl_gene_id"] == g}, key=int)
        go_pmids = sorted({s.split(":")[1] for x in gene_go for s in x["evidence_sources"].split("; ") if s.startswith("PubMed:")}, key=int)
        r = {
            "record_order": order.index(g) + 1,
            "fact_comparison": comp, "fact_comparison_description": COMPARISONS[comp][0],
            "fact_target_symbol": a["target_symbol"], "fact_record_role": a["record_role"],
            "fact_ensembl_gene_id": g, "fact_versioned_gene_id": a["versioned_gene_id"],
            "fact_ensembl_description": a["description"], "fact_description_source": a["description_source"] or "none",
            "fact_hgnc_id": a["hgnc_id"] or "none", "fact_gene_biotype": a["gene_biotype"],
            "fact_location_grch38": f'chr{a["chromosome"]}:{a["start"]}-{a["end"]} ({a["strand"]})',
            "fact_n_transcripts": a["n_transcripts"], "fact_transcript_biotypes": a["transcript_biotypes"],
            "fact_gene_model_note": "; ".join(note) or "none",
            "fact_havana_curation_terms": r_hv,
            "fact_baseMean": f'{float(d["baseMean"]):.1f}',
            "fact_log2FC_MLE": d["log2FoldChange_MLE"], "fact_lfcSE_MLE": d["lfcSE_MLE"],
            "fact_approx_fold_change": f"{2 ** abs(lfc):.2f}x {'higher' if lfc > 0 else 'lower'} in treated",
            "fact_padj": d["padj"], "fact_log2FC_shrunk": d["log2FoldChange_shrunk"],
            "fact_shrinkage_collapsed": d["shrinkage_collapsed"], "fact_significant": d["significant"],
            "fact_direction": {"up": "higher in treated", "down": "lower in treated"}.get(d["direction"], d["direction"]),
            "fact_normalised_counts": "; ".join(f'{k[5:]}={v}' for k, v in d.items() if k.startswith("norm_")),
        }
        # Group means (normalised counts) and an approximate 95% Wald interval from the saved MLE and SE
        ctrl = [float(d["norm_" + s]) for s in groups[comp]["control"]]
        trt = [float(d["norm_" + s]) for s in groups[comp]["treatment"]]
        se = float(d["lfcSE_MLE"])
        lo, hi = lfc - 1.96 * se, lfc + 1.96 * se
        fc_int = (f"{2 ** lo:.1f}x to {2 ** hi:.1f}x higher" if lo > 0 else
                  f"{2 ** -hi:.1f}x to {2 ** -lo:.1f}x lower" if hi < 0 else
                  f"{2 ** -lo:.1f}x lower to {2 ** hi:.1f}x higher")
        r.update({
            "fact_control_samples": ", ".join(groups[comp]["control"]),
            "fact_control_mean_normalised": f"{sum(ctrl) / len(ctrl):.1f}",
            "fact_treated_mean_normalised": f"{sum(trt) / len(trt):.1f}",
            "fact_log2FC_95pct_wald_interval": f"{lo:.2f} to {hi:.2f}",
            "fact_fold_change_95pct_wald_interval": fc_int + " (approximate, log2FC +/- 1.96 x SE; Phase 6 values)",
        })
        if comp == "C1":
            n, s = de_noA3[g], sens[g]
            r.update({"fact_A3_excluded_log2FC_MLE": n["log2FoldChange_MLE"], "fact_A3_excluded_padj": n["padj"],
                      "fact_A3_excluded_significant": n["significant"],
                      "fact_A3_excluded_shrinkage_collapsed": n["shrinkage_collapsed"],
                      "fact_A3_sensitivity_category": s["significance_category"],
                      "fact_A3_direction_concordant": s["direction_concordant"]})
        else:
            for k in ("fact_A3_excluded_log2FC_MLE", "fact_A3_excluded_padj", "fact_A3_excluded_significant",
                      "fact_A3_excluded_shrinkage_collapsed", "fact_A3_sensitivity_category", "fact_A3_direction_concordant"):
                r[k] = "not applicable (C2)"
        r.update({
            "fact_uniprot_accession": f.get("uniprot_accession", ""),
            "fact_uniprot_review_status": f.get("review_status", ""),
            "fact_uniprot_other_used_accessions": f.get("other_used_accessions", ""),
            "fact_uniprot_protein_name": (f.get("protein_name", "") + (f' ({f["protein_name_type"]})' if f.get("protein_name_type") else "")),
            "fact_uniprot_gene_names": f.get("uniprot_gene_names", ""),
            "fact_uniprot_function": f.get("function", "") or "none",
            "fact_function_evidence_class": evidence_class(f.get("function_evidence_codes", "")),
            "fact_uniprot_subcellular_location": f.get("subcellular_location", "") or "none",
            "fact_location_evidence_class": evidence_class(f.get("subcellular_location_evidence_codes", "")),
            "fact_go_terms": "; ".join(f'{x["go_id"]} {x["go_term"]} [{x["aspect"]}; {x["go_evidence"]}]' for x in gene_go) or "none",
            "fact_electronic_or_automatic_only": "yes" if auto else ("not applicable" if not f.get("uniprot_accession") else "no"),
            "fact_reactome_ensembl97": "; ".join(reactome) or "none",
            "fact_omim_ensembl97": "; ".join(omim) or "none",
            "fact_supporting_pubmed_ids": "; ".join(f"PubMed:{p}" for p in sorted(set(pmids) | set(go_pmids), key=int)) or "none",
            "fact_annotation_records": (f"Ensembl 97 core db gene {a['versioned_gene_id']} (results/tables/phase8_raw/ensembl97_*.csv)"
                                        + (f"; UniProtKB 2026_04 {f['uniprot_accession']} (results/tables/phase8_raw/uniprot_2026_04/entries/{f['uniprot_accession']}.json)"
                                           if f.get("uniprot_accession") else "; no UniProtKB record mapped")),
        })
        de_conf, de_issues, bio_conf = assign_confidence(r)
        prec, prec_issues = effect_size_precision(r)
        r.update({"conf_expression_signal": de_conf, "conf_expression_signal_issues": de_issues,
                  "conf_effect_size_precision": prec, "conf_effect_size_precision_issues": prec_issues,
                  "conf_biological_interpretation": bio_conf})
        r.update(INTERPRETATION[g])
        rows.append(r)

    # MATR3 descriptive note: sum of normalised counts of the two overlapping models (C2 samples). Not tested.
    m1 = de["C2"]["ENSG00000280987"]; m2 = de["C2"]["ENSG00000015479"]
    sums = {k[5:]: round(float(m1[k]) + float(m2[k])) for k in m1 if k.startswith("norm_")}
    matr3_sum = "; ".join(f"{k}={v}" for k, v in sums.items())
    for r in rows:
        r["fact_matr3_combined_normalised_counts"] = matr3_sum + " (descriptive sum of both MATR3 models; not tested)" \
            if r["fact_target_symbol"] == "MATR3" else "not applicable"

    # ---------------- validation ----------------
    checks = []

    def check(rec, name, ok, detail=""):
        checks.append({"record": rec, "check": name, "result": "PASS" if ok else "FAIL", "detail": detail})

    check("all", "Nine records, unique gene IDs, both MATR3 IDs present",
          len(rows) == 9 and len({r["fact_ensembl_gene_id"] for r in rows}) == 9
          and {"ENSG00000280987", "ENSG00000015479"} <= {r["fact_ensembl_gene_id"] for r in rows})
    known_pmids = {x["pubmed_id"] for x in refs} | {s.split(":")[1] for x in go for s in x["evidence_sources"].split("; ") if s.startswith("PubMed:")}
    for r in rows:
        g, sym, comp = r["fact_ensembl_gene_id"], r["fact_target_symbol"], r["fact_comparison"]
        rec = f"{sym} {g}"
        d = de[comp][g]
        check(rec, "Phase 6 values copied exactly from the record's own comparison",
              r["fact_log2FC_MLE"] == d["log2FoldChange_MLE"] and r["fact_padj"] == d["padj"]
              and r["fact_log2FC_shrunk"] == d["log2FoldChange_shrunk"], f"{comp}: LFC {float(d['log2FoldChange_MLE']):.3f}, padj {fmt_p(d['padj'])}")
        check(rec, "Direction matches the sign of the MLE fold change and Phase 6 direction",
              (float(r["fact_log2FC_MLE"]) > 0) == (d["direction"] == "up") and d["significant"] == "TRUE")
        check(rec, "Gene ID/version matches the Ensembl 97 annotation table", r["fact_versioned_gene_id"] == ann[g]["versioned_gene_id"])
        check(rec, "UniProt accession and status match the Step 3 functional table",
              r["fact_uniprot_accession"] == fun[g].get("uniprot_accession", "") and r["fact_uniprot_review_status"] == fun[g]["review_status"])
        interp_text = " ".join(v for k, v in r.items() if k.startswith("interp_"))
        cited = set(re.findall(r"PubMed:(\d+)", interp_text)) | {p.split(":")[1] for p in r["fact_supporting_pubmed_ids"].split("; ") if p.startswith("PubMed:")}
        check(rec, "Every cited PubMed ID is in the saved UniProt records", cited <= known_pmids, ", ".join(sorted(cited)) or "none cited")
        if comp == "C1":
            n = de_noA3[g]
            check(rec, "A3 sensitivity values present and copied exactly",
                  r["fact_A3_excluded_log2FC_MLE"] == n["log2FoldChange_MLE"] and r["fact_A3_excluded_padj"] == n["padj"]
                  and r["fact_A3_sensitivity_category"] == sens[g]["significance_category"], r["fact_A3_sensitivity_category"])
            check(rec, "Limitations state that C1 without A3 is exploratory (two control replicates)",
                  "two control replicates" in r["interp_limitations"] and "exploratory" in r["interp_limitations"])
        else:
            check(rec, "C2 record carries no C1 or A3 statistics", r["fact_A3_excluded_log2FC_MLE"] == "not applicable (C2)")
        if "pseudogene" in r["fact_gene_biotype"]:
            check(rec, "Pseudogene: no UniProt record, no protein function, interpretation does not claim a protein",
                  not r["fact_uniprot_accession"] and r["fact_uniprot_function"] == "none"
                  and "No protein function" in r["interp_known_function"]
                  and "That a protein is made" in r["interp_cannot_conclude"])
        if g == "ENSG00000286192":
            check(rec, "AC069288.1: no MRM2 annotation used in the interpretation",
                  "methyltransferase" not in interp_text.lower() and "rRNA" not in interp_text
                  and r["fact_uniprot_function"] == "none", "unreviewed record flagged; its GO term appears only in fact columns")
            check(rec, "AC069288.1: automatic/electronic-only annotation flagged", r["fact_electronic_or_automatic_only"] == "yes")
        if g == "ENSG00000280987":
            check(rec, "MATR3 primary: reviewed P43243 not used; unreviewed record flagged automatic-only",
                  r["fact_uniprot_accession"] != "P43243" and r["fact_electronic_or_automatic_only"] == "yes"
                  and any(m["ensembl_gene_id"] == g and m["uniprot_accession"] == "P43243" and m["decision"].startswith("ambiguous") for m in mapping))
        if sym == "MATR3":
            other = "ENSG00000015479" if g == "ENSG00000280987" else "ENSG00000280987"
            check(rec, "MATR3: kept as a separate record; text refers to the other model without claiming independence",
                  other in interp_text and "independent" in r["interp_cannot_conclude"])
        if g == "ENSG00000213064":
            check(rec, "SFT2D2: function evidence is by similarity only and is stated as such",
                  r["fact_function_evidence_class"] == "by similarity only" and "similarity" in r["interp_known_function"]
                  and "not by experiments" in r["interp_known_function"])
        check(rec, "Interpretation does not claim proof of function or mechanism",
              not re.search(r"\b(proves?|demonstrates?|confirms?|establishes)\b", interp_text, re.I))
        check(rec, "Biological-interpretation confidence is never 'High'", not r["conf_biological_interpretation"].startswith("High"),
              r["conf_biological_interpretation"])
        check(rec, "All interpretation fields are filled",
              all(r.get(k) for k in ("interp_known_function", "interp_process_pathway", "interp_expression_change",
                                     "interp_hypotheses", "interp_cannot_conclude", "interp_limitations")))
        # --- Revision 1: signal confidence vs effect-size precision ---
        check(rec, "Signal confidence and effect-size precision recomputed from the stated rules",
              (r["conf_expression_signal"], r["conf_expression_signal_issues"])[0] == assign_confidence(r)[0]
              and r["conf_effect_size_precision"] == effect_size_precision(r)[0],
              f'signal {r["conf_expression_signal"]}; precision {r["conf_effect_size_precision"]}')
        check(rec, "Shrinkage collapse affects effect-size precision only, not signal confidence",
              "collapse" not in r["conf_expression_signal_issues"]
              and (("collapse" in r["conf_effect_size_precision_issues"]) ==
                   (r["fact_shrinkage_collapsed"] == "TRUE" or r.get("fact_A3_excluded_shrinkage_collapsed") == "TRUE")))
        lfc, se = float(r["fact_log2FC_MLE"]), float(r["fact_lfcSE_MLE"])
        check(rec, "Wald interval recomputes from the Phase 6 log2FC and SE",
              r["fact_log2FC_95pct_wald_interval"] == f"{lfc - 1.96 * se:.2f} to {lfc + 1.96 * se:.2f}",
              r["fact_log2FC_95pct_wald_interval"])
        # --- Revision 2: overlapping-gene wording ---
        check(rec, "No claim that reads are divided or redistributed between gene models",
              not OVERLAP_BANNED.search(interp_text))
        if r["fact_gene_model_note"] != "none":
            check(rec, "Overlapping/readthrough record: workflow note and 'may not be independent' wording present",
                  WORKFLOW_NOTE in r["interp_limitations"] and "may not be independent" in interp_text)
    tub = next(r for r in rows if r["fact_ensembl_gene_id"] == "ENSG00000127589")
    check("TUBBP1 ENSG00000127589", "TUBBP1: High signal confidence kept, effect-size precision Low, text separates direction from size",
          tub["conf_expression_signal"] == "High" and tub["conf_effect_size_precision"] == "Low"
          and "direction of the change is well supported" in tub["interp_expression_change"]
          and "exact size of the change is not well determined" in tub["interp_expression_change"],
          tub["conf_effect_size_precision_issues"])

    # --- Report consistency (results/PHASE8_BIOLOGICAL_INTERPRETATION.md) ---
    md_path = REPO_ROOT / "results/PHASE8_BIOLOGICAL_INTERPRETATION.md"
    md = md_path.read_text(encoding="utf-8")
    check("report", "Report has no claim that reads are divided or redistributed between gene models",
          not OVERLAP_BANNED.search(md), "; ".join(m.group(0) for m in OVERLAP_BANNED.finditer(md)) or "none found")
    check("report", "Report describes the actual counting workflow (StringTie, prepDE.py, Phase 5 sum)",
          all(s in md for s in ("StringTie", "prepDE.py", "Phase 5", "may not be independent")))
    check("report", "Report defines expression-signal confidence separately from effect-size precision",
          "Expression-signal confidence" in md and "Effect-size precision" in md
          and "not certainty about the exact" in md)
    summary = md[md.index("## 3. Summary"):md.index("### A3 sensitivity")]
    for r in rows:
        line = next((l for l in summary.splitlines() if r["fact_versioned_gene_id"] in l), "")
        lfc_s = f'{float(r["fact_log2FC_MLE"]):.2f}'.replace("-", "−")
        p = f'{float(r["fact_padj"]):.2e}'
        check("report", f'Summary row for {r["fact_target_symbol"]} {r["fact_ensembl_gene_id"]} matches the table',
              lfc_s in line and p in line and f'| {r["conf_expression_signal"]} |' in line
              and f'| {r["conf_effect_size_precision"]} |' in line, line[:60] + "...")
    tub_sec = md[md.index("### 4.2 TUBBP1"):md.index("### 4.3 INHBE")]
    tub_fc = tub["fact_fold_change_95pct_wald_interval"].split(" (")[0]
    lo_fc, hi_fc = re.findall(r"([\d.]+)x", tub_fc)
    check("report", "TUBBP1 section states High signal, Low precision and the table's interval",
          "Expression signal: High" in tub_sec and "Effect-size precision: Low" in tub_sec
          and f"{float(lo_fc):.0f}" in tub_sec and f"{float(hi_fc):.0f}" in tub_sec, tub_fc)

    fields = list(rows[0].keys())
    fields = [k for k in fields if not k.startswith(("conf_", "interp_"))] + \
             [k for k in fields if k.startswith("conf_")] + [k for k in fields if k.startswith("interp_")]
    with open(OUT, "x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    with open(OUT_VAL, "x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["record", "check", "result", "detail"])
        w.writeheader()
        w.writerows(checks)
    n_fail = sum(c["result"] == "FAIL" for c in checks)
    for c in checks:
        if c["result"] == "FAIL":
            print("FAIL:", c["record"], "-", c["check"], c["detail"])
    print(f"wrote {OUT.relative_to(REPO_ROOT)} ({len(rows)} rows, {len(fields)} columns)")
    print(f"wrote {OUT_VAL.relative_to(REPO_ROOT)}: {len(checks) - n_fail}/{len(checks)} PASS")
    for r in rows:
        print(f'  {r["fact_comparison"]} {r["fact_target_symbol"]:<11} {r["fact_ensembl_gene_id"]}  LFC {float(r["fact_log2FC_MLE"]):6.2f}  '
              f'padj {fmt_p(r["fact_padj"]):>9}  signal: {r["conf_expression_signal"]:<8} precision: '
              f'{r["conf_effect_size_precision"]:<8} bio: {r["conf_biological_interpretation"]}\n'
              f'      signal issues: {r["conf_expression_signal_issues"]} | precision issues: '
              f'{r["conf_effect_size_precision_issues"]} | 95% FC: {r["fact_fold_change_95pct_wald_interval"].split(" (")[0]}')


if __name__ == "__main__":
    main()
