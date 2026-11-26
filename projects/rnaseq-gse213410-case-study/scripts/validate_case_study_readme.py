"""
Validate the case-study README against the machine-readable outputs.

Checks (read only; writes nothing):
  - every number quoted in the README is rebuilt from the saved tables and
    must appear in the README text exactly as formatted;
  - the Phase 8 summary table matches phase8_gene_interpretation.csv;
  - the GSEA theme tables cover every significant Phase 7 gene set exactly
    once, agree with each other, and keep the reviewed classifications;
  - every relative link and image resolves to a file in the repository, and
    every in-page anchor matches a heading;
  - original case-study figures are byte-identical to the committed versions;
  - the README and case-study files contain no local computer paths.

Run from anywhere:
    python projects/rnaseq-gse213410-case-study/scripts/validate_case_study_readme.py
Exit code 0 only if every check passes.
"""

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
CASE = REPO / "projects/rnaseq-gse213410-case-study"
README = CASE / "README.md"
T = REPO / "results/tables"
text = README.read_text(encoding="utf-8")
flat = re.sub(r"\s+", " ", text)
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail else ""))


def quoted(name, *fragments):
    missing = [f for f in fragments if f not in flat]
    check(name, not missing, "missing: " + "; ".join(missing) if missing else f"{len(fragments)} values")


def n(x):
    return f"{int(x):,}"


# ---------------------------------------------------------------- Phase 6 table
de = pd.read_csv(T / "phase6_de_summary.csv").set_index("analysis")
rows = {"C1 (all samples, primary)": "C1 - 5-OXO-ETE vs Ethanol, A3 included",
        "C1 without A3 (sensitivity)": "C1 - 5-OXO-ETE vs Ethanol, excluding A3",
        "C2 (all samples)": "C2 - Calcitriol vs DMSO"}
for key, label in rows.items():
    r = de.loc[key]
    quoted(f"DE table row: {label}",
           f"| {label} | {n(r.genes_after_prefilter)} | {n(r.genes_with_pvalue)} | "
           f"{n(r.genes_with_padj)} | {n(r['significant_padj_0.05'])} |")
    check(f"{label}: valid p-values = retained - Cook's outliers",
          r.genes_after_prefilter - r.cooks_outliers == r.genes_with_pvalue)
    check(f"{label}: adjusted p-values = valid p-values - independent filter",
          r.genes_with_pvalue - r.independent_filter_removed == r.genes_with_padj)
quoted("Cook's outlier and independent-filter counts",
       f"There were {int(de.iloc[0].cooks_outliers)}, {int(de.iloc[1].cooks_outliers)} and "
       f"{int(de.iloc[2].cooks_outliers)} such genes",
       f"removed {int(de.iloc[1].independent_filter_removed)} low-mean genes")
c1, c2 = de.iloc[0], de.iloc[2]
quoted("Volcano up/down counts",
       f"found {int(c1['significant_padj_0.05'])} genes with adjusted p < 0.05: {int(c1.up)} higher and {int(c1.down)} lower",
       f"found {int(c2['significant_padj_0.05'])} genes with adjusted p < 0.05: {int(c2.up)} higher and {int(c2.down)} lower")

# ---------------------------------------------------------------- A3 sensitivity
s = pd.read_csv(T / "C1_A3_sensitivity_summary.csv").set_index("metric")["value"]
comp = pd.read_csv(T / "C1_A3_sensitivity_comparison.csv")
both = comp[comp.tested_in == "both"]
plotted_wo = int((both.significance_category == "significant_without_A3_only").sum())
not_plotted_wo = int(s["significant without A3 only"]) - plotted_wo
quoted("Sensitivity denominators",
       f"{n(s['genes tested in both'])} genes were retained after prefiltering in **both** analyses",
       f"{n(both[['padj_with', 'padj_without']].notna().all(axis=1).sum())} have an adjusted p-value in both",
       f"The {n((comp.tested_in == 'without_A3_only').sum())} genes retained only")
quoted("Sensitivity significance counts",
       f"{int(s['significant with A3'])} with A3 included", f"{int(s['significant without A3'])} with A3 excluded",
       f"{int(s['significant in both'])} in both", f"{int(s['significant with A3 only'])} with A3 included only",
       f"{int(s['significant without A3 only'])} with A3 excluded only, of which {plotted_wo} are among",
       f"and {not_plotted_wo} were retained only without A3")
quoted("Sensitivity Pearson r", f"**{s['Pearson r of MLE LFC (genes tested in both)']:.3f}**",
       f"Pearson correlation was {s['Pearson r of MLE LFC (genes tested in both)']:.3f}")
check("Direction concordance among genes significant in both is 1",
      s["direction concordance among genes significant in both"] == 1)
p7 = (REPO / "results/PHASE7_REPORT.md").read_text(encoding="utf-8")
m = re.search(r"NES correlation across the ([\d,]+) shared sets: Pearson r = (\d\.\d+)", p7)
check("Phase 7 GSEA sensitivity correlation quoted as in PHASE7_REPORT.md",
      bool(m) and f"Pearson r = {m.group(2)} across {m.group(1)} shared gene sets" in flat,
      m.group(0) if m else "pattern not found in report")

# ---------------------------------------------------------------- QC and PCA
prov = pd.read_csv(T / "count_provenance_check.csv").set_index("sample_id")
a3 = prov.loc["A3"]
med = prov.transcript_level_counts.median()
quoted("A3 QC figures",
       f"approximately {a3.spots / 1e6:.1f} million read pairs were sequenced for A3",
       f"approximately {a3.transcript_level_counts / 1e6:.1f} million counts were assigned",
       f"{a3.transcript_level_counts / med:.0%} of the {med / 1e6:.1f} million median",
       f"range across all samples {prov.spots.min() / 1e6:.1f}–{prov.spots.max() / 1e6:.1f} million")
others = prov.drop("A3").transcript_level_counts
check("Other samples' assigned counts are roughly 38 to 55 million",
      37.5e6 <= others.min() < 38.5e6 and 54.5e6 <= others.max() < 55.5e6,
      f"{others.min() / 1e6:.1f}–{others.max() / 1e6:.1f} M")
pca = pd.read_csv(T / "gene_level_A3_sensitivity.csv")
quoted("PCA variance percentages",
       f"which explains {pca.iloc[0].PC1_pct_variance:.1f}% of the variance",
       f"which explain {pca.iloc[2].PC1_pct_variance:.1f}% and {pca.iloc[2].PC2_pct_variance:.1f}% of the variance")

# ---------------------------------------------------------------- count provenance
tm = pd.read_csv(T / "transcript_mapping_summary.csv").set_index("metric")["value"]
quoted("Transcript mapping figures",
       f"The {n(tm['ENST_transcripts'])} Ensembl transcripts", f"The {n(tm['TCONS_transcripts'])} novel StringTie transcripts",
       f"holding {float(tm['pct_counts_TCONS']):.1f}% of counts", f"{n(tm['total_transcripts'])} transcript rows")
check("Every ENST maps to exactly one gene (plain sum)",
      tm["ENST_mapped_to_one_gene"] == tm["ENST_transcripts"] and int(tm["ENST_ambiguous_multiple_genes"]) == 0)
p6 = (REPO / "results/PHASE6_REPORT.md").read_text(encoding="utf-8")
check("StringTie 1.3.3b and prepDE.py documented in PHASE6_REPORT.md",
      "StringTie 1.3.3b" in p6 and "prepDE.py" in p6)
check("README does not claim a raw-read pipeline was run here",
      "did **not** run a raw-read-to-count pipeline" in text)

# ---------------------------------------------------------------- Phase 8 table
interp = pd.read_csv(T / "phase8_gene_interpretation.csv")
for _, r in interp.iterrows():
    gid = r.fact_ensembl_gene_id
    line = next((ln for ln in text.splitlines() if ln.startswith(f"| {r.fact_comparison} |") and f"({gid})" in ln), None)
    if line is None:
        check(f"Phase 8 row {r.fact_target_symbol} {gid}", False, "row missing")
        continue
    cells = [c.strip() for c in line.strip("|").split("|")]
    lfc = f"{r.fact_log2FC_MLE:.2f}".replace("-", "−")
    want = [lfc, f"{r.fact_padj:.2e}", r.conf_expression_signal.split(" (")[0],
            r.conf_effect_size_precision.split(" (")[0], r.conf_biological_interpretation.split(" (")[0]]
    check(f"Phase 8 row {r.fact_target_symbol} {gid}", cells[3:8] == want, f"README {cells[3:8]} vs table {want}")
check("Phase 8 table has nine records", len(interp) == 9 and interp.fact_ensembl_gene_id.nunique() == 9)
for gid in ("ENSG00000280987", "ENSG00000015479"):
    r = interp[interp.fact_ensembl_gene_id == gid].iloc[0]
    quoted(f"MATR3 section row {gid}",
           f"| {r.fact_log2FC_MLE:.2f}".replace("-", "−") + f" | {r.fact_padj:.2e} |")
check("Primary MATR3 record has no P43243 assignment",
      interp.set_index("fact_ensembl_gene_id").at["ENSG00000280987", "fact_uniprot_accession"] != "P43243")

# ---------------------------------------------------------------- GSEA themes
ts = pd.read_csv(CASE / "tables/gsea_theme_summary.csv")
for comp_id in ("C1", "C2"):
    g = pd.read_csv(T / f"phase7_{comp_id}_GSEA_results.csv")
    sig = g[g.significant == True]          # noqa: E712
    t = ts[ts.comparison == comp_id].set_index("theme").drop(columns="comparison")
    check(f"{comp_id} theme totals equal Phase 7 significant sets",
          int(t.values.sum()) == len(sig) and int(t.sets_higher_in_treatment.sum())
          == int((sig.direction == "higher in treatment").sum()))
    up, dn = int(t.sets_higher_in_treatment.sum()), int(t.sets_lower_in_treatment.sum())
    quoted(f"{comp_id} GSEA totals", f"{len(sig)} significant gene sets", f"{up} higher and {dn} lower")
t1 = ts[ts.comparison == "C1"].set_index("theme")
t2 = ts[ts.comparison == "C2"].set_index("theme")
H = "Histone-gene dominated (shared histone genes)"
quoted("C1 theme counts", f"{t1.at[H, 'sets_lower_in_treatment']} of the "
       f"{int(t1.sets_lower_in_treatment.sum())} \"lower in treatment\" sets are histone-gene dominated",
       f"RNA processing, splicing and transport ({t1.at['RNA processing, splicing and transport', 'sets_higher_in_treatment']} sets)",
       f"infection, immune and disease-labelled sets ({t1.at['Infection, immune and disease-labelled sets', 'sets_higher_in_treatment']})",
       f"A further {t1.at['Other', 'sets_higher_in_treatment']} \"higher in treatment\" sets matched no theme rule")
INF_T, RNA_T = "Infection, immune and disease-labelled sets", "RNA processing, splicing and transport"
named = t1.drop(index="Other").sets_higher_in_treatment.sort_values(ascending=False)
check("C1 infection/immune 'higher' sets = 23 (after the reviewed viral-pathway reclassification)",
      t1.at[INF_T, "sets_higher_in_treatment"] == 23, f"{t1.at[INF_T, 'sets_higher_in_treatment']}")
check("C1 RNA processing 'higher' sets = 18", t1.at[RNA_T, "sets_higher_in_treatment"] == 18,
      f"{t1.at[RNA_T, 'sets_higher_in_treatment']}")
check("C1: infection/immune and RNA processing are the two largest named 'higher' themes, in that order",
      list(named.index[:2]) == [INF_T, RNA_T] and named.iloc[1] > named.iloc[2],
      ", ".join(f"{k}: {v}" for k, v in named.head(3).items()))
pos_inf = flat.find(f"infection, immune and disease-labelled sets ({t1.at[INF_T, 'sets_higher_in_treatment']})")
pos_rna = flat.find(f"RNA processing, splicing and transport ({t1.at[RNA_T, 'sets_higher_in_treatment']} sets)")
check("README reports the two largest C1 'higher' themes in rank order", -1 < pos_inf < pos_rna,
      f"positions {pos_inf}, {pos_rna}")
lookup_c1 = {r.Description: r.theme for r in pd.read_csv(CASE / "tables/gsea_theme_assignments.csv").itertuples()
             if r.comparison == "C1"}
viral = ["Nuclear import of Rev protein", "Interactions of Rev with host cellular proteins",
         "Transport of Ribonucleoproteins into the Host Nucleus"]
check("Reviewed viral-context pathways are grouped with infection/immune",
      all(lookup_c1.get(v_) == INF_T for v_ in viral), "; ".join(f"{v_}: {lookup_c1.get(v_)}" for v_ in viral))
quoted("C2 theme counts",
       f"histone-gene dominated sets ({t2.at[H, 'sets_higher_in_treatment']})",
       f"cell cycle, mitosis and chromosome organisation ({t2.at['Cell cycle, mitosis and chromosome organisation', 'sets_higher_in_treatment']})",
       f"DNA replication and repair ({t2.at['DNA replication and repair', 'sets_higher_in_treatment']})",
       f"membrane channels, transporters and receptors ({t2.at['Membrane channels, transporters and receptors', 'sets_lower_in_treatment']})",
       f"extracellular matrix, collagen and adhesion ({t2.at['Extracellular matrix, collagen and adhesion', 'sets_lower_in_treatment']})")

# Assignment table: every significant Phase 7 set exactly once, and consistent with the summary table
asg = pd.read_csv(CASE / "tables/gsea_theme_assignments.csv")
for comp_id in ("C1", "C2"):
    g = pd.read_csv(T / f"phase7_{comp_id}_GSEA_results.csv")
    sig = g[g.significant == True]          # noqa: E712
    a = asg[asg.comparison == comp_id]
    keys_sig = set(zip(sig.database, sig.ID, sig.direction))
    keys_a = set(zip(a.database, a.ID, a.direction))
    check(f"{comp_id} assignments: every significant Phase 7 set exactly once",
          len(a) == len(sig) == int((g["p.adjust"] < 0.05).sum()) and keys_a == keys_sig
          and not a.duplicated(["database", "ID"]).any() and a.theme.notna().all(),
          f"{len(a)} rows, {len(keys_sig - keys_a)} missing, {len(keys_a - keys_sig)} extra")
    check(f"{comp_id} assignments match Phase 7 counts by database and direction",
          sig.groupby(["database", "direction"]).size().equals(a.groupby(["database", "direction"]).size()))
    xt = a.groupby(["theme", "direction"]).size().unstack(fill_value=0)
    s = ts[ts.comparison == comp_id].set_index("theme")
    same = all(int(xt.get("higher in treatment", pd.Series(dtype=int)).get(th, 0)) == r.sets_higher_in_treatment
               and int(xt.get("lower in treatment", pd.Series(dtype=int)).get(th, 0)) == r.sets_lower_in_treatment
               for th, r in s.iterrows()) and set(xt.index) <= set(s.index)
    check(f"{comp_id} summary table equals the assignment table", same)

# Corrected classifications (regression checks for the 12 rule errors and the reviewed ambiguous sets)
expected = {
    ("C1", "external encapsulating structure"): "Extracellular matrix, collagen and adhesion",
    ("C2", "external encapsulating structure"): "Extracellular matrix, collagen and adhesion",
    ("C2", "external side of plasma membrane"): "Membrane channels, transporters and receptors",
    ("C2", "detection of external stimulus"): "Membrane channels, transporters and receptors",
    ("C1", "RORA,B,C and NR1D1 (REV-ERBA) regulate gene expression"): "Other",
    ("C2", "DNA-binding transcription activator activity, RNA polymerase II-specific"): "Other",
    ("C1", "Regulation of CDH1 mRNA translation by microRNAs"): "RNA processing, splicing and transport",
}
for c in ("C1", "C2"):
    for n_ in ("signaling receptor regulator activity", "receptor ligand activity",
               "signaling receptor activator activity", "helicase activity",
               "nucleic acid conformation isomerase activity"):
        expected[(c, n_)] = "Other"
for n_ in ("septin ring", "septin complex"):
    expected[("C1", n_)] = "Other"
expected[("C2", "Post-chaperonin tubulin folding pathway")] = "Other"
lookup = {(r.comparison, r.Description): r.theme for r in asg.itertuples()}
wrong = [f"{k[0]} {k[1]}: {lookup.get(k)}" for k, v in expected.items() if lookup.get(k) != v]
check("Reviewed classifications have their corrected themes", not wrong,
      "; ".join(wrong) if wrong else f"{len(expected)} rows")
check("No RNA-processing set matched only through 'RNA' inside another word",
      not asg[(asg.theme == "RNA processing, splicing and transport")
              & asg.Description.str.contains("external|REV-ERB|RNA polymerase", case=False)].shape[0])

# Histone threshold and its sensitivity, as quoted
hs = pd.read_csv(CASE / "tables/gsea_histone_threshold_sensitivity.csv")
v = {(r.comparison, round(r.threshold, 2)): r.histone_dominated_sets for r in hs.itertuples()}
check("Histone-dominated counts at 50% equal the summary table",
      v[("C1", 0.5)] == int(t1.loc[H, ["sets_higher_in_treatment", "sets_lower_in_treatment"]].sum())
      and v[("C2", 0.5)] == int(t2.loc[H, ["sets_higher_in_treatment", "sets_lower_in_treatment"]].sum()))
quoted("Histone threshold and sensitivity",
       "at least 50% of its leading-edge genes are histone genes",
       f"{v[('C1', 0.4)]}, {v[('C1', 0.5)]} and {v[('C1', 0.6)]} C1 sets at 40%, 50% and 60%",
       f"{v[('C2', 0.4)]}, {v[('C2', 0.5)]} and {v[('C2', 0.6)]} C2 sets",
       "approximate groupings for presentation, not definitive biological classifications")
nuc = asg[(asg.comparison == "C2") & (asg.Description == "nucleosome organization")].iloc[0]
check("Near-threshold example is below 50% and keyword-assigned", nuc.histone_fraction_leading_edge < 0.5
      and nuc.theme == "Chromatin and epigenetic regulation", f"{nuc.histone_fraction_leading_edge:.4f}")
quoted("Near-threshold example", f"has {nuc.histone_fraction_leading_edge:.0%} and is counted under chromatin")

# Wording corrections from the review
quoted("INHBE wording limited to C1", "the only C1 target with experimentally supported function")
check("No unqualified 'only target with experimentally supported function'",
      "the only target with experimentally supported function" not in flat)
quoted("TCONS limitation wording",
       f"The {n(tm['TCONS_transcripts'])} novel StringTie transcripts (TCONS), holding "
       f"{float(tm['pct_counts_TCONS']):.1f}% of counts, have no Ensembl annotation")
quoted("Phase 8 selection described honestly", "chose these genes judgementally",
       "No formal selection criteria were recorded")

# ---------------------------------------------------------------- status statement
check("Phase-completion statement updated", "Eight analytical phases completed" in text
      and "Seven analytical phases" not in text and "completed in seven phases" not in text)

# ---------------------------------------------------------------- links and anchors
def slug(heading):
    s_ = heading.strip().lower()
    s_ = re.sub(r"[^\w\- ]", "", s_)
    return s_.replace(" ", "-")


anchors = {slug(h) for h in re.findall(r"^#{1,6} (.+)$", text, flags=re.M)}
targets = re.findall(r"!?\[[^\]]*\]\(([^)\s]+)\)", text)
bad_files, bad_anchors, n_files = [], [], 0
for tgt in targets:
    if tgt.startswith(("http://", "https://")):
        continue
    path, _, frag = tgt.partition("#")
    if path:
        n_files += 1
        resolved = (CASE / path).resolve()
        if not resolved.exists() or REPO.resolve() not in resolved.parents:
            bad_files.append(tgt)
    elif frag and frag not in anchors:
        bad_anchors.append(tgt)
check("All relative links and images resolve inside the repository", not bad_files,
      ", ".join(bad_files) if bad_files else f"{n_files} links")
check("All in-page anchors match a heading", not bad_anchors, ", ".join(bad_anchors))
images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
check("Every image has descriptive alt text", all(len(a) > 20 for a in re.findall(r"!\[([^\]]*)\]", text)),
      f"{len(images)} images")

# ---------------------------------------------------------------- originals unchanged
originals = ["01_library_sizes.png", "11_gene_level_pca.png", "13_C1_volcano.png", "14_C2_volcano.png",
             "25_C1_A3_sensitivity_LFC_comparison.png", "27_C1_GSEA.png", "30_C2_GSEA.png"]
changed = []
for f in originals:
    rel = f"projects/rnaseq-gse213410-case-study/figures/{f}"
    head = subprocess.run(["git", "rev-parse", f"HEAD:{rel}"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    work = subprocess.run(["git", "hash-object", rel], cwd=REPO, capture_output=True, text=True).stdout.strip()
    if not head or head != work:
        changed.append(f)
check("Original case-study figures identical to committed versions", not changed,
      ", ".join(changed) if changed else f"{len(originals)} figures")

# ---------------------------------------------------------------- privacy
pattern = re.compile(r"[A-Za-z]:\\\\?Users|/c/Users|AppData|@gmail\.com", re.IGNORECASE)
scan = [README, CASE / "tables/gsea_theme_summary.csv", CASE / "tables/gsea_theme_assignments.csv",
        CASE / "tables/gsea_histone_threshold_sensitivity.csv"]
hits = [p.name for p in scan if pattern.search(p.read_text(encoding="utf-8"))]
check("No local paths or e-mail addresses in README or case-study tables", not hits, ", ".join(hits))

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
