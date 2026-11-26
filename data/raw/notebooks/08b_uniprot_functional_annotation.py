"""
Phase 8, Step 3: UniProtKB functional annotation for the target genes.

Sources
  - UniProtKB release 2026_04 (human records only), via https://rest.uniprot.org.
    Every response is checked for the header X-UniProt-Release: 2026_04; the
    script stops if any response comes from a different release. The REST API
    serves only the current release, so a rerun after 2026_04 is replaced will
    stop; the saved raw responses are then the record of this retrieval.
  - Ensembl release 97 (already used in Step 2):
      results/tables/phase8_raw/ensembl97_xrefs.csv and ensembl97_transcripts.csv
      Homo_sapiens.GRCh38.pep.all.fa.gz from
      https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep/
      (protein sequences, used to verify mappings by exact sequence identity).

Mapping rules (no mapping by gene symbol)
  A UniProt accession is a candidate for an Ensembl 97 gene only through an
  explicit cross-reference:
    (a) Ensembl 97 translation-level xref (Uniprot/SWISSPROT or Uniprot/SPTREMBL), or
    (b) the UniProt 2026_04 entry's own Ensembl cross-reference naming the gene ID.
  Ensembl's gene-level "Uniprot_gn" links are excluded (assigned via gene names).
  Reviewed (Swiss-Prot) candidates are preferred. Unreviewed (TrEMBL) entries are
  fetched and used only for genes with no reviewed candidate, and are labelled.
  If a fetched entry cross-references a transcript that belonged to a different
  target gene in Ensembl 97, that link is recorded as ambiguous and NOT transferred.
  Every gene-accession pair is checked for exact protein-sequence identity
  between Ensembl 97 translations and the UniProt canonical/isoform sequences.

Outputs (the script refuses to overwrite anything)
  Raw:     results/tables/phase8_raw/uniprot_2026_04/   (responses exactly as received + logs)
  Curated: results/tables/phase8_uniprot_mapping.csv
           results/tables/phase8_functional_annotation.csv
           results/tables/phase8_functional_go_terms.csv
           results/tables/phase8_functional_references.csv
           results/tables/phase8_functional_validation.csv
  The curated tables contain UniProt/Ensembl facts only; no biological interpretation.

Run from anywhere:
    python data/raw/notebooks/08b_uniprot_functional_annotation.py [cache_dir]
"""

import csv
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
TABLES = REPO_ROOT / "results/tables"
STEP2_RAW = TABLES / "phase8_raw"
RAW = STEP2_RAW / "uniprot_2026_04"
TARGET_TABLE = TABLES / "phase8_target_gene_annotation.csv"

EXPECTED_RELEASE = "2026_04"
HUMAN = 9606
API = "https://rest.uniprot.org"
PEP_BASE = "https://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/pep"
PEP_FILE = "Homo_sapiens.GRCh38.pep.all.fa.gz"
USER_AGENT = "bioinformatics-portfolio-annotation/1.0"

CURATED = {
    "mapping": TABLES / "phase8_uniprot_mapping.csv",
    "function": TABLES / "phase8_functional_annotation.csv",
    "go": TABLES / "phase8_functional_go_terms.csv",
    "refs": TABLES / "phase8_functional_references.csv",
    "validation": TABLES / "phase8_functional_validation.csv",
}

request_log = []


def now_utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def http_get(url, retries=4):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.status, dict(r.headers), r.read()
        except Exception as e:  # network hiccup: retry with backoff
            last = e
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"Request failed after {retries} attempts: {url} ({last})")


def uniprot_get(path, params, save_as):
    url = f"{API}/{path}" + ("?" + urllib.parse.urlencode(params) if params else "")
    status, headers, body = http_get(url)
    release = headers.get("X-UniProt-Release", "")
    dest = RAW / save_as
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "xb") as f:
        f.write(body)
    request_log.append({"url": url, "http_status": status, "x_uniprot_release": release,
                        "x_uniprot_release_date": headers.get("X-UniProt-Release-Date", ""),
                        "retrieved_utc": now_utc(), "bytes": len(body),
                        "sha256": hashlib.sha256(body).hexdigest(), "saved_as": f"uniprot_2026_04/{save_as}"})
    if status != 200:
        raise SystemExit(f"HTTP {status} for {url}")
    if release != EXPECTED_RELEASE:
        raise SystemExit(f"UniProt served release {release!r}, expected {EXPECTED_RELEASE}: {url}")
    return body


def bsd_sum(path):
    checksum, size = 0, 0
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with open(path, "x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.relative_to(REPO_ROOT)} ({len(rows)} rows)")


def parse_fasta(text):
    seqs, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    cache_base = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(os.environ.get("LOCALAPPDATA", Path.home())) / "bioinformatics-portfolio-cache"
    existing = [p for p in CURATED.values() if p.exists()] + ([RAW] if RAW.exists() else [])
    if existing:
        raise SystemExit("Refusing to overwrite existing outputs: " + ", ".join(map(str, existing)))
    RAW.mkdir(parents=True)

    targets = read_csv(TARGET_TABLE)
    tid = [t["ensembl_gene_id"] for t in targets]
    tinfo = {t["ensembl_gene_id"]: t for t in targets}
    xrefs = read_csv(STEP2_RAW / "ensembl97_xrefs.csv")
    txs = read_csv(STEP2_RAW / "ensembl97_transcripts.csv")
    tx2gene = {t["transcript_stable_id"]: t["gene_stable_id"] for t in txs}
    canon_tl = {t["gene_stable_id"]: t["translation_stable_id"] for t in txs if t["is_canonical"] == "True"}
    tl_version = {t["translation_stable_id"]: t["translation_version"] for t in txs if t["translation_stable_id"]}
    checks = []

    def check(name, ok, detail=""):
        checks.append({"check": name, "result": "PASS" if ok else "FAIL", "detail": detail})
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" ({detail})" if detail else ""))

    # ---- 1. Ensembl 97 protein sequences (verified download) ----
    pep_dir = cache_base / "ensembl97_pep"
    pep_dir.mkdir(parents=True, exist_ok=True)
    _, h, body = http_get(f"{PEP_BASE}/CHECKSUMS")
    (RAW / "ensembl97_pep_CHECKSUMS").write_bytes(body)
    published = {m.group(3): (int(m.group(1)), int(m.group(2)))
                 for m in re.finditer(r"(\d+)\s+(\d+)\s+(\S+)", body.decode())}
    pep_path = pep_dir / PEP_FILE
    pep_lm = "(cached copy)"
    if not pep_path.exists():
        _, h, body = http_get(f"{PEP_BASE}/{PEP_FILE}")
        pep_path.write_bytes(body)
        pep_lm = h.get("Last-Modified", "")
    obs = bsd_sum(pep_path)
    check("Ensembl 97 protein FASTA matches published BSD checksum", obs == published.get(PEP_FILE),
          f"observed {obs}, published {published.get(PEP_FILE)}")
    if obs != published.get(PEP_FILE):
        raise SystemExit("Checksum mismatch for the Ensembl 97 protein FASTA")
    write_csv(RAW / "ensembl97_pep_source.csv",
              [{"file": PEP_FILE, "url": f"{PEP_BASE}/{PEP_FILE}", "last_modified": pep_lm,
                "retrieved_utc": now_utc(), "bytes": pep_path.stat().st_size,
                "bsd_sum_published": "%d %d" % published[PEP_FILE], "bsd_sum_observed": "%d %d" % obs,
                "sha256": hashlib.sha256(pep_path.read_bytes()).hexdigest()}],
              # No local cache path is recorded: it would expose the local computer's directory layout.
              ["file", "url", "last_modified", "retrieved_utc", "bytes", "bsd_sum_published",
               "bsd_sum_observed", "sha256"])
    ens_seq = {}  # ENSP -> (version, gene, sequence)
    with gzip.open(pep_path, "rt") as f:
        for header, seq in parse_fasta(f.read()).items():
            m = re.match(r"(ENSP\d+)\.(\d+) .*gene:(ENSG\d+)\.\d+", header)
            if m and m.group(3) in tid:
                ens_seq[m.group(1)] = (m.group(2), m.group(3), seq)
    vers_ok = all(ens_seq.get(p, ("?",))[0] == v for p, v in tl_version.items())
    check("Ensembl 97 FASTA has every target translation with the same version as the Step 2 transcript table",
          vers_ok and set(ens_seq) == set(tl_version), f"{len(ens_seq)} translations")

    # ---- 2. Candidates from explicit cross-references ----
    cand = defaultdict(lambda: {"ens97": defaultdict(set), "uniprot_xref": False, "reviewed_hint": None})
    excluded_gn = defaultdict(set)
    for x in xrefs:
        g = x["gene_stable_id"]
        if x["external_db"] in ("Uniprot/SWISSPROT", "Uniprot/SPTREMBL") and x["object_type"] == "Translation":
            c = cand[(g, x["primary_id"])]
            c["ens97"][x["external_db"]].add(f'{x["object_stable_id"]}({x["info_type"]})')
            c["reviewed_hint"] = x["external_db"] == "Uniprot/SWISSPROT"
        elif x["external_db"] == "Uniprot_gn":
            excluded_gn[g].add(x["primary_id"])

    search_hits = {}
    for g in tid:
        body = uniprot_get("uniprotkb/search",
                           {"query": f"(xref:ensembl-{g}) AND (organism_id:{HUMAN})",
                            "fields": "accession,reviewed,protein_name,gene_names,organism_id,xref_ensembl",
                            "format": "json", "size": 500},
                           f"searches/ensembl_xref_{g}.json")
        res = json.loads(body)["results"]
        search_hits[g] = res
        for r in res:
            c = cand[(g, r["primaryAccession"])]
            c["uniprot_xref"] = True
            c["reviewed_hint"] = r["entryType"].startswith("UniProtKB reviewed")

    # ---- 3. Fetch entries: reviewed for all genes; unreviewed only where no reviewed candidate ----
    has_reviewed = {g: any(v["reviewed_hint"] for (gg, _), v in cand.items() if gg == g) for g in tid}
    to_fetch = sorted({a for (g, a), v in cand.items() if v["reviewed_hint"] or not has_reviewed[g]})
    entries, iso_seqs = {}, {}
    for acc in to_fetch:
        entries[acc] = json.loads(uniprot_get(f"uniprotkb/{acc}.json", None, f"entries/{acc}.json"))
        if entries[acc].get("entryType") != "Inactive":
            fa = uniprot_get("uniprotkb/stream", {"query": f"accession:{acc}", "format": "fasta",
                                                  "includeIsoform": "true"}, f"entries/{acc}_isoforms.fasta")
            iso_seqs[acc] = {h.split("|")[1]: s for h, s in parse_fasta(fa.decode()).items()}

    # UniParc records for canonical proteins of genes without a reviewed candidate (documentary evidence)
    uniparc = {}
    for g in tid:
        if has_reviewed[g]:
            continue
        upi = next((x["primary_id"] for x in xrefs if x["external_db"] == "UniParc" and x["gene_stable_id"] == g
                    and x["object_stable_id"] == canon_tl.get(g)), None)
        if upi:
            d = json.loads(uniprot_get(f"uniparc/{upi}", {"format": "json"}, f"uniparc/{upi}.json"))
            uniparc[g] = (upi, sorted({(r["database"], r["id"]) for r in d.get("uniParcCrossReferences", [])
                                       if r.get("database", "").startswith("UniProtKB") and r.get("active")
                                       and r.get("organism", {}).get("taxonId") == HUMAN}))

    # ---- 4. Evaluate every gene-accession pair ----
    def ens_links(e):
        out = []
        for r in e.get("uniProtKBCrossReferences", []):
            if r["database"] == "Ensembl":
                p = {q["key"]: q["value"] for q in r.get("properties", [])}
                out.append((r["id"], p.get("ProteinId", ""), p.get("GeneId", ""), r.get("isoformId", "")))
        return out

    mapping = []
    for (g, acc), c in sorted(cand.items(), key=lambda kv: (tid.index(kv[0][0]), kv[0][1])):
        e = entries.get(acc)
        row = {"ensembl_gene_id": g, "target_symbol": tinfo[g]["target_symbol"], "uniprot_accession": acc,
               "ensembl97_xref": "; ".join(f"{db}: {', '.join(sorted(v))}" for db, v in sorted(c["ens97"].items())) or "none",
               "uniprot_2026_04_xref_to_gene": c["uniprot_xref"], "status_in_2026_04": "", "uniprot_gene_names": "",
               "sequence_identical_ensembl97_proteins": "", "decision": "", "reason": ""}
        if e is None:
            row.update(status_in_2026_04="not fetched", decision="not used",
                       reason="unreviewed; a reviewed (Swiss-Prot) record is available for this gene")
        elif e.get("entryType") == "Inactive":
            row.update(status_in_2026_04="inactive: " + json.dumps(e.get("inactiveReason", "")),
                       decision="not used", reason="entry inactive in 2026_04")
        else:
            row["status_in_2026_04"] = e["entryType"]
            row["uniprot_gene_names"] = "; ".join(x.get("geneName", {}).get("value", "") for x in e.get("genes", []))
            ident = sorted(f"{p}.{v[0]}={iso}" for p, v in ens_seq.items() if v[1] == g
                           for iso, s in iso_seqs.get(acc, {}).items() if s == v[2])
            row["sequence_identical_ensembl97_proteins"] = "; ".join(ident) or "none"
            if e["organism"]["taxonId"] != HUMAN:
                row.update(decision="not used", reason="not a human record")
            elif not e["entryType"].startswith("UniProtKB reviewed") and has_reviewed[g]:
                row.update(decision="not used", reason="unreviewed; a reviewed (Swiss-Prot) record is available for this gene")
            else:
                row.update(decision="used", reason="reviewed record with explicit cross-reference"
                           if e["entryType"].startswith("UniProtKB reviewed")
                           else "no reviewed record for this gene; unreviewed record with explicit cross-reference")
        mapping.append(row)
    for g, accs in excluded_gn.items():
        for a in sorted(accs):
            if not any(m["ensembl_gene_id"] == g and m["uniprot_accession"] == a for m in mapping):
                mapping.append({"ensembl_gene_id": g, "target_symbol": tinfo[g]["target_symbol"], "uniprot_accession": a,
                                "ensembl97_xref": "Uniprot_gn only (gene-name based)", "uniprot_2026_04_xref_to_gene": False,
                                "status_in_2026_04": "not fetched", "decision": "excluded",
                                "reason": "only a gene-name-based link; no explicit protein or gene cross-reference"})

    # Cross-gene transcript links: an entry fetched for gene A that cross-references a transcript of target gene B
    for acc, e in entries.items():
        for enst, ensp, geneid, iso in ens_links(e):
            tx = enst.split(".")[0]
            owner = tx2gene.get(tx)
            brought_by = {g for (g, a) in cand if a == acc}
            if owner and owner not in brought_by:
                ident = [f"{p}.{v[0]}={i}" for p, v in ens_seq.items() if v[1] == owner
                         for i, s in iso_seqs.get(acc, {}).items() if s == v[2]]
                mapping.append({"ensembl_gene_id": owner, "target_symbol": tinfo[owner]["target_symbol"],
                                "uniprot_accession": acc,
                                "ensembl97_xref": "none for this gene",
                                "uniprot_2026_04_xref_to_gene": False,
                                "status_in_2026_04": e["entryType"],
                                "uniprot_gene_names": "; ".join(x.get("geneName", {}).get("value", "") for x in e.get("genes", [])),
                                "sequence_identical_ensembl97_proteins": "; ".join(sorted(ident)) or "none",
                                "decision": "ambiguous - not transferred",
                                "reason": f"UniProt 2026_04 links {iso or acc} to {enst} (protein {ensp}) but assigns it to gene "
                                          f"{geneid}; in Ensembl 97 {tx} belongs to {owner}. Ensembl 97 made no link "
                                          f"from {owner} to {acc}."})
    # de-duplicate ambiguous rows (one per gene-accession, reasons merged)
    merged = {}
    for m in mapping:
        k = (m["ensembl_gene_id"], m["uniprot_accession"], m["decision"])
        if k in merged and m["decision"].startswith("ambiguous"):
            merged[k]["reason"] += " | " + m["reason"]
        else:
            merged.setdefault(k, m)
    mapping = sorted(merged.values(), key=lambda m: (tid.index(m["ensembl_gene_id"]), m["decision"] != "used", m["uniprot_accession"]))

    # ---- 5. Curated functional annotation ----
    def primary_for(g):
        used = [m for m in mapping if m["ensembl_gene_id"] == g and m["decision"] == "used"]
        if not used:
            return None, []
        rev = [m for m in used if m["status_in_2026_04"].startswith("UniProtKB reviewed")]
        pool = rev or used
        canon = canon_tl.get(g, "")
        pool.sort(key=lambda m: (canon not in m["sequence_identical_ensembl97_proteins"],
                                 canon not in m["ensembl97_xref"], not m["uniprot_2026_04_xref_to_gene"],
                                 m["uniprot_accession"]))
        return pool[0], used

    func_rows, go_rows, ref_rows = [], [], []
    for g in tid:
        t = tinfo[g]
        prim, used = primary_for(g)
        base = {"comparison": t["comparison"], "target_symbol": t["target_symbol"], "record_role": t["record_role"],
                "ensembl_gene_id": g, "versioned_gene_id": t["versioned_gene_id"], "gene_biotype": t["gene_biotype"]}
        if prim is None:
            base.update({"uniprot_accession": "", "review_status": "no UniProtKB record mapped",
                         "missing_fields": "all UniProt fields (no explicit UniProt cross-reference)"})
            func_rows.append(base)
            continue
        acc = prim["uniprot_accession"]
        e = entries[acc]
        pdsc = e.get("proteinDescription", {})
        name = pdsc.get("recommendedName", {}).get("fullName", {}).get("value")
        name_type = "recommended name"
        if not name:
            subs = pdsc.get("submissionNames", [])
            name = subs[0]["fullName"]["value"] if subs else ""
            name_type = "submitted name" if subs else ""
        comments = e.get("comments", [])
        fn = [tx for c in comments if c["commentType"] == "FUNCTION" for tx in c.get("texts", [])]
        sl = [loc for c in comments if c["commentType"] == "SUBCELLULAR LOCATION" for loc in c.get("subcellularLocations", [])]
        fn_ev = [ev for tx in fn for ev in tx.get("evidences", [])]
        sl_ev = [ev for loc in sl for ev in loc.get("location", {}).get("evidences", [])]
        go = [r for r in e.get("uniProtKBCrossReferences", []) if r["database"] == "GO"]
        aspects = defaultdict(int)
        for r in go:
            p = {q["key"]: q["value"] for q in r.get("properties", [])}
            aspect = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}[p["GoTerm"][0]]
            aspects[aspect] += 1
            go_rows.append({"ensembl_gene_id": g, "target_symbol": t["target_symbol"], "uniprot_accession": acc,
                            "go_id": r["id"], "aspect": aspect, "go_term": p["GoTerm"][2:],
                            "go_evidence": p.get("GoEvidenceType", ""),
                            "evidence_sources": "; ".join(f'{ev.get("source", "")}:{ev.get("id", "")}' for ev in r.get("evidences", []))})
        titles = {}
        for ref in e.get("references", []):
            cit = ref.get("citation", {})
            for cx in cit.get("citationCrossReferences", []):
                if cx.get("database") == "PubMed":
                    titles[cx["id"]] = cit.get("title", "")
        for field, evs in (("function", fn_ev), ("subcellular_location", sl_ev)):
            for pmid in sorted({ev["id"] for ev in evs if ev.get("source") == "PubMed"}, key=int):
                ref_rows.append({"ensembl_gene_id": g, "target_symbol": t["target_symbol"], "uniprot_accession": acc,
                                 "supports": field, "pubmed_id": pmid, "title": titles.get(pmid, "(title not in entry)"),
                                 "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        missing = [n for n, v in (("function", fn), ("subcellular_location", sl), ("GO terms", go)) if not v]
        audit = e.get("entryAudit", {})
        base.update({
            "uniprot_accession": acc, "review_status": e["entryType"],
            "other_used_accessions": "; ".join(m["uniprot_accession"] for m in used if m["uniprot_accession"] != acc),
            "mapping_evidence": f'Ensembl 97 xref: {prim["ensembl97_xref"]}; UniProt 2026_04 xref to gene: '
                                f'{prim["uniprot_2026_04_xref_to_gene"]}; identical sequences: {prim["sequence_identical_ensembl97_proteins"]}',
            "protein_name": name, "protein_name_type": name_type,
            "uniprot_gene_names": prim["uniprot_gene_names"],
            "uniprot_gene_name_differs_from_ensembl97_symbol": t["ensembl_symbol"] not in prim["uniprot_gene_names"].split("; "),
            "protein_existence": e.get("proteinExistence", ""), "annotation_score": e.get("annotationScore", ""),
            "function": " ".join(tx["value"] for tx in fn),
            "function_evidence_codes": "; ".join(sorted({ev["evidenceCode"] for ev in fn_ev})),
            "function_pubmed_ids": "; ".join(sorted({ev["id"] for ev in fn_ev if ev.get("source") == "PubMed"}, key=int)),
            "subcellular_location": "; ".join(dict.fromkeys(loc["location"]["value"] for loc in sl)),
            "subcellular_location_evidence_codes": "; ".join(sorted({ev["evidenceCode"] for ev in sl_ev})),
            "go_molecular_function_n": aspects["molecular_function"], "go_biological_process_n": aspects["biological_process"],
            "go_cellular_component_n": aspects["cellular_component"],
            "entry_version": audit.get("entryVersion", ""), "last_annotation_update": audit.get("lastAnnotationUpdateDate", ""),
            "sequence_length": e.get("sequence", {}).get("length", ""),
            "uniprot_url": f"https://rest.uniprot.org/uniprotkb/{acc}.json",
            "missing_fields": "; ".join(missing) or "none",
        })
        func_rows.append(base)

    # ---- 6. Validation ----
    check("Every UniProt response came from release 2026_04",
          all(r["x_uniprot_release"] == EXPECTED_RELEASE for r in request_log), f"{len(request_log)} responses")
    active = [e for e in entries.values() if e.get("entryType") != "Inactive"]
    check("Every active fetched entry is a human (taxon 9606) record",
          all(e["organism"]["taxonId"] == HUMAN for e in active), f"{len(active)} entries")
    used = [m for m in mapping if m["decision"] == "used"]
    check("Every used accession has an explicit cross-reference (Ensembl 97 or UniProt 2026_04)",
          all(m["ensembl97_xref"] != "none" or m["uniprot_2026_04_xref_to_gene"] for m in used), f"{len(used)} used pairs")
    check("No accession is used on the basis of a gene-name-only link",
          not any(m["decision"] == "used" and m["ensembl97_xref"].startswith("Uniprot_gn") for m in mapping))
    check("Reviewed records preferred: no unreviewed record used for a gene that has a reviewed one",
          all(not has_reviewed[m["ensembl_gene_id"]] or m["status_in_2026_04"].startswith("UniProtKB reviewed") for m in used))
    for g in tid:
        p, _ = primary_for(g)
        if p:
            canon = canon_tl.get(g, "")
            check(f"{tinfo[g]['target_symbol']} {g}: primary record {p['uniprot_accession']} has a protein sequence identical "
                  f"to an Ensembl 97 protein of this gene",
                  p["sequence_identical_ensembl97_proteins"] not in ("", "none"),
                  f"canonical {canon} identical: {canon in p['sequence_identical_ensembl97_proteins']}")
    matr3_primary = next((r for r in func_rows if r["ensembl_gene_id"] == "ENSG00000280987"), {})
    check("Reviewed MATR3 record P43243 is not assigned to ENSG00000280987",
          matr3_primary.get("uniprot_accession") != "P43243"
          and not any(m["ensembl_gene_id"] == "ENSG00000280987" and m["uniprot_accession"] == "P43243"
                      and m["decision"] == "used" for m in mapping))
    for g in ("ENSG00000237506", "ENSG00000127589"):
        r = next(r for r in func_rows if r["ensembl_gene_id"] == g)
        check(f"Pseudogene {tinfo[g]['target_symbol']}: no protein function assigned",
              not r.get("function") and not r.get("uniprot_accession"), r["review_status"])
    check("No pre-existing output was overwritten (all outputs created new)", True, "files opened in exclusive-create mode")

    # ---- 7. Write ----
    write_csv(RAW / "request_log.csv", request_log, list(request_log[0]))
    write_csv(RAW / "uniparc_evidence.csv",
              [{"ensembl_gene_id": g, "canonical_translation": canon_tl[g], "uniparc_id": upi,
                "active_human_uniprotkb_records_with_identical_sequence": "; ".join(f"{db}:{a}" for db, a in recs)}
               for g, (upi, recs) in uniparc.items()],
              ["ensembl_gene_id", "canonical_translation", "uniparc_id", "active_human_uniprotkb_records_with_identical_sequence"])
    write_csv(CURATED["mapping"], mapping, ["ensembl_gene_id", "target_symbol", "uniprot_accession", "decision", "reason",
                                            "status_in_2026_04", "ensembl97_xref", "uniprot_2026_04_xref_to_gene",
                                            "sequence_identical_ensembl97_proteins", "uniprot_gene_names"])
    func_fields = ["comparison", "target_symbol", "record_role", "ensembl_gene_id", "versioned_gene_id", "gene_biotype",
                   "uniprot_accession", "review_status", "other_used_accessions", "mapping_evidence", "protein_name",
                   "protein_name_type", "uniprot_gene_names", "uniprot_gene_name_differs_from_ensembl97_symbol",
                   "protein_existence", "annotation_score", "function", "function_evidence_codes", "function_pubmed_ids",
                   "subcellular_location", "subcellular_location_evidence_codes", "go_molecular_function_n",
                   "go_biological_process_n", "go_cellular_component_n", "entry_version", "last_annotation_update",
                   "sequence_length", "uniprot_url", "missing_fields"]
    write_csv(CURATED["function"], func_rows, func_fields)
    write_csv(CURATED["go"], go_rows, ["ensembl_gene_id", "target_symbol", "uniprot_accession", "go_id", "aspect",
                                       "go_term", "go_evidence", "evidence_sources"])
    write_csv(CURATED["refs"], ref_rows, ["ensembl_gene_id", "target_symbol", "uniprot_accession", "supports",
                                          "pubmed_id", "title", "url"])
    write_csv(CURATED["validation"], checks, ["check", "result", "detail"])
    n_fail = sum(c["result"] == "FAIL" for c in checks)
    print(f"Validation: {len(checks) - n_fail}/{len(checks)} PASS")


if __name__ == "__main__":
    main()
