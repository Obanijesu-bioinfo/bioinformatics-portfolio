"""
Phase 8, Step 2: official Ensembl release 97 annotation for the target genes.

Source
  The Ensembl 97 website archive (jul2019.archive.ensembl.org) has been retired:
  every URL on it, including its BioMart service, now redirects to
  https://www.ensembl.org/help/articles/archives. This script therefore reads the
  same release from Ensembl's own release-97 database dump, which is the data
  the archive site served:
      https://ftp.ensembl.org/pub/release-97/mysql/homo_sapiens_core_97_38/
  No current Ensembl record is used.

What it does
  1. Downloads the needed core-database tables into a local cache outside the
     repository (re-used if present) and verifies each file against Ensembl's
     published CHECKSUMS (BSD `sum`); also records SHA-256.
  2. Reads column names from the release's schema file (homo_sapiens_core_97_38.sql.gz).
  3. Extracts, for each target gene ID: the gene record, its transcripts and
     translations, gene and transcript attributes, and all external references
     (xrefs) attached to the gene, its transcripts and its translations.
  4. Checks the gene IDs, versions, biotypes and coordinates against the local
     release-97 GTF and stops on any mismatch.
  5. Writes the raw extracted rows (results/tables/phase8_raw/) separately from
     a curated table of annotation facts (results/tables/phase8_target_gene_annotation.csv).
     The curated table holds official annotation only; it contains no
     biological interpretation.

It refuses to overwrite any existing output file.

Run from anywhere:
    python data/raw/notebooks/08_gene_annotation_retrieval.py [cache_dir]
Default cache: %LOCALAPPDATA%/bioinformatics-portfolio-cache/ensembl97_core_db
"""

import csv
import gzip
import hashlib
import os
import re
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
GTF_FILE = REPO_ROOT / "data/reference/Homo_sapiens.GRCh38.97.gtf.gz"
RAW_DIR = REPO_ROOT / "results/tables/phase8_raw"
CURATED_FILE = REPO_ROOT / "results/tables/phase8_target_gene_annotation.csv"

BASE_URL = "https://ftp.ensembl.org/pub/release-97/mysql/homo_sapiens_core_97_38"
SCHEMA_FILE = "homo_sapiens_core_97_38.sql.gz"
TABLES = ["gene", "transcript", "translation", "seq_region", "attrib_type",
          "gene_attrib", "transcript_attrib", "external_db", "xref",
          "object_xref", "external_synonym"]
USER_AGENT = "bioinformatics-portfolio-annotation/1.0"

# Targets: (comparison, target symbol, Ensembl gene ID, record role)
TARGETS = [
    ("C1", "RPSAP15", "ENSG00000237506", "target"),
    ("C1", "TUBBP1", "ENSG00000127589", "target"),
    ("C1", "INHBE", "ENSG00000139269", "target"),
    ("C1", "AC069288.1", "ENSG00000286192", "target"),
    ("C2", "SFT2D2", "ENSG00000213064", "target"),
    ("C2", "MATR3", "ENSG00000280987", "target (primary MATR3 record)"),
    ("C2", "MATR3", "ENSG00000015479", "overlapping MATR3 gene model (reported separately)"),
    ("C2", "STEAP2", "ENSG00000157214", "target"),
    ("C2", "EOGT", "ENSG00000163378", "target"),
]
TARGET_IDS = [t[2] for t in TARGETS]

OUT_FILES = {
    "genes": RAW_DIR / "ensembl97_gene_records.csv",
    "transcripts": RAW_DIR / "ensembl97_transcripts.csv",
    "attributes": RAW_DIR / "ensembl97_attributes.csv",
    "xrefs": RAW_DIR / "ensembl97_xrefs.csv",
    "sources": RAW_DIR / "ensembl97_source_files.csv",
}


def now_utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def bsd_sum(path):
    """BSD `sum` (the algorithm of Ensembl's CHECKSUMS files): (checksum, 1K blocks)."""
    checksum, size = 0, 0
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as r:
        last_modified = r.headers.get("Last-Modified", "")
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as out:
            while chunk := r.read(1 << 20):
                out.write(chunk)
    tmp.replace(dest)
    return last_modified


def split_mysql_row(line):
    """Split one mysqldump --tab row; handles backslash escapes and \\N (NULL)."""
    line = line.rstrip("\n")
    if "\\" not in line:
        return line.split("\t")
    fields, cur, i = [], [], 0
    while i < len(line):
        ch = line[i]
        if ch == "\\" and i + 1 < len(line):
            nxt = line[i + 1]
            if nxt == "N" and not cur and (i + 2 == len(line) or line[i + 2] == "\t"):
                cur.append("\\N")
            else:
                cur.append({"t": "\t", "n": "\n", "0": "\0"}.get(nxt, nxt))
            i += 2
            continue
        if ch == "\t":
            fields.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
        i += 1
    fields.append("".join(cur))
    return fields


def read_table(path, columns, keep=None):
    """Yield dict rows; `keep(row)` filters early. NULL becomes ''."""
    with gzip.open(path, "rt", encoding="utf-8", newline="\n") as f:
        pending = ""
        for line in f:
            if pending:
                line = pending + line
                pending = ""
            body = line.rstrip("\n")
            # a record continues on the next line if it ends with an odd number of backslashes
            if (len(body) - len(body.rstrip("\\"))) % 2 == 1:
                pending = line
                continue
            vals = split_mysql_row(line)
            if len(vals) != len(columns):
                raise SystemExit(f"{path.name}: expected {len(columns)} fields, got {len(vals)}")
            row = {c: ("" if v == "\\N" else v) for c, v in zip(columns, vals)}
            if keep is None or keep(row):
                yield row


def parse_schema(path):
    cols, table = {}, None
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            m = re.match(r"CREATE TABLE `(\w+)`", line)
            if m:
                table = m.group(1)
                cols[table] = []
                continue
            if table:
                m = re.match(r"\s+`(\w+)`\s", line)
                if m:
                    cols[table].append(m.group(1))
                elif line.startswith(")"):
                    table = None
    return cols


def write_csv(path, rows, fields):
    with open(path, "x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.relative_to(REPO_ROOT)} ({len(rows)} rows)")


def main():
    cache = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(os.environ.get("LOCALAPPDATA", Path.home())) / "bioinformatics-portfolio-cache/ensembl97_core_db"
    existing = [p for p in list(OUT_FILES.values()) + [CURATED_FILE] if p.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing outputs: " + ", ".join(str(p) for p in existing))
    cache.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Cache: {cache}")

    # 1. Download and verify
    checks_path = cache / "CHECKSUMS"
    checks_lm = fetch(f"{BASE_URL}/CHECKSUMS", checks_path)
    published = {}
    for line in checks_path.read_text().splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\S+)$", line)
        if m:
            published[m.group(3)] = (int(m.group(1)), int(m.group(2)))
    sources = [{"file": "CHECKSUMS", "url": f"{BASE_URL}/CHECKSUMS", "last_modified": checks_lm,
                "retrieved_utc": now_utc(), "bytes": checks_path.stat().st_size,
                "bsd_sum_published": "", "bsd_sum_observed": "", "checksum_match": "n/a (checksum list)",
                "sha256": sha256(checks_path)}]
    for name in [SCHEMA_FILE] + [f"{t}.txt.gz" for t in TABLES]:
        dest = cache / name
        if dest.exists():
            lm, when = "(cached copy)", f"cached; file time {datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc):%Y-%m-%dT%H:%M:%SZ}"
        else:
            lm = fetch(f"{BASE_URL}/{name}", dest)
            when = now_utc()
        obs = bsd_sum(dest)
        pub = published.get(name)
        ok = pub == obs
        print(f"  {name}: {dest.stat().st_size:,} bytes, BSD sum {obs} vs published {pub} -> {'OK' if ok else 'MISMATCH'}")
        sources.append({"file": name, "url": f"{BASE_URL}/{name}", "last_modified": lm, "retrieved_utc": when,
                        "bytes": dest.stat().st_size, "bsd_sum_published": f"{pub[0]} {pub[1]}" if pub else "",
                        "bsd_sum_observed": f"{obs[0]} {obs[1]}", "checksum_match": ok, "sha256": sha256(dest)})
        if not ok:
            raise SystemExit(f"Checksum mismatch for {name}; delete it from the cache and rerun.")

    # 2. Schema
    cols = parse_schema(cache / SCHEMA_FILE)
    T = lambda t: (cache / f"{t}.txt.gz", cols[t])

    # 3. Extract
    genes = {r["stable_id"]: r for r in read_table(*T("gene"), keep=lambda r: r["stable_id"] in TARGET_IDS)}
    missing = set(TARGET_IDS) - set(genes)
    if missing:
        raise SystemExit(f"Gene IDs not found in the release-97 gene table: {sorted(missing)}")
    gid2stable = {g["gene_id"]: s for s, g in genes.items()}
    seq_ids = {g["seq_region_id"] for g in genes.values()}
    seq_names = {r["seq_region_id"]: r["name"] for r in read_table(*T("seq_region"), keep=lambda r: r["seq_region_id"] in seq_ids)}
    transcripts = list(read_table(*T("transcript"), keep=lambda r: r["gene_id"] in gid2stable))
    tid2tx = {t["transcript_id"]: t for t in transcripts}
    translations = list(read_table(*T("translation"), keep=lambda r: r["transcript_id"] in tid2tx))
    tlid2tl = {t["translation_id"]: t for t in translations}
    attrib_types = {r["attrib_type_id"]: r for r in read_table(*T("attrib_type"))}
    gene_attr = list(read_table(*T("gene_attrib"), keep=lambda r: r["gene_id"] in gid2stable))
    tx_attr = list(read_table(*T("transcript_attrib"), keep=lambda r: r["transcript_id"] in tid2tx))
    ext_db = {r["external_db_id"]: r for r in read_table(*T("external_db"))}

    wanted = {("Gene", k) for k in gid2stable} | {("Transcript", k) for k in tid2tx} | {("Translation", k) for k in tlid2tl}
    oxrefs = list(read_table(*T("object_xref"), keep=lambda r: (r["ensembl_object_type"], r["ensembl_id"]) in wanted))
    display_ids = {g["display_xref_id"] for g in genes.values()} | {t["display_xref_id"] for t in transcripts}
    xref_ids = {o["xref_id"] for o in oxrefs} | display_ids
    xrefs = {r["xref_id"]: r for r in read_table(*T("xref"), keep=lambda r: r["xref_id"] in xref_ids)}
    synonyms = defaultdict(list)
    for r in read_table(*T("external_synonym"), keep=lambda r: r["xref_id"] in xref_ids):
        synonyms[r["xref_id"]].append(r["synonym"])

    # 4. Cross-check against the local GTF
    gtf = {}
    with gzip.open(GTF_FILE, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] == "gene":
                a = dict(re.findall(r'(\S+) "([^"]*)"', p[8]))
                if a["gene_id"] in TARGET_IDS:
                    gtf[a["gene_id"]] = {"version": a.get("gene_version"), "name": a.get("gene_name"),
                                         "biotype": a.get("gene_biotype"), "chr": p[0], "start": p[3],
                                         "end": p[4], "strand": p[6]}
    for s, g in genes.items():
        db = {"version": g["version"], "biotype": g["biotype"], "chr": seq_names[g["seq_region_id"]],
              "start": g["seq_region_start"], "end": g["seq_region_end"],
              "strand": "+" if g["seq_region_strand"] == "1" else "-"}
        diff = {k: (v, gtf[s][k]) for k, v in db.items() if gtf[s][k] != v}
        if diff:
            raise SystemExit(f"{s}: database dump and local GTF disagree: {diff}")
    print("  GTF cross-check: IDs, versions, biotypes and coordinates match for all", len(genes), "genes")

    # 5a. Raw outputs
    role = {t[2]: t for t in TARGETS}
    gene_rows = []
    for s in TARGET_IDS:
        g = dict(genes[s])
        dx = xrefs.get(g["display_xref_id"], {})
        g.update({"chromosome": seq_names[g["seq_region_id"]],
                  "display_xref_label": dx.get("display_label", ""),
                  "display_xref_db": ext_db.get(dx.get("external_db_id", ""), {}).get("db_name", "")})
        gene_rows.append(g)
    write_csv(OUT_FILES["genes"], gene_rows, ["stable_id", "version"] + [c for c in cols["gene"] if c not in ("stable_id", "version")]
              + ["chromosome", "display_xref_label", "display_xref_db"])

    tl_by_tx = {t["transcript_id"]: t for t in translations}
    tx_rows = []
    for t in sorted(transcripts, key=lambda t: (TARGET_IDS.index(gid2stable[t["gene_id"]]), t["stable_id"])):
        tl = tl_by_tx.get(t["transcript_id"], {})
        dx = xrefs.get(t["display_xref_id"], {})
        tx_rows.append({"gene_stable_id": gid2stable[t["gene_id"]], "transcript_stable_id": t["stable_id"],
                        "transcript_version": t["version"], "transcript_name": dx.get("display_label", ""),
                        "biotype": t["biotype"], "source": t["source"], "start": t["seq_region_start"],
                        "end": t["seq_region_end"], "strand": t["seq_region_strand"],
                        "is_canonical": t["transcript_id"] == genes[gid2stable[t["gene_id"]]]["canonical_transcript_id"],
                        "translation_stable_id": tl.get("stable_id", ""), "translation_version": tl.get("version", "")})
    write_csv(OUT_FILES["transcripts"], tx_rows, list(tx_rows[0]))

    attr_rows = []
    for a in gene_attr:
        at = attrib_types[a["attrib_type_id"]]
        attr_rows.append({"gene_stable_id": gid2stable[a["gene_id"]], "object_type": "Gene", "object_stable_id": gid2stable[a["gene_id"]],
                          "attrib_code": at["code"], "attrib_name": at["name"], "value": a["value"]})
    for a in tx_attr:
        at = attrib_types[a["attrib_type_id"]]
        t = tid2tx[a["transcript_id"]]
        attr_rows.append({"gene_stable_id": gid2stable[t["gene_id"]], "object_type": "Transcript", "object_stable_id": t["stable_id"],
                          "attrib_code": at["code"], "attrib_name": at["name"], "value": a["value"]})
    attr_rows.sort(key=lambda r: (TARGET_IDS.index(r["gene_stable_id"]), r["object_type"], r["object_stable_id"], r["attrib_code"]))
    write_csv(OUT_FILES["attributes"], attr_rows, list(attr_rows[0]))

    xref_rows = []
    for o in oxrefs:
        typ, eid = o["ensembl_object_type"], o["ensembl_id"]
        if typ == "Gene":
            gs, obj = gid2stable[eid], gid2stable[eid]
        elif typ == "Transcript":
            gs, obj = gid2stable[tid2tx[eid]["gene_id"]], tid2tx[eid]["stable_id"]
        else:
            tx = tid2tx[tlid2tl[eid]["transcript_id"]]
            gs, obj = gid2stable[tx["gene_id"]], tlid2tl[eid]["stable_id"]
        x = xrefs[o["xref_id"]]
        db = ext_db[x["external_db_id"]]
        xref_rows.append({"gene_stable_id": gs, "object_type": typ, "object_stable_id": obj,
                          "external_db": db["db_name"], "external_db_display": db["db_display_name"],
                          "primary_id": x["dbprimary_acc"], "display_label": x["display_label"],
                          "xref_version": x["version"], "description": x["description"],
                          "info_type": x["info_type"], "info_text": x["info_text"],
                          "synonyms": ";".join(sorted(synonyms.get(o["xref_id"], [])))})
    xref_rows.sort(key=lambda r: (TARGET_IDS.index(r["gene_stable_id"]), r["object_type"], r["object_stable_id"], r["external_db"], r["primary_id"]))
    write_csv(OUT_FILES["xrefs"], xref_rows, list(xref_rows[0]))
    write_csv(OUT_FILES["sources"], sources, list(sources[0]))

    # 5b. Curated annotation facts (no interpretation)
    curated = []
    for comp, sym, s, rl in TARGETS:
        g = genes[s]
        gx = [x for x in xref_rows if x["gene_stable_id"] == s]
        pick = lambda db, typ="Gene": sorted({x["primary_id"] for x in gx if x["external_db"] == db and x["object_type"] == typ})
        hgnc = [x for x in gx if x["external_db"] == "HGNC" and x["object_type"] == "Gene"]
        desc = g["description"]
        m = re.match(r"^(.*?)\s*\[Source:([^;\]]+);Acc:([^\]]+)\]$", desc)
        txs = [t for t in tx_rows if t["gene_stable_id"] == s]
        canon = next((t for t in txs if t["is_canonical"]), {})
        bt = defaultdict(int)
        for t in txs:
            bt[t["biotype"]] += 1
        rt = sorted({a["object_stable_id"] for a in attr_rows if a["gene_stable_id"] == s and "readthrough" in a["attrib_code"].lower()})
        dx = xrefs.get(g["display_xref_id"], {})
        curated.append({
            "comparison": comp, "target_symbol": sym, "record_role": rl,
            "ensembl_gene_id": s, "gene_version": g["version"], "versioned_gene_id": f"{s}.{g['version']}",
            "gtf_match": True, "ensembl_symbol": dx.get("display_label", ""),
            "symbol_source": ext_db.get(dx.get("external_db_id", ""), {}).get("db_name", ""),
            "description": m.group(1) if m else desc,
            "description_source": f"{m.group(2)}; {m.group(3)}" if m else ("" if not desc else "(no source tag)"),
            "hgnc_id": ";".join(sorted({x["primary_id"] for x in hgnc})),
            "hgnc_symbol": ";".join(sorted({x["display_label"] for x in hgnc})),
            "hgnc_synonyms": ";".join(sorted({y for x in hgnc for y in x["synonyms"].split(";") if y})),
            "entrez_gene_id": ";".join(pick("EntrezGene")),
            "uniprot_swissprot": ";".join(sorted({x["primary_id"] for x in gx if x["external_db"] == "Uniprot/SWISSPROT"})),
            "gene_biotype": g["biotype"], "gene_source": g["source"],
            "chromosome": seq_names[g["seq_region_id"]], "start": g["seq_region_start"], "end": g["seq_region_end"],
            "strand": "+" if g["seq_region_strand"] == "1" else "-",
            "n_transcripts": len(txs),
            "transcript_biotypes": "; ".join(f"{k}={v}" for k, v in sorted(bt.items())),
            "canonical_transcript": f"{canon.get('transcript_stable_id', '')}.{canon.get('transcript_version', '')}",
            "canonical_transcript_biotype": canon.get("biotype", ""),
            "readthrough_attribute_transcripts": ";".join(rt) if rt else "none",
            "annotation_source": "Ensembl 97 core database (homo_sapiens_core_97_38)",
        })
    write_csv(CURATED_FILE, curated, list(curated[0]))
    print("Done.")


if __name__ == "__main__":
    main()
