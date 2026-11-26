"""
Phase 8, Step 2 validation: re-checks the Ensembl release 97 annotation retrieval.

Every check below is performed when the script runs; the CSV records the expected
value, the actual value, PASS/FAIL and the time it was performed. Nothing is
modified except the output CSV, which is never overwritten.

Inputs (read only)
  results/tables/phase8_raw/ensembl97_*.csv, results/tables/phase8_target_gene_annotation.csv
  results/tables/gene_annotation.csv (Phase 5)
  data/reference/Homo_sapiens.GRCh38.97.gtf.gz and data/reference/ANNOTATION_SOURCE.md
  The downloaded release-97 database files in the local cache (see 08_gene_annotation_retrieval.py)
  Two read-only web requests: Ensembl's published CHECKSUMS file, and the retired archive URL.
Output
  results/tables/phase8_annotation_validation.csv

Run from anywhere:
    python data/raw/notebooks/08_validate_gene_annotation.py [cache_dir]
"""

import csv
import gzip
import hashlib
import os
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
T = REPO_ROOT / "results/tables"
RAW = T / "phase8_raw"
OUT = T / "phase8_annotation_validation.csv"
GTF = REPO_ROOT / "data/reference/Homo_sapiens.GRCh38.97.gtf.gz"
BASE_URL = "https://ftp.ensembl.org/pub/release-97/mysql/homo_sapiens_core_97_38"
ARCHIVE_URL = "https://jul2019.archive.ensembl.org/"
USER_AGENT = "bioinformatics-portfolio-annotation/1.0"
DOCUMENTED_COUNTS = {"ensembl97_gene_records.csv": 9, "ensembl97_transcripts.csv": 71,
                     "ensembl97_attributes.csv": 660, "ensembl97_xrefs.csv": 1183}

rows = []


def now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def record(category, check, expected, actual, method):
    ok = str(expected) == str(actual)
    rows.append({"check_id": f"S2-{len(rows) + 1:02d}", "category": category, "check": check,
                 "expected": expected, "actual": actual, "result": "PASS" if ok else "FAIL",
                 "method": method, "performed_utc": now()})


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def bsd_sum(path):
    checksum, size = 0, 0
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + byte) & 0xFFFF
    return f"{checksum} {(size + 1023) // 1024}"


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.geturl(), r.read()


def main():
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite {OUT}")
    cache = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(os.environ.get("LOCALAPPDATA", Path.home())) / "bioinformatics-portfolio-cache/ensembl97_core_db"

    # A. Downloaded database files: SHA-256 vs recorded; BSD sum vs Ensembl's published CHECKSUMS
    src = read(RAW / "ensembl97_source_files.csv")
    published = {}
    for line in (cache / "CHECKSUMS").read_text().splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+(\S+)$", line)
        if m:
            published[m.group(3)] = f"{int(m.group(1))} {int(m.group(2))}"
    for s in src:
        p = cache / s["file"]
        if s["file"] == "CHECKSUMS":
            record("Ensembl 97 database files", "CHECKSUMS: cached copy matches recorded SHA-256",
                   s["sha256"], sha256(p) if p.exists() else "file missing", "SHA-256 of cached file")
            continue
        actual_sha = sha256(p) if p.exists() else "file missing"
        actual_bsd = bsd_sum(p) if p.exists() else "file missing"
        record("Ensembl 97 database files",
               f"{s['file']}: SHA-256 matches record and BSD sum matches Ensembl's published CHECKSUMS",
               f"sha256 {s['sha256'][:16]}...; bsd {published.get(s['file'])}",
               f"sha256 {actual_sha[:16] if actual_sha != s['sha256'] else s['sha256'][:16]}...; bsd {actual_bsd}"
               if actual_sha == s["sha256"] else f"sha256 MISMATCH {actual_sha[:16]}...; bsd {actual_bsd}",
               "SHA-256 and BSD sum recomputed from cached file")

    # B. Published CHECKSUMS on the Ensembl server still identical to the copy used
    _, body = fetch(f"{BASE_URL}/CHECKSUMS")
    record("Ensembl 97 database files", "Ensembl's CHECKSUMS file on the server is identical to the copy used",
           sha256(cache / "CHECKSUMS"), hashlib.sha256(body).hexdigest(), "Fresh download compared by SHA-256")

    # C. Archive website retired (reason for using the database files)
    final_url, _ = fetch(ARCHIVE_URL)
    record("Source", "jul2019.archive.ensembl.org redirects to Ensembl's archives help page (site retired)",
           "https://www.ensembl.org/help/articles/archives", final_url, "HTTP request, final URL after redirects")

    # D. Reference GTF checksum
    doc = (REPO_ROOT / "data/reference/ANNOTATION_SOURCE.md").read_text(encoding="utf-8")
    documented = re.search(r"SHA-256 \| `([0-9a-f]{64})`", doc).group(1)
    record("Reference GTF", "Homo_sapiens.GRCh38.97.gtf.gz matches the SHA-256 in ANNOTATION_SOURCE.md",
           documented, sha256(GTF), "SHA-256 of GTF")

    # E-F. Gene and transcript records vs GTF
    genes = read(RAW / "ensembl97_gene_records.csv")
    txs = read(RAW / "ensembl97_transcripts.csv")
    ids = {g["stable_id"] for g in genes}
    gtf_genes, gtf_tx = {}, {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] not in ("gene", "transcript"):
                continue
            gid = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
            if gid not in ids:
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', p[8]))
            strand = "1" if p[6] == "+" else "-1"
            if p[2] == "gene":
                gtf_genes[gid] = f'{a.get("gene_version")}|{a.get("gene_name")}|{a.get("gene_biotype")}|{p[0]}|{p[3]}|{p[4]}|{strand}'
            else:
                gtf_tx[a["transcript_id"]] = f'{gid}|{a.get("transcript_version")}|{a.get("transcript_biotype")}'
    for g in genes:
        exp = f'{g["version"]}|{g["display_xref_label"]}|{g["biotype"]}|{g["chromosome"]}|{g["seq_region_start"]}|{g["seq_region_end"]}|{g["seq_region_strand"]}'
        record("Gene records vs GTF", f'{g["stable_id"]}: version, name, biotype, chromosome, start, end, strand',
               exp, gtf_genes.get(g["stable_id"], "not in GTF"), "Field-by-field comparison with GTF gene record")
    tx_ok = sum(1 for t in txs if gtf_tx.get(t["transcript_stable_id"]) ==
                f'{t["gene_stable_id"]}|{t["transcript_version"]}|{t["biotype"]}')
    record("Transcript records vs GTF", "Every extracted transcript matches the GTF (gene, version, biotype)",
           len(txs), tx_ok, "Count of transcripts whose gene ID, version and biotype equal the GTF")
    record("Transcript records vs GTF", "No GTF transcript of the nine genes is missing from the extract",
           len(gtf_tx), len({t["transcript_stable_id"] for t in txs} & set(gtf_tx)), "Set comparison")

    # G. Raw extract row counts vs documented counts
    for name, n in DOCUMENTED_COUNTS.items():
        record("Raw extracts", f"{name}: row count equals the count documented in PHASE8_ANNOTATION_PROVENANCE.md",
               n, len(read(RAW / name)), "Row count")

    # H. Curated table consistent with raw extracts
    cur = read(T / "phase8_target_gene_annotation.csv")
    attrs = read(RAW / "ensembl97_attributes.csv")
    xrefs = read(RAW / "ensembl97_xrefs.csv")
    gmap = {g["stable_id"]: g for g in genes}
    record("Curated table", "Nine rows, one per unique target gene ID, all present in the raw gene records",
           "9 unique; all in raw", f'{len({c["ensembl_gene_id"] for c in cur})} unique; '
           f'{"all in raw" if all(c["ensembl_gene_id"] in gmap for c in cur) else "some missing"}', "Set comparison")
    mism = []
    for c in cur:
        g = gmap[c["ensembl_gene_id"]]
        gx = [t for t in txs if t["gene_stable_id"] == g["stable_id"]]
        canon = next(t for t in gx if t["is_canonical"] == "True")
        rt = sorted({a["object_stable_id"] for a in attrs if a["gene_stable_id"] == g["stable_id"] and a["attrib_code"] == "readthrough_tra"})
        hg = sorted({x["primary_id"] for x in xrefs if x["gene_stable_id"] == g["stable_id"] and x["external_db"] == "HGNC" and x["object_type"] == "Gene"})
        checks = {
            "versioned_gene_id": f'{g["stable_id"]}.{g["version"]}', "gene_biotype": g["biotype"],
            "chromosome": g["chromosome"], "start": g["seq_region_start"], "end": g["seq_region_end"],
            "n_transcripts": str(len(gx)), "canonical_transcript": f'{canon["transcript_stable_id"]}.{canon["transcript_version"]}',
            "readthrough_attribute_transcripts": ";".join(rt) or "none", "hgnc_id": ";".join(hg),
            "ensembl_symbol": g["display_xref_label"],
        }
        mism += [f'{c["ensembl_gene_id"]}:{k}' for k, v in checks.items() if c[k] != v]
    record("Curated table", "Curated fields (ID/version, biotype, location, transcripts, canonical transcript, "
           "readthrough flags, HGNC ID, symbol) equal the raw extracts", "0 mismatches",
           f"{len(mism)} mismatches" + (f" ({'; '.join(mism)})" if mism else ""), "Recomputed from raw extracts")

    # I. Cross-phase consistency with the Phase 5 annotation table
    p5 = {r["gene_id"]: r for r in read(T / "gene_annotation.csv")}
    diffs = [c["ensembl_gene_id"] for c in cur if p5[c["ensembl_gene_id"]]["gene_name"] != c["ensembl_symbol"]
             or p5[c["ensembl_gene_id"]]["gene_biotype"] != c["gene_biotype"]
             or p5[c["ensembl_gene_id"]]["n_transcripts"] != c["n_transcripts"]]
    record("Cross-phase", "Symbol, biotype and transcript count agree with Phase 5 gene_annotation.csv",
           "0 differences", f"{len(diffs)} differences" + (f" ({', '.join(diffs)})" if diffs else ""), "Field comparison")

    # J. No local computer paths in Step 2 outputs
    files = [RAW / n for n in list(DOCUMENTED_COUNTS) + ["ensembl97_source_files.csv"]] + [T / "phase8_target_gene_annotation.csv"]
    pat = re.compile(r"[A-Za-z]:[\\/]+Users[\\/]|/home/\w+|/Users/\w+|AppData")
    hits = [f.name for f in files if pat.search(f.read_text(encoding="utf-8"))]
    record("Privacy", "No local computer paths in Step 2 output files", "0 files", f"{len(hits)} files" +
           (f" ({', '.join(hits)})" if hits else ""), "Pattern search for user-directory paths")

    with open(OUT, "x", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["check_id", "category", "check", "expected", "actual", "result",
                                          "method", "performed_utc"])
        w.writeheader()
        w.writerows(rows)
    n_fail = sum(r["result"] == "FAIL" for r in rows)
    for r in rows:
        print(f'  [{r["result"]}] {r["check_id"]} {r["check"]}' + ("" if r["result"] == "PASS" else f' | expected {r["expected"]} | actual {r["actual"]}'))
    print(f"wrote {OUT.relative_to(REPO_ROOT)}: {len(rows) - n_fail}/{len(rows)} PASS")


if __name__ == "__main__":
    main()
