#!/usr/bin/env python3
"""
01_fetch_isolates.py  --  Stage 1: data acquisition.

Downloads the Salmonella enterica genome assemblies listed in
data/metadata/isolates.csv and writes them to data/raw/.

WHY THIS SCRIPT EXISTS
    The 18 genome FASTA files are NOT committed to this repository. They total
    roughly 85 MB and every one of them is fully regenerable from a single
    accession number. Committing them would bloat the repo and freeze a
    snapshot; running this script instead re-derives them from the
    authoritative source on demand. This is what makes the analysis
    reproducible from scratch.

METHOD
    NCBI Datasets v2 REST API -- the supported, stable download route:

        GET https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<ACC>/download
            ?include_annotation_type=GENOME_FASTA

    Assemblies are requested by their GCA_ accession. Do NOT try to guess an
    NCBI FTP directory path for an assembly: the folder segment is not
    derivable from the accession and the URL will 404.

HOW THE ISOLATES WERE CHOSEN
    Source: NCBI Pathogen Detection -- a curated surveillance database that
    carries real isolation metadata (country, host, collection date) alongside
    NCBI's own pre-computed AMR genotype for each isolate. Selection aimed for
    six isolates per host group (Human, Poultry, Cattle) so the three groups
    are comparable in size. See data/metadata/isolates.csv for the final set
    and the README for the full rationale.

USAGE
    python scripts/01_fetch_isolates.py
    python scripts/01_fetch_isolates.py --force       # re-download everything
    python scripts/01_fetch_isolates.py --metadata path/to/isolates.csv

OUTPUT
    data/raw/<accession>.fna      one FASTA per genome
    data/raw/fetch_report.csv     accession, file, size_bytes, n_contigs, status

DEPENDENCIES
    pip install requests pandas
"""
from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

import pandas as pd
import requests

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# Candidate column names for the accession field. The metadata table has been
# through a few revisions, so we detect the column rather than hard-coding it.
ACCESSION_CANDIDATES = ["asm_acc", "assembly_accession", "assembly", "accession",
                        "asm_accn", "gca", "Assembly Accession"]


def repo_root() -> Path:
    """The repository root -- the parent of this script's directory."""
    return Path(__file__).resolve().parent.parent


def find_accession_column(df: pd.DataFrame) -> str:
    """Locate the column holding the assembly accession, whatever it is called."""
    for name in ACCESSION_CANDIDATES:
        if name in df.columns:
            return name
    # fall back: any column whose values look like a GCA_/GCF_ accession
    for name in df.columns:
        sample = df[name].dropna().astype(str).head(20)
        if sample.str.match(r"^GC[AF]_\d+\.\d+$").any():
            return name
    raise SystemExit(
        "Could not find an accession column in the metadata table.\n"
        f"Columns present: {list(df.columns)}\n"
        "Expected one of: " + ", ".join(ACCESSION_CANDIDATES)
    )


def fetch_one(acc: str, dest: Path) -> dict:
    """Download one assembly, extract its FASTA, and report what arrived."""
    url = f"{API}/genome/accession/{acc}/download"
    params = {"include_annotation_type": "GENOME_FASTA"}

    try:
        r = requests.get(url, params=params, timeout=600)
        r.raise_for_status()
    except requests.HTTPError as e:
        return {"status": f"HTTP {r.status_code}: {e}", "n_contigs": None, "size_bytes": None}
    except requests.RequestException as e:
        return {"status": f"network error: {e}", "n_contigs": None, "size_bytes": None}

    # NCBI returns a zip; the FASTA sits at
    # ncbi_dataset/data/<ACC>/<ACC>_ASM<...>_genomic.fna
    try:
        zf = zipfile.ZipFile(io.BytesIO(r.content))
    except zipfile.BadZipFile:
        return {"status": "response was not a zip archive", "n_contigs": None, "size_bytes": None}

    fna_members = [n for n in zf.namelist() if n.endswith(".fna")]
    if not fna_members:
        return {"status": "no .fna found inside the archive", "n_contigs": None, "size_bytes": None}

    seq = zf.read(fna_members[0])
    dest.write_bytes(seq)

    # count contigs by counting FASTA headers
    n_contigs = sum(1 for line in seq.splitlines() if line.startswith(b">"))

    return {"status": "ok", "n_contigs": n_contigs, "size_bytes": len(seq)}


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Fetch Salmonella genome assemblies from NCBI.")
    ap.add_argument("--metadata", type=Path, default=root / "data" / "metadata" / "isolates.csv",
                    help="CSV listing the isolates (default: data/metadata/isolates.csv)")
    ap.add_argument("--outdir", type=Path, default=root / "data" / "raw",
                    help="where to write the FASTA files (default: data/raw)")
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the FASTA already exists")
    args = ap.parse_args()

    if not args.metadata.exists():
        print(f"ERROR: metadata table not found: {args.metadata}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.metadata)
    acc_col = find_accession_column(df)
    accessions = df[acc_col].dropna().astype(str).str.strip().tolist()

    print(f"Metadata : {args.metadata}")
    print(f"Accession column: '{acc_col}'")
    print(f"Isolates : {len(accessions)}")
    print(f"Output   : {args.outdir}")
    print("-" * 72)

    report = []
    for i, acc in enumerate(accessions, 1):
        dest = args.outdir / f"{acc}.fna"

        if dest.exists() and not args.force:
            size = dest.stat().st_size
            n = sum(1 for line in dest.open("rb") if line.startswith(b">"))
            print(f"[{i:>2}/{len(accessions)}] {acc}  cached  "
                  f"{size/1e6:5.1f} MB  {n} contigs")
            report.append({"accession": acc, "file": dest.name, "size_bytes": size,
                           "n_contigs": n, "status": "cached"})
            continue

        print(f"[{i:>2}/{len(accessions)}] {acc}  downloading ...", end="", flush=True)
        res = fetch_one(acc, dest)
        if res["status"] == "ok":
            print(f"\r[{i:>2}/{len(accessions)}] {acc}  ok       "
                  f"{res['size_bytes']/1e6:5.1f} MB  {res['n_contigs']} contigs")
        else:
            print(f"\r[{i:>2}/{len(accessions)}] {acc}  FAILED   {res['status']}")

        report.append({"accession": acc, "file": dest.name,
                       "size_bytes": res["size_bytes"], "n_contigs": res["n_contigs"],
                       "status": res["status"]})

    rep = pd.DataFrame(report)
    rep_path = args.outdir / "fetch_report.csv"
    rep.to_csv(rep_path, index=False)

    print("-" * 72)
    ok = (rep["status"].isin(["ok", "cached"])).sum()
    print(f"{ok}/{len(rep)} genomes available")
    print(f"Total size: {rep['size_bytes'].fillna(0).sum()/1e6:.1f} MB")
    print(f"Report written: {rep_path}")

    if ok < len(rep):
        print("\nSome downloads failed. NCBI occasionally rate-limits bulk requests;")
        print("simply re-run this script -- completed files are skipped.")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
