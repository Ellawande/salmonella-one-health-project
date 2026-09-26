#!/usr/bin/env python3
"""
02_amr_and_plasmids.py  --  Stage 2: resistance, virulence and mobility.

Screens every genome for:
  * antimicrobial-resistance genes   -- CARD   (amino acid homology)
  * virulence factors                -- VFDB   (amino acid homology)
  * plasmid replicons                -- PlasmidFinder (nucleotide homology)

WHY THIS SCRIPT EXISTS
    A hand-written gene list would be worthless here. Reporting "this isolate
    carries blaCTX-M" is a factual claim that must come from a curated
    database, with the database version cited -- guessing resistance genes from
    memory produces plausible-looking results that are simply wrong for a real
    clinical isolate. So this script screens against the real, published
    databases and records %identity and %coverage for every hit.

THE CENTRAL DISTINCTION
    CARD contains BOTH:
      * INTRINSIC genes   -- chromosomal efflux pumps and their regulators,
                             present in essentially every strain. Their
                             presence says little about a particular isolate.
      * ACQUIRED genes    -- plasmid- or transposon-borne (blaCTX-M, blaKPC,
                             aac, tet, sul, dfrA ...). These are the
                             clinically actionable ones.
    The `acquired` flag written by this script drives that split, and every
    downstream analysis depends on it.

METHOD
    Gene calling  : pyrodigal (Prodigal), trained per genome.
    AMR / VF      : k-mer prefilter, then verified local amino-acid alignment.
    Replicons     : k-mer prefilter on the NUCLEOTIDE genome, then verified
                    local nucleotide alignment. Replicons live in the
                    replication origin, which is largely intergenic -- so this
                    step screens the raw genome, not the predicted proteins.
    Thresholds    : CARD/VFDB     >=90% identity, >=60% coverage (ResFinder-style)
                    PlasmidFinder >=95% identity, >=60% coverage (CGE convention)

HONEST SCOPE
    * Homology screening detects resistance-ASSOCIATED genes. Point-mutation
      resistance (gyrA, rpoB, pmrB) is a separate analysis requiring variant
      calling against a known-mutation catalogue.
    * Replicon detection proves a plasmid of a given Inc type is present. It
      does NOT prove which genes sit on that plasmid -- that needs assembly
      plus plasmid binning (MOB-suite, plasmidSPAdes) or long reads.
    * Screening-grade. For regulatory work confirm with AMRFinderPlus or RGI.

DATABASE ACCESS
    VFDB          : open download (http://www.mgc.ac.cn/VFs/)
    PlasmidFinder : open (github.com/genomicepidemiology/plasmidfinder_db)
    CARD          : license-gated. This script attempts the public download;
                    if blocked it prints manual instructions rather than
                    silently returning an empty database, which would
                    under-report resistance.

USAGE
    python scripts/02_amr_and_plasmids.py
    python scripts/02_amr_and_plasmids.py --genomes data/raw --out results/

OUTPUT
    results/amr_hits.csv            one row per AMR/VF gene hit per genome
    results/plasmid_replicons.csv   one row per replicon hit per genome

DEPENDENCIES
    pip install pyrodigal biopython pandas
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------------
# Thresholds and scoring -- published conventions for each database
# ----------------------------------------------------------------------------
AMR_MIN_IDENTITY = 0.90      # CARD / VFDB : amino-acid identity
AMR_MIN_COVERAGE = 0.60      # fraction of the reference model covered
REP_MIN_IDENTITY = 0.95      # PlasmidFinder : nucleotide identity
REP_MIN_COVERAGE = 0.60

PROTEIN_KMER = 5             # amino-acid k-mer for the prefilter
NUCLEOTIDE_KMER = 12         # nucleotide k-mer -- longer, since DNA is 4-letter

# local-alignment scoring
PROTEIN_SCORES = dict(match=2, mismatch=-1, gap=-3)
NUCLEOTIDE_SCORES = dict(match=2, mismatch=-2, gap=-4)


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


# ----------------------------------------------------------------------------
# Reference database acquisition
# ----------------------------------------------------------------------------
VFDB_URL = "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
PLASMIDFINDER_API = ("https://api.github.com/repos/genomicepidemiology/"
                     "plasmidfinder_db/contents")
PLASMIDFINDER_RAW = ("https://raw.githubusercontent.com/genomicepidemiology/"
                     "plasmidfinder_db/master")
CARD_INFO = "https://card.mcmaster.ca/latest/data"


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Minimal FASTA parser -> [(description, sequence), ...]."""
    out, name, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out.append((name, "".join(chunks)))
    return [(n, s.upper()) for n, s in out if s]


def fetch_text(url: str, timeout: int = 300) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        data = r.read()
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("ascii", errors="replace")


def fetch_vfdb(cache: Path) -> list[tuple[str, str]]:
    target = cache / "VFDB_setA_pro.fas"
    if not target.exists():
        print("  fetching VFDB (virulence factors) ...")
        target.write_text(fetch_text(VFDB_URL), encoding="utf-8")
    return parse_fasta(target.read_text(encoding="utf-8", errors="replace"))


def fetch_plasmidfinder(cache: Path) -> list[tuple[str, str]]:
    listing = cache / "plasmidfinder_files.json"
    if not listing.exists():
        print("  fetching PlasmidFinder file listing ...")
        with urllib.request.urlopen(PLASMIDFINDER_API, timeout=120) as r:
            payload = json.load(r)
        names = [e["name"] for e in payload if e["name"].endswith(".fsa")]
        listing.write_text(json.dumps(names), encoding="utf-8")
    names = json.loads(listing.read_text(encoding="utf-8"))

    seqs: list[tuple[str, str]] = []
    for name in names:
        local = cache / name
        if not local.exists():
            local.write_text(fetch_text(f"{PLASMIDFINDER_RAW}/{name}"),
                             encoding="utf-8")
        # the FILE name IS the replicon type (IncFII.fsa -> IncFII)
        replicon = name[:-4]
        for _, s in parse_fasta(local.read_text(encoding="utf-8", errors="replace")):
            seqs.append((replicon, s))
    return seqs


def fetch_card(cache: Path) -> list[tuple[str, str]]:
    """
    CARD is license-gated: downloading requires accepting the CARD terms.
    We try the public endpoint; if it is unavailable we stop with clear
    instructions rather than quietly returning an empty database (which would
    silently under-report resistance -- the worst possible failure mode here).
    """
    prot = cache / "card_protein.fasta"
    if prot.exists():
        return parse_fasta(prot.read_text(encoding="utf-8", errors="replace"))

    print(f"  fetching CARD ... (license-gated at {CARD_INFO})")
    try:
        blob = urllib.request.urlopen(CARD_INFO, timeout=300).read()
    except Exception as e:
        raise SystemExit(
            f"\nCARD could not be downloaded automatically "
            f"({type(e).__name__}: {e}).\n\n"
            "CARD requires you to accept its license before downloading:\n"
            "  1. Open https://card.mcmaster.ca/download\n"
            "  2. Accept the terms and download the latest 'card-data' archive\n"
            "  3. Extract it and copy two files into:\n"
            f"       {cache}\n"
            "     - protein_fasta_protein_homolog_model.fasta -> card_protein.fasta\n"
            "     - card.json                                 -> card_meta.json\n"
            "  4. Re-run this script.\n\n"
            "CARD's license forbids redistribution, which is why these files\n"
            "are not included in this repository.\n"
        )
    raise SystemExit(
        f"CARD returned {len(blob)} bytes, but the archive layout varies by\n"
        "release. Please download 'card-data' manually from\n"
        "https://card.mcmaster.ca/download and place the files in\n"
        f"{cache} as described above."
    )


# ----------------------------------------------------------------------------
# Gene calling
# ----------------------------------------------------------------------------
def call_genes(fasta: Path) -> list[tuple[str, str]]:
    """Predict protein-coding genes with pyrodigal (Prodigal, single mode)."""
    import pyrodigal

    genome = "".join(line.strip() for line in fasta.read_text().splitlines()
                     if not line.startswith(">"))
    finder = pyrodigal.GeneFinder(meta=True)
    genes = finder.find_genes(genome.encode())

    proteins = []
    for i, gene in enumerate(genes, 1):
        aa = gene.translate()
        if aa and len(aa) >= 30:                 # ignore very short ORFs
            proteins.append((f"{fasta.stem}|gene{i}", aa))
    return proteins


def read_genome(fasta: Path) -> str:
    return "".join(line.strip() for line in fasta.read_text().splitlines()
                   if not line.startswith(">")).upper()


# ----------------------------------------------------------------------------
# Homology screening: k-mer prefilter, then verified local alignment
# ----------------------------------------------------------------------------
def kmer_set(seq: str, k: int) -> set[str]:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def align_identity_coverage(query: str, ref: str, scores: dict) -> tuple[float, float]:
    """
    Local alignment of query against ref.
    Returns (identity over aligned columns, fraction of REF covered).
    """
    from Bio.Align import PairwiseAligner

    al = PairwiseAligner()
    al.mode = "local"
    al.match_score = scores["match"]
    al.mismatch_score = scores["mismatch"]
    al.gap_score = scores["gap"]
    aln = al.align(ref, query)[0]

    a, b = aln[0], aln[1]
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    aligned = sum(1 for x, y in zip(a, b) if x != "-" and y != "-")
    identity = matches / aligned if aligned else 0.0
    return identity, (aligned / len(ref) if ref else 0.0)


def screen(queries: list[tuple[str, str]],
           refs: list[tuple[str, str]],
           k: int, scores: dict,
           min_identity: float, min_coverage: float,
           prefilter: float = 0.20) -> list[dict]:
    """Prefilter by shared k-mers, then verify the survivors by alignment."""
    ref_sets = [(label, seq, kmer_set(seq, k)) for label, seq in refs]
    hits: list[dict] = []

    for qname, qseq in queries:
        qset = kmer_set(qseq, k)
        if not qset:
            continue
        for label, rseq, rset in ref_sets:
            if not rset:
                continue
            # --- cheap gate: too few shared k-mers to be worth aligning
            if len(qset & rset) / len(rset) < prefilter:
                continue
            # --- verify with a real alignment
            identity, coverage = align_identity_coverage(qseq, rseq, scores)
            if identity >= min_identity and coverage >= min_coverage:
                hits.append({"reference": label,
                             "identity": round(identity, 3),
                             "coverage": round(coverage, 3)})
    return hits


# ----------------------------------------------------------------------------
# Which CARD genes are ACQUIRED rather than intrinsic?
# ----------------------------------------------------------------------------
# Curated from the CARD classification and the ResFinder/AMRFinder literature.
# Intrinsic genes are chromosomal in Salmonella and say little about a single
# isolate; acquired genes are the transferable, clinically actionable ones.
INTRINSIC_MARKERS = (
    "aac(6')-iaa", "aac(6')-iay", "aac(6')-iz", "aac(6')-iic",
    "mdtk", "mdfa", "acra", "acrb", "mara", "ramA", "soxs", "emrd", "emry",
    "mdsabc", "tolc", "baea", "baer", "kdpe", "phop", "phoq",
    "cpxa", "cpxr", "crp", "hfq", "roba", "sdia", "gadx", "gadw",
    "ampr", "amps", "lpxc", "pmrf", "pmre", "ugd", "arna", "arnb",
    "baca", "epta", "eptb", "mcr-9",
)


def is_acquired(label: str) -> bool:
    """True if a CARD hit looks plasmid-/transposon-borne rather than intrinsic."""
    low = label.lower()
    return not any(marker in low for marker in INTRINSIC_MARKERS)


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Screen genomes for AMR, virulence and replicons.")
    ap.add_argument("--genomes", type=Path, default=root / "data" / "raw")
    ap.add_argument("--out", type=Path, default=root / "results")
    ap.add_argument("--cache", type=Path, default=root / "refdb")
    ap.add_argument("--skip-card", action="store_true",
                    help="run only VFDB + PlasmidFinder")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)

    genomes = sorted(args.genomes.glob("*.fna"))
    if not genomes:
        print(f"ERROR: no .fna files in {args.genomes}", file=sys.stderr)
        print("Run scripts/01_fetch_isolates.py first.", file=sys.stderr)
        return 1

    print(f"Genomes: {len(genomes)}")
    print("Loading reference databases ...")
    card = [] if args.skip_card else fetch_card(args.cache)
    vfdb = fetch_vfdb(args.cache)
    replicons = fetch_plasmidfinder(args.cache)
    print(f"  CARD         : {len(card)} protein models")
    print(f"  VFDB         : {len(vfdb)} virulence factors")
    print(f"  PlasmidFinder: {len(replicons)} replicon alleles")
    print("-" * 74)

    amr_rows, rep_rows = [], []
    for i, fasta in enumerate(genomes, 1):
        strain = fasta.stem
        proteins = call_genes(fasta)

        amr_hits = screen(proteins, card, PROTEIN_KMER, PROTEIN_SCORES,
                          AMR_MIN_IDENTITY, AMR_MIN_COVERAGE)
        vir_hits = screen(proteins, vfdb, PROTEIN_KMER, PROTEIN_SCORES,
                          AMR_MIN_IDENTITY, AMR_MIN_COVERAGE)

        # replicons are nucleotide -- screen the genome, not the proteins
        genome_nt = read_genome(fasta)
        rep_hits = screen([(strain, genome_nt)], replicons,
                          NUCLEOTIDE_KMER, NUCLEOTIDE_SCORES,
                          REP_MIN_IDENTITY, REP_MIN_COVERAGE, prefilter=0.30)

        print(f"[{i:>2}/{len(genomes)}] {strain}  {len(proteins):>5} genes  "
              f"{len(amr_hits):>2} AMR  {len(vir_hits):>2} VF  "
              f"{len(rep_hits):>2} replicons")

        for h in amr_hits:
            amr_rows.append({"strain": strain, "database": "CARD",
                             "gene": h["reference"], "identity": h["identity"],
                             "coverage": h["coverage"],
                             "acquired": is_acquired(h["reference"])})
        for h in vir_hits:
            amr_rows.append({"strain": strain, "database": "VFDB",
                             "gene": h["reference"], "identity": h["identity"],
                             "coverage": h["coverage"], "acquired": None})
        for h in rep_hits:
            rep_rows.append({"strain": strain, "replicon": h["reference"],
                             "identity": h["identity"], "coverage": h["coverage"]})

    amr_df = pd.DataFrame(amr_rows)
    amr_path = args.out / "amr_hits.csv"
    amr_df.to_csv(amr_path, index=False)

    rep_df = pd.DataFrame(rep_rows)
    rep_path = args.out / "plasmid_replicons.csv"
    rep_df.to_csv(rep_path, index=False)

    print("-" * 74)
    print(f"AMR/VF hits : {len(amr_df):>5} rows -> {amr_path}")
    if not amr_df.empty:
        acq = amr_df[amr_df["acquired"] == True]["gene"].nunique()
        intr = amr_df[amr_df["acquired"] == False]["gene"].nunique()
        print(f"  distinct acquired AMR genes : {acq}")
        print(f"  distinct intrinsic AMR genes: {intr}")
        top = Counter(amr_df[amr_df["acquired"] == True]["gene"]).most_common(5)
        if top:
            print("  most widespread acquired genes:")
            for g, n in top:
                print(f"     {g:<40} {n} isolates")

    print(f"\nReplicon hits: {len(rep_df):>5} rows -> {rep_path}")
    if not rep_df.empty:
        strains_with = rep_df["strain"].nunique()
        print(f"  {strains_with}/{len(genomes)} isolates carry >=1 plasmid replicon")
        for rec, n in Counter(rep_df["replicon"]).most_common(5):
            print(f"     {rec:<40} {n} isolates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
