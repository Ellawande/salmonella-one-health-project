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
import re
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

# Shared AMR classification lives in amr_common.py.  The path insert makes
# the import work whether this file is run from the repo root or from scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_common import INTRINSIC_MARKERS, is_acquired  # noqa: E402
import tarfile

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
# Database hosts (mgc.ac.cn among them) answer 403 to the default urllib
# user-agent, so every outbound request identifies itself as a browser.
USER_AGENT = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")

VFDB_URL = "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
PLASMIDFINDER_TARBALL = ("https://bitbucket.org/genomicepidemiology/"
                         "plasmidfinder_db/get/master.tar.gz")
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
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
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


# PlasmidFinder headers are <replicon>_<cluster>_<description>_<accession>, e.g.
# IncC_1__JN157804, so the replicon type is the leading token up to the cluster
# number. Entries with no cluster field (ColpEC648__CP008718) keep the whole
# header as the name. Verified against all 10 replicon types in the published
# results.
REPLICON_NAME_RE = re.compile(r"^(?P<name>.+?)_(?P<cluster>\d+)_")


def replicon_name(header: str) -> str:
    m = REPLICON_NAME_RE.match(header)
    return m.group("name") if m else header


def fetch_plasmidfinder(cache: Path) -> list[tuple[str, str]]:
    """PlasmidFinder replicon sequences.

    The database moved from GitHub to Bitbucket: the old GitHub contents
    listing now returns 404. We therefore take the whole repository tarball in
    one request rather than listing files and fetching them one by one.
    """
    archive = cache / "plasmidfinder_db.tar.gz"
    if not archive.exists():
        print("  fetching PlasmidFinder database ...")
        try:
            req = urllib.request.Request(PLASMIDFINDER_TARBALL,
                                         headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r:
                archive.write_bytes(r.read())
        except Exception as exc:
            raise SystemExit(
                "ERROR: could not download the PlasmidFinder database.\n"
                f"  {PLASMIDFINDER_TARBALL}\n  -> {exc}\n"
                f"  Download it manually and save it as {archive}"
            ) from exc

    seqs: list[tuple[str, str]] = []
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.endswith(".fsa"):
                continue
            text = tar.extractfile(member).read().decode("utf-8", errors="replace")
            # The replicon type comes from the HEADER (see replicon_name above),
            # not the file name: the database now groups sequences into a handful
            # of files such as enterobacteriales.fsa.
            for header, s in parse_fasta(text):
                seqs.append((replicon_name(header), s))
    return seqs


# --- annotation helpers -------------------------------------------------------
# CARD's model FASTA carries only sequences; the fields that say WHAT a hit
# confers live in aro_index.tsv. VFDB encodes its own annotation in the FASTA
# header. Neither was being read, so the published tables could not be
# regenerated. These helpers recover both, plus host_group from the metadata.

HIT_COLUMNS = ("strain", "query", "database", "gene", "identity_%", "coverage_%",
               "category", "product", "drug", "mech", "family")

AMR_COLUMNS = HIT_COLUMNS + ("acquired",)

REPLICON_COLUMNS = ("strain", "host_group", "replicon", "identity_%", "coverage_%")


def resolve_card_dir(cache: Path) -> Path:
    """CARD files: the project cache if populated, else the Omicsboard cache."""
    if (cache / "aro_index.tsv").exists():
        return cache
    alt = Path.home() / ".omicsboard" / "amrdb" / "card"
    if (alt / "aro_index.tsv").exists():
        return alt
    return cache


def card_protein_fasta(cache: Path) -> Path:
    """CARD ships the protein-homolog model under either of two filenames."""
    for d in (cache, resolve_card_dir(cache)):
        for name in ("card_protein.fasta",
                     "protein_fasta_protein_homolog_model.fasta"):
            if (d / name).exists():
                return d / name
    return cache / "card_protein.fasta"


def load_card_meta(cache: Path) -> dict:
    """ARO accession -> (drug class, resistance mechanism, AMR gene family).

    Returns {} rather than raising when aro_index.tsv is missing, so the
    script still runs -- with empty annotation columns and a loud warning.
    """
    idx = resolve_card_dir(cache) / "aro_index.tsv"
    if not idx.exists():
        return {}
    want = ("ARO Accession", "Drug Class", "Resistance Mechanism", "AMR Gene Family")
    meta = {}
    with idx.open(encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if not all(w in header for w in want):
            return {}
        i_acc, i_drug, i_mech, i_fam = (header.index(w) for w in want)
        widest = max(i_acc, i_drug, i_mech, i_fam)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= widest:
                continue
            meta[f[i_acc].strip()] = (f[i_drug].strip(), f[i_mech].strip(),
                                      f[i_fam].strip())
    return meta


def card_fields(reference: str, card_meta: dict) -> dict:
    """Gene name + CARD annotation for one hit, from its FASTA header."""
    parts = reference.split("|")
    gene = parts[3].split()[0] if len(parts) > 3 else reference.split()[0]
    m = re.search(r"ARO:\d+", reference)
    drug, mech, family = card_meta.get(m.group(0), ("", "", "")) if m else ("", "", "")
    return {"gene": gene, "category": "", "product": "",
            "drug": drug, "mech": mech, "family": family}


def vfdb_fields(reference: str) -> dict:
    """VFDB annotates inside the FASTA header, e.g.
    VFG037176(gb|WP_001081735) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [organism]
    -> gene 'gb|WP_001081735', product '(plc1) phospholipase C', category 'Exotoxin'.
    """
    m = re.search(r"\((gb\|[A-Za-z0-9._]+)\)", reference)
    gene = m.group(1) if m else reference.split()[0]
    rest = reference[m.end():].strip() if m else ""
    product = rest.split("[")[0].strip()
    category = ""
    cm = re.search(r"\[([^\]]+)\]", rest)
    if cm:
        c = cm.group(1).strip()
        if " - " in c:
            c = c.split(" - ")[-1]
        category = re.sub(r"\s*\(VFC\d+\)\s*$", "", c).strip()
    return {"gene": gene, "category": category, "product": product,
            "drug": "", "mech": "", "family": ""}


def load_host_groups(metadata: Path, out: Path) -> dict:
    """asm_acc -> host_group, so replicons can be compared across host groups."""
    for p in (metadata, out / "isolates_selected.csv"):
        if not p.exists():
            continue
        with p.open(encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split(",")
            if "asm_acc" not in header or "host_group" not in header:
                continue
            i_acc, i_host = header.index("asm_acc"), header.index("host_group")
            widest = max(i_acc, i_host)
            groups = {}
            for line in fh:
                f = line.rstrip("\n").split(",")
                if len(f) > widest:
                    groups[f[i_acc].strip()] = f[i_host].strip()
            return groups
    return {}


def fetch_card(cache: Path) -> list[tuple[str, str]]:
    """
    CARD is license-gated: downloading requires accepting the CARD terms.
    We try the public endpoint; if it is unavailable we stop with clear
    instructions rather than quietly returning an empty database (which would
    silently under-report resistance -- the worst possible failure mode here).
    """
    prot = card_protein_fasta(cache)
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
    """Prefilter by shared k-mers, then verify the survivors by alignment.

    The gate is |qset & rset| / |rset| >= prefilter. Testing it with a set
    intersection for every query/reference pair is O(queries x refs) -- about
    28 million intersections for one genome, which is ~6 minutes per genome
    and over three hours for all eighteen.

    An inverted index (k-mer -> the references carrying it) computes the SAME
    shared count, because a pair can only pass the gate if it shares at least
    one k-mer; pairs sharing none are zero and never get visited. Candidate
    references are visited in their original ascending order, so the output is
    identical to the naive version, only faster.
    """
    ref_sets = [(label, seq, kmer_set(seq, k)) for label, seq in refs]
    ref_len = [len(rset) for _, _, rset in ref_sets]
    hits: list[dict] = []

    # prefilter <= 0 admits every pair, so no index could prune anything
    if prefilter <= 0:
        for qname, qseq in queries:
            qset = kmer_set(qseq, k)
            if not qset:
                continue
            for label, rseq, rset in ref_sets:
                if not rset or len(qset & rset) / len(rset) < prefilter:
                    continue
                identity, coverage = align_identity_coverage(qseq, rseq, scores)
                if identity >= min_identity and coverage >= min_coverage:
                    hits.append({"reference": label,
                                 "identity": round(identity, 3),
                                 "coverage": round(coverage, 3)})
        return hits

    index: dict[str, list[int]] = {}
    for i, (_, _, rset) in enumerate(ref_sets):
        for km in rset:
            index.setdefault(km, []).append(i)

    for qname, qseq in queries:
        qset = kmer_set(qseq, k)
        if not qset:
            continue
        shared: dict[int, int] = {}
        for km in qset:
            for i in index.get(km, ()):
                shared[i] = shared.get(i, 0) + 1
        for i in sorted(shared):
            if shared[i] / ref_len[i] < prefilter:
                continue
            label, rseq, _ = ref_sets[i]
            identity, coverage = align_identity_coverage(qseq, rseq, scores)
            if identity >= min_identity and coverage >= min_coverage:
                hits.append({"reference": label,
                             "identity": round(identity, 3),
                             "coverage": round(coverage, 3)})
    return hits

def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Screen genomes for AMR, virulence and replicons.")
    ap.add_argument("--genomes", type=Path, default=root / "data" / "raw")
    ap.add_argument("--out", type=Path, default=root / "results")
    ap.add_argument("--cache", type=Path, default=root / "refdb")
    ap.add_argument("--metadata", type=Path,
                    default=root / "data" / "metadata" / "isolates.csv",
                    help="isolate metadata carrying host_group")
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

    card_meta = load_card_meta(args.cache)
    host_groups = load_host_groups(args.metadata, args.out)
    if card and not card_meta:
        print("  WARNING: aro_index.tsv not found -- CARD drug/mechanism/"
              "family columns will be empty.")
    print(f"  CARD         : {len(card)} protein models")
    print(f"  VFDB         : {len(vfdb)} virulence factors")
    print(f"  PlasmidFinder: {len(replicons)} replicon alleles")
    print("-" * 74)

    amr_rows, vir_rows, rep_rows = [], [], []
    for i, fasta in enumerate(genomes, 1):
        strain = fasta.stem

        # Per-genome checkpoints: screening is ~4 min per genome, so a long run
        # must be resumable rather than restartable. Completed genomes load back
        # from disk instead of being recomputed.
        ckpt = root / "data" / "amr_checkpoints"
        ck_paths = {t: ckpt / f"{strain}.{t}.pkl" for t in ("amr", "vir", "rep")}
        if all(q.exists() for q in ck_paths.values()):
            amr_rows += pd.read_pickle(ck_paths["amr"]).to_dict("records")
            vir_rows += pd.read_pickle(ck_paths["vir"]).to_dict("records")
            rep_rows += pd.read_pickle(ck_paths["rep"]).to_dict("records")
            print(f"genome {i} of {len(genomes)}: {strain} -- resumed from checkpoint")
            continue
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
            f = card_fields(h["reference"], card_meta)
            amr_rows.append({"strain": strain, "query": h["reference"],
                             "database": "CARD", "gene": f["gene"],
                             "identity_%": h["identity"], "coverage_%": h["coverage"],
                             "category": f["category"], "product": f["product"],
                             "drug": f["drug"], "mech": f["mech"], "family": f["family"]})
        for h in vir_hits:
            f = vfdb_fields(h["reference"])
            vir_rows.append({"strain": strain, "query": h["reference"],
                             "database": "VFDB", "gene": f["gene"],
                             "identity_%": h["identity"], "coverage_%": h["coverage"],
                             "category": f["category"], "product": f["product"],
                             "drug": f["drug"], "mech": f["mech"], "family": f["family"]})
        for h in rep_hits:
            rep_rows.append({"strain": strain, "host_group": host_groups.get(strain, ""),
                             "replicon": h["reference"],
                             "identity_%": h["identity"], "coverage_%": h["coverage"]})


        # write this genome's rows to disk immediately so a crash costs one genome
        ckpt.mkdir(parents=True, exist_ok=True)
        for t, rows in (("amr", amr_rows), ("vir", vir_rows), ("rep", rep_rows)):
            sub = [r for r in rows if r["strain"] == strain]
            pd.DataFrame(sub).to_pickle(ck_paths[t])
    amr_df = pd.DataFrame(amr_rows).reindex(columns=HIT_COLUMNS)
    # CARD-derived classification: acquired = transferable / plasmid-borne,
    # intrinsic = chromosomal. Written into the table so consumers (script 06)
    # read one source of truth instead of each re-deriving the flag privately.
    amr_df["acquired"] = amr_df["gene"].map(is_acquired)
    amr_df = amr_df.reindex(columns=AMR_COLUMNS)
    amr_path = args.out / "amr_genes.csv"
    amr_df.to_csv(amr_path, index=False)

    vir_df = pd.DataFrame(vir_rows).reindex(columns=HIT_COLUMNS)
    vir_path = args.out / "virulence_genes.csv"
    vir_df.to_csv(vir_path, index=False)

    rep_df = pd.DataFrame(rep_rows).reindex(columns=REPLICON_COLUMNS)
    rep_path = args.out / "plasmid_replicons_all.csv"
    rep_df.to_csv(rep_path, index=False)

    print("-" * 74)
    print(f"CARD AMR hits: {len(amr_df)} rows -> {amr_path}")
    if not amr_df.empty:
        flags = amr_df["acquired"]
        print(f"  distinct acquired AMR genes : {amr_df.loc[flags, 'gene'].nunique()}")
        print(f"  distinct intrinsic AMR genes: {amr_df.loc[~flags, 'gene'].nunique()}")
        top = Counter(amr_df.loc[flags, "gene"]).most_common(5)
        if top:
            print("  most widespread acquired genes:")
            for g, n in top:
                print(f"     {g} -- {n} isolates")

    print(f"VFDB hits    : {len(vir_df)} rows -> {vir_path}")
    if not vir_df.empty:
        print(f"  distinct virulence genes: {vir_df['gene'].nunique()}")

    print(f"Replicons    : {len(rep_df)} rows -> {rep_path}")
    if not rep_df.empty:
        print(f"  {rep_df['strain'].nunique()}/{len(genomes)} isolates carry a plasmid replicon")
        for rec, n in Counter(rep_df["replicon"]).most_common(5):
            print(f"     {rec} -- {n} isolates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
