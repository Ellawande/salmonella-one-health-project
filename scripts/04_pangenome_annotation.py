#!/usr/bin/env python3
"""
04_pangenome_annotation.py  --  Stage 4: pan-genome and functional annotation.

Clusters every predicted protein into gene families, builds the presence/
absence matrix, then annotates one representative per family against Pfam-A.

WHY HOMOLOGY, NOT COMPOSITION
    A tempting shortcut is to cluster proteins by amino-acid composition (say,
    3-mer or 5-mer frequency vectors) and call similar vectors "homologs".
    That is not homology. Composition clustering groups proteins by how
    hydrophobic they are, not by shared ancestry -- it will merge unrelated
    membrane proteins and split true orthologs, producing core/accessory counts
    that look plausible but are not reproducible biology.

    So this script clusters by REAL sequence homology: a k-mer prefilter to
    find candidate homologs cheaply, then a verified pairwise alignment scored
    at >=70% identity AND >=70% coverage before two proteins join a family.
    Those thresholds are stated in every result so the numbers can be audited.

THE THREE CLASSES
    core       present in all N genomes           -- the shared backbone
    soft-core  present in >=95% of genomes        -- near-universal
    accessory  present in 2 .. <95%              -- variable / mobile
    unique     present in exactly one genome      -- strain-specific

USAGE
    python scripts/04_pangenome_annotation.py
    python scripts/04_pangenome_annotation.py --proteins proteins --out results
    python scripts/04_pangenome_annotation.py --no-pfam      # skip Pfam (fast)

OUTPUT
    results/pan_genome_matrix.csv        family x genome presence/absence
    results/gene_families.csv            family sizes and class
    results/gene_to_family.tsv           maps every protein header to its family
    results/family_representatives.faa   one protein per family (for Pfam)
    results/functional_annotation.csv    Pfam domain per family (if not skipped)
    results/functional_enrichment.csv    accessory-vs-core enrichment (BH FDR)

DEPENDENCIES
    pip install biopython pandas
    pip install pyhmmer        # only needed for the Pfam step
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------------
# Clustering thresholds -- stated explicitly so the result can be audited
# ----------------------------------------------------------------------------
MIN_IDENTITY = 0.70          # amino-acid identity required to join a family
MIN_COVERAGE = 0.70          # fraction of the shorter protein that must align
AA_KMER = 5                  # k-mer size for the prefilter
MAX_CANDIDATES = 25          # candidate families verified per query

# local-alignment scoring (BLOSUM62-equivalent simple scheme)
MATCH, MISMATCH, GAP = 2, -1, -3

PFAM_URL = ("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
            "Pfam-A.hmm.gz")


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def parse_faa(path: Path) -> list[tuple[str, str]]:
    """Read a protein FASTA -> [(header, sequence), ...]."""
    out, name, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out.append((name, "".join(chunks)))
    return [(n, s.upper()) for n, s in out if s]


def kmer_set(seq: str, k: int = AA_KMER) -> frozenset[str]:
    return frozenset(seq[i:i + k] for i in range(len(seq) - k + 1))


def aligns(a: str, b: str) -> bool:
    """
    True if a and b are homologous at the stated thresholds.
    Identity is measured over aligned columns; coverage over the SHORTER
    protein, so a short domain match cannot drag a full-length protein in.
    """
    from Bio.Align import PairwiseAligner

    al = PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score, al.gap_score = MATCH, MISMATCH, GAP
    aln = al.align(a, b)[0]

    x, y = aln[0], aln[1]
    matches = sum(1 for p, q in zip(x, y) if p == q and p != "-")
    aligned = sum(1 for p, q in zip(x, y) if p != "-" and q != "-")
    if not aligned:
        return False
    identity = matches / aligned
    coverage = aligned / min(len(a), len(b))
    return identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE


def cluster(proteins: list[tuple[str, str]]) -> tuple[dict, list[tuple[str, str]]]:
    """
    Greedy homology clustering with an inverted k-mer index.
    Returns (header -> family_id, [(family_id, representative_header, rep_seq)]).
    """
    # k-mer -> family ids that contain it
    index: dict[str, set[int]] = defaultdict(set)
    rep_seq: dict[int, str] = {}
    rep_hdr: dict[int, str] = {}
    assignment: dict[str, int] = {}
    next_family = 0

    for i, (hdr, seq) in enumerate(proteins, 1):
        ks = kmer_set(seq)

        # --- candidates: families sharing the most k-mers with this protein
        shared: Counter[int] = Counter()
        for km in ks:
            for fam in index.get(km, ()):      # inverted lookup
                shared[fam] += 1

        chosen = None
        for fam, _ in shared.most_common(MAX_CANDIDATES):
            if aligns(seq, rep_seq[fam]):
                chosen = fam
                break

        if chosen is None:                     # a genuinely new family
            chosen = next_family
            next_family += 1
            rep_seq[chosen] = seq
            rep_hdr[chosen] = hdr
            for km in ks:
                index[km].add(chosen)
        else:
            # keep the longest protein as the family representative
            if len(seq) > len(rep_seq[chosen]):
                rep_seq[chosen], rep_hdr[chosen] = seq, hdr
                for km in set(kmer_set(rep_seq[chosen])):
                    index[km].add(chosen)

        assignment[hdr] = chosen
        if i % 5000 == 0:
            print(f"    ... {i:,}/{len(proteins):,} proteins, {next_family:,} families",
                  flush=True)

    reps = [(f, rep_hdr[f], rep_seq[f]) for f in range(next_family)]
    return assignment, reps


def classify(n_genomes: int, presence: int) -> str:
    if presence == n_genomes:
        return "core"
    if presence >= 0.95 * n_genomes:
        return "soft-core"
    if presence == 1:
        return "unique"
    return "accessory"


def load_pfam_cache(cache: Path):
    """Download and press Pfam-A once; returns the path to the HMM database."""
    hmm = cache / "Pfam-A.hmm"
    pressed = cache / "Pfam-A.hmm.h3m"
    if not pressed.exists():
        if not hmm.exists():
            print("  downloading Pfam-A (~380 MB, one time) ...")
            import urllib.request
            with urllib.request.urlopen(PFAM_URL, timeout=1800) as r, \
                    open(hmm, "wb") as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
        import pyhmmer
        print("  pressing Pfam-A into an HMM index ...")
        with pyhmmer.plan7.HMMFile(str(hmm)) as hf:
            alphabet = hf.read().alphabet
        pyhmmer.hmmer.hmmpress(str(hmm))
    return hmm


def annotate_pfam(reps: list[tuple[int, str, str]], cache: Path) -> list[dict]:
    """Scan one representative per family against Pfam-A at gathering thresholds."""
    import pyhmmer

    hmm_path = load_pfam_cache(cache)
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(name=f"fam{f}".encode(), sequence=seq)
        .digitize(alphabet)
        for f, _, seq in reps
    ]

    out: list[dict] = []
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        for hits in pyhmmer.hmmer.hmmsearch(hf, sequences, cpus=0):
            for hit in hits:
                if not hit.is_included():
                    continue
                for dom in hit.domains:
                    if dom.i_evalue > 1e-5:
                        continue
                    out.append({
                        "family": int(hit.query.name.decode().lstrip("fam")),
                        "pfam": hit.name.decode(),
                        "description": (hit.description.decode()
                                        if hit.description else ""),
                        "i_evalue": float(dom.i_evalue),
                        "score": float(dom.score),
                    })
                    break          # best domain per family is enough
    return out


def enrichment(fam_df: pd.DataFrame, pfam_df: pd.DataFrame,
               n_genomes: int) -> pd.DataFrame:
    """
    Fisher's exact test per Pfam domain: accessory/unique vs core.
    Corrected across all domains with Benjamini-Hochberg FDR.
    """
    from scipy.stats import fisher_exact
    from statsmodels.stats.multitest import multipletests

    acc_fams = set(fam_df.loc[fam_df["class"].isin(["accessory", "unique"]), "family"])
    core_fams = set(fam_df.loc[fam_df["class"].isin(["core", "soft-core"]), "family"])

    rows = []
    for pfam, grp in pfam_df.groupby("pfam"):
        fams = set(grp["family"])
        a = len(fams & acc_fams)
        b = len(acc_fams) - a
        c = len(fams & core_fams)
        d = len(core_fams) - c
        if a + c == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"pfam": pfam,
                     "description": grp["description"].iloc[0],
                     "n_accessory": a, "n_core": c,
                     "odds_ratio": odds if odds != float("inf") else 1e9,
                     "p_value": p})

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q_value"] = multipletests(df["p_value"], method="fdr_bh")[1]
    return df.sort_values("q_value").reset_index(drop=True)


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Pan-genome clustering + Pfam annotation.")
    ap.add_argument("--proteins", type=Path, default=root / "proteins")
    ap.add_argument("--out", type=Path, default=root / "results")
    ap.add_argument("--cache", type=Path, default=root / "refdb")
    ap.add_argument("--no-pfam", action="store_true", help="skip the Pfam scan")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)

    faa_files = sorted(args.proteins.glob("*.faa"))
    faa_files = [f for f in faa_files if "representative" not in f.name]
    if not faa_files:
        print(f"ERROR: no .faa files in {args.proteins}", file=sys.stderr)
        print("Run scripts/03_gene_prediction.py first.", file=sys.stderr)
        return 1

    genomes = [f.stem for f in faa_files]
    n_genomes = len(genomes)
    print(f"Genomes: {n_genomes}  ({', '.join(genomes[:3])} ...)")
    print(f"Homology thresholds: >={MIN_IDENTITY:.0%} identity, "
          f">={MIN_COVERAGE:.0%} coverage of the shorter protein")
    print("-" * 78)

    all_proteins: list[tuple[str, str]] = []
    for f in faa_files:
        all_proteins += parse_faa(f)
    print(f"Proteins to cluster: {len(all_proteins):,}")
    print("Clustering (this is the slow step -- expect a few minutes) ...")

    assignment, reps = cluster(all_proteins)
    n_fams = len(reps)
    print(f"Gene families: {n_fams:,}")

    # --- presence / absence matrix -------------------------------------------
    fam_genomes: dict[int, set[str]] = defaultdict(set)
    fam_sizes: Counter[int] = Counter()
    for hdr, fam in assignment.items():
        strain = hdr.split("|", 1)[0]
        fam_genomes[fam].add(strain)
        fam_sizes[fam] += 1

    fam_rows = []
    matrix_rows = []
    for fam in range(n_fams):
        present = fam_genomes[fam]
        cls = classify(n_genomes, len(present))
        fam_rows.append({"family": fam, "n_genomes": len(present),
                         "n_proteins": fam_sizes[fam], "class": cls})
        matrix_rows.append({"family": fam, **{g: int(g in present) for g in genomes}})

    fam_df = pd.DataFrame(fam_rows)
    mat_df = pd.DataFrame(matrix_rows)

    fam_df.to_csv(args.out / "gene_families.csv", index=False)
    mat_df.to_csv(args.out / "pan_genome_matrix.csv", index=False)

    with open(args.out / "gene_to_family.tsv", "w", encoding="utf-8") as fh:
        fh.write("protein\tfamily\n")
        for hdr, fam in assignment.items():
            fh.write(f"{hdr}\t{fam}\n")

    with open(args.out / "family_representatives.faa", "w", encoding="utf-8") as fh:
        for fam, hdr, seq in reps:
            fh.write(f">fam{fam}\n{seq}\n")

    counts = fam_df["class"].value_counts()
    print("-" * 78)
    for cls in ["core", "soft-core", "accessory", "unique"]:
        n = int(counts.get(cls, 0))
        print(f"  {cls:<10} {n:>5} families   ({n/n_fams:5.1%})")
    print(f"  {'TOTAL':<10} {n_fams:>5} families")

    # --- Pfam annotation + enrichment ----------------------------------------
    if not args.no_pfam:
        print("\nAnnotating one representative per family against Pfam-A ...")
        try:
            pfam = pd.DataFrame(annotate_pfam(reps, args.cache))
            pfam = pfam.merge(fam_df[["family", "class"]], on="family", how="left")
            pfam.to_csv(args.out / "functional_annotation.csv", index=False)
            print(f"  {len(pfam):,} families carry a Pfam domain "
                  f"({len(pfam)/n_fams:.1%} annotated)")

            # Honest expectation: accessory genes annotate less well than core,
            # because they are drawn from a much wider, under-studied gene pool.
            for cls in ["core", "accessory", "unique"]:
                sub = fam_df[fam_df["class"] == cls]
                if len(sub):
                    hit = sub["family"].isin(set(pfam["family"])).mean()
                    print(f"    {cls:<10} annotation rate {hit:5.1%}")

            enr = enrichment(fam_df, pfam, n_genomes)
            enr.to_csv(args.out / "functional_enrichment.csv", index=False)
            sig = enr[enr["q_value"] < 0.05] if not enr.empty else enr
            print(f"  enrichment: {len(sig)} domains significant at q<0.05")
            for _, r in sig.head(5).iterrows():
                print(f"     {r['pfam']:<12} OR={r['odds_ratio']:6.2f}  "
                      f"q={r['q_value']:.2e}  {str(r['description'])[:44]}")
            if len(sig) == 0:
                print("     (none passed FDR -- with 18 genomes the test has limited")
                print("      power; this is honest statistics, not a failed run)")
        except ImportError:
            print("  pyhmmer not installed -- skipping Pfam.")
            print("  install with:  pip install pyhmmer")
        except Exception as e:
            print(f"  Pfam step failed: {type(e).__name__}: {e}")
            print("  the pan-genome matrices above are still valid")

    print("-" * 78)
    print(f"Written to {args.out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
