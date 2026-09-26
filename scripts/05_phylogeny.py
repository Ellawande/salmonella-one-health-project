#!/usr/bin/env python3
"""
05_phylogeny.py  --  Stage 5: core-genome phylogeny.

Builds a tree from the SINGLE-COPY CORE genome: the gene families present
exactly once in every genome. Each family is aligned individually with FAMSA,
the alignments are concatenated into one supermatrix, and a Neighbour-Joining
tree is built on pairwise p-distances.

WHY SINGLE-COPY, NOT JUST "CORE"
    Two different sets are easy to confuse:

      core families           present in ALL genomes
      single-copy core        present EXACTLY ONCE in all genomes

    A family present twice in one genome means paralogs. Concatenating a
    paralogous family aligns gene A from one genome against its paralog B in
    another, which invents differences that have nothing to do with the
    organism's history. Single-copy filtering removes that artefact, which is
    why the tree uses the smaller set. This script reports both counts so the
    distinction is explicit rather than implied.

METHOD AND HONEST SCOPE
    Alignment : FAMSA (via pyfamsa) -- a real multiple sequence alignment
    Distance  : pairwise p-distance (proportion of differing columns)
    Tree      : Neighbour-Joining
    This is a distance-based NJ tree. It is a real, defensible tree, but it is
    NOT a maximum-likelihood tree and it carries NO bootstrap support. For a
    publication-grade tree with ModelFinder model selection and ultrafast
    bootstrap, take the concatenated alignment written here to IQ-TREE or
    RAxML (e.g. on Galaxy). Say so plainly if this tree is used in a paper.

USAGE
    python scripts/05_phylogeny.py
    python scripts/05_phylogeny.py --out results --proteins proteins

OUTPUT
    results/core_supermatrix.faa   the concatenated alignment (feed this to IQ-TREE)
    results/core_genome.nwk        the Neighbour-Joining tree (Newick)
    results/core_gene_stats.csv    how many families were used

DEPENDENCIES
    pip install pyfamsa scikit-bio pandas
"""
from __future__ import annotations

from pathlib import Path
import argparse
import sys

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def parse_faa(path: Path) -> list[tuple[str, str]]:
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


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Core-genome NJ phylogeny.")
    ap.add_argument("--proteins", type=Path, default=root / "proteins")
    ap.add_argument("--results", type=Path, default=root / "results")
    args = ap.parse_args()

    out = args.results
    out.mkdir(parents=True, exist_ok=True)

    for needed in ["pan_genome_matrix.csv", "gene_families.csv", "gene_to_family.tsv"]:
        if not (out / needed).exists():
            print(f"ERROR: {out/needed} not found.", file=sys.stderr)
            print("Run scripts/04_pangenome_annotation.py first.", file=sys.stderr)
            return 1

    matrix = pd.read_csv(out / "pan_genome_matrix.csv")
    fam_df = pd.read_csv(out / "gene_families.csv")
    genomes = [c for c in matrix.columns if c != "family"]
    n_genomes = len(genomes)
    print(f"Genomes: {n_genomes}")

    # ---- the two counts, kept explicitly distinct ---------------------------
    core = fam_df[fam_df["n_genomes"] == n_genomes]
    single_copy = core[core["n_proteins"] == n_genomes]
    print(f"  core families        (present in all {n_genomes}): {len(core):,}")
    print(f"  single-copy core     (exactly once in all {n_genomes}): "
          f"{len(single_copy):,}")
    print(f"  -> the tree uses the single-copy set "
          f"({len(core)-len(single_copy):,} paralogous families excluded)")

    # ---- load every protein, index by (strain, family) -----------------------
    fam_of: dict[str, int] = {}
    for line in (out / "gene_to_family.tsv").read_text().splitlines()[1:]:
        hdr, fam = line.split("\t")
        fam_of[hdr] = int(fam)

    seq_of: dict[tuple[str, int], str] = {}
    for faa in sorted(args.proteins.glob("*.faa")):
        if "representative" in faa.name:
            continue
        for hdr, seq in parse_faa(faa):
            fam = fam_of.get(hdr)
            if fam is not None:
                seq_of[(hdr.split("|", 1)[0], fam)] = seq
    print(f"  indexed {len(seq_of):,} protein sequences")

    # ---- align each single-copy core family, then concatenate ---------------
    from pyfamsa import Aligner, Sequence
    aligner = Aligner()

    single_copy_ids = single_copy["family"].tolist()
    blocks: list[dict[str, str]] = []
    skipped = 0

    print(f"\nAligning {len(single_copy_ids):,} single-copy core families with FAMSA ...")
    for i, fam in enumerate(single_copy_ids, 1):
        seqs = []
        for g in genomes:
            s = seq_of.get((g, fam))
            if s is None:
                break
            seqs.append(Sequence(g.encode(), s.encode()))
        if len(seqs) != n_genomes:
            skipped += 1
            continue

        msa = aligner.align(seqs)
        block = {aln.id.decode(): aln.sequence.decode() for aln in msa}
        blocks.append(block)

        if i % 500 == 0:
            print(f"    ... {i:,}/{len(single_copy_ids):,}", flush=True)

    print(f"  aligned families used : {len(blocks):,}")
    if skipped:
        print(f"  skipped (missing copy): {skipped:,}")

    # ---- concatenate into the supermatrix -----------------------------------
    concat = {g: [] for g in genomes}
    total_columns = 0
    for block in blocks:
        width = max(len(v) for v in block.values())
        for g in genomes:
            # pad short (gapped) alignments so every genome has equal width
            concat[g].append(block[g].ljust(width, "-"))
        total_columns += width

    supermatrix = {g: "".join(concat[g]) for g in genomes}
    print(f"\nSupermatrix: {total_columns:,} columns x {n_genomes} taxa")

    with open(out / "core_supermatrix.faa", "w", encoding="utf-8") as fh:
        for g in genomes:
            fh.write(f">{g}\n")
            s = supermatrix[g]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ---- pairwise p-distance ------------------------------------------------
    import numpy as np
    dist = np.zeros((n_genomes, n_genomes))
    for i in range(n_genomes):
        for j in range(i + 1, n_genomes):
            a, b = supermatrix[genomes[i]], supermatrix[genomes[j]]
            informative = [(x, y) for x, y in zip(a, b)
                           if x != "-" and y != "-"]
            if not informative:
                d = 1.0
            else:
                d = sum(1 for x, y in informative if x != y) / len(informative)
            dist[i, j] = dist[j, i] = d

    print(f"p-distance range: {dist[dist > 0].min():.4f} - {dist.max():.4f}")

    # ---- Neighbour-Joining tree --------------------------------------------
    from skbio import DistanceMatrix
    from skbio.tree import nj

    dm = DistanceMatrix(dist, ids=genomes)
    tree = nj(dm)

    nwk = out / "core_genome.nwk"
    nwk.write_text(str(tree) + "\n", encoding="utf-8")

    pd.DataFrame([{
        "n_genomes": n_genomes,
        "core_families": len(core),
        "single_copy_core_families": len(single_copy),
        "families_used_in_tree": len(blocks),
        "supermatrix_columns": total_columns,
        "min_p_distance": round(float(dist[dist > 0].min()), 5),
        "max_p_distance": round(float(dist.max()), 5),
    }]).to_csv(out / "core_gene_stats.csv", index=False)

    print("-" * 78)
    print(f"Tree written : {nwk}")
    print(f"Supermatrix  : {out/'core_supermatrix.faa'}")
    print("\nREMINDER: this is a distance-based NJ tree without bootstrap support.")
    print("For a paper, load core_supermatrix.faa into IQ-TREE or RAxML")
    print("(ModelFinder + 1000 ultrafast bootstrap) and report those support values.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
