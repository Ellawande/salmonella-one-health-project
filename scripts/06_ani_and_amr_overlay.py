#!/usr/bin/env python3
"""
06_ani_and_amr_overlay.py  --  Stage 6: relatedness and the signal test.

This is the scientific centrepiece. It asks, for every acquired AMR gene:

    Is this gene's distribution across human, poultry and cattle isolates
    explained by shared ancestry (inherited), or does it cut across lineages
    (consistent with horizontal transfer)?

THE IDEA
    If resistance were only ever inherited, two isolates with the same gene
    would tend to be close relatives. So the gene's presence/absence pattern
    across the tips of a phylogeny should require FEW changes of state -- the
    gene appears once and is passed down, a "clumped" pattern.

    If a gene moves horizontally, it will show up in distantly related
    isolates, requiring MANY changes of state -- a "dispersed" pattern.

    To decide whether "many" or "few" is meaningful we need a null. Shuffling
    the presence/absence labels across the tree tips 1,000 times gives the
    distribution of parsimony scores expected if the gene were placed at
    random with respect to phylogeny. Comparing the observed score to that
    null yields p_clumped.

METHOD
    Relatedness : sourmash MinHash signatures, k=31, scaled=1000, compared
                  pairwise; Jaccard -> ANI via the Mash containment formula
    Signal test : Fitch parsimony (changes of state) on the core-genome tree,
                  against a 1,000-permutation null
    p_clumped   : fraction of permutations whose parsimony score is <= the
                  observed score -- i.e. how often random placement looks at
                  least as clustered as the real data. Small p_clumped means
                  the real gene is MORE clumped than chance, which is the
                  signature of inheritance.

HONEST SCOPE -- READ THIS BEFORE QUOTING A RESULT
    * Clumped is evidence CONSISTENT WITH inheritance. It is not proof. A gene
      can also be clumped because it entered one lineage once and then spread
      within it.
    * Dispersed CLAIMS NOTHING about mechanism. It means "not explained by the
      tree" -- the cause could be horizontal transfer, but equally could be
      convergent loss, poor assembly, or an imperfect phylogeny.
    * The biggest confounder is clonality: if the 18 isolates contain several
      near-identical clones, those clones share everything, and shared genes
      will look inherited whether or not transfer occurred. Proper correction
      requires population-structure-aware models (see the README limitations).
    * Correlation between a gene and a plasmid replicon in the same isolate is
      CO-OCCURRENCE, not linkage. Proving a gene sits ON a plasmid needs
      assembly and plasmid binning.

USAGE
    python scripts/06_ani_and_amr_overlay.py
    python scripts/06_ani_and_amr_overlay.py --perms 1000

OUTPUT
    results/ani_matrix.csv                   pairwise ANI (%)
    results/amr_sharing_matrix.csv           gene x host-group counts
    results/amr_phylogenetic_signal.csv      the Fitch/permutation test
    results/amr_inherited_vs_acquired.csv    the test + a clumped/dispersed call
    results/onehealth_summary.csv            per-isolate one-line summary

DEPENDENCIES
    pip install sourmash scikit-bio pandas numpy
"""
from __future__ import annotations

import argparse
import io
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ANI_K = 31
ANI_SCALED = 1000
DEFAULT_PERMS = 1000
CLUMPED_ALPHA = 0.05


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


# ----------------------------------------------------------------------------
# ANI from sourmash MinHash sketches
# ----------------------------------------------------------------------------
def ani_matrix(genomes: list[Path]) -> pd.DataFrame:
    """
    Pairwise ANI. Jaccard similarity J is converted with the Mash formula:
        ANI = 1 + (1/k) * ln(2J / (1 + J))
    """
    import sourmash

    sketches = []
    for g in genomes:
        # sourmash wants DNA; read the FASTA into memory (a few MB per genome)
        seqs = "".join(line.strip() for line in g.read_text().splitlines()
                       if not line.startswith(">"))
        sig = sourmash.MinHash(n=0, ksize=ANI_K, scaled=ANI_SCALED)
        sig.add_sequence(seqs, force=True)
        sketches.append(sig)

    names = [g.stem for g in genomes]
    n = len(names)
    m = np.zeros((n, n))

    for i in range(n):
        m[i, i] = 100.0
        for j in range(i + 1, n):
            jaccard = sketches[i].jaccard(sketches[j])
            if jaccard <= 0:
                ani = 0.0
            else:
                ani = 100.0 * (1 + (1.0 / ANI_K) * np.log(2 * jaccard / (1 + jaccard)))
            m[i, j] = m[j, i] = ani

    return pd.DataFrame(m, index=names, columns=names)


# ----------------------------------------------------------------------------
# Fitch parsimony on a tree, for a binary character
# ----------------------------------------------------------------------------
def fitch_changes(tree, states: dict[str, int]) -> int:
    """
    Fitch's algorithm for a binary character.
    Returns the number of state changes required on the tree.
    """
    def visit(node):
        if node.is_tip():
            return {states.get(node.name, 0)}, 0
        child_results = [visit(c) for c in node.children]
        sets = [s for s, _ in child_results]
        changes = sum(c for _, c in child_results)

        intersection = set.intersection(*sets)
        if intersection:                       # shared state -> no extra change
            return intersection, changes
        return set.union(*sets), changes + 1   # conflict -> one change

    _, total = visit(tree.root())
    return total


def permutation_test(tree, tip_names: list[str], observed: dict[str, int],
                     n_perms: int, rng: np.random.Generator):
    """
    Returns (observed_changes, p_clumped, mean_null_changes).

    p_clumped = P(a random placement is at least as clumped as the real one).
    A SMALL p_clumped means the real gene is more clustered than chance, which
    is the signature of inheritance. The null mean is returned so the observed
    score can be reported against what randomness would give.
    """
    obs = fitch_changes(tree, observed)
    k = sum(observed.values())                  # how many tips carry the gene

    # no variation across tips -> the test is undefined, not "not significant"
    if k == 0 or k == len(tip_names):
        return obs, float("nan"), float("nan")

    states = np.zeros(len(tip_names), dtype=int)
    states[:k] = 1
    as_clumped_or_more = 0
    null_scores = np.empty(n_perms, dtype=int)

    for i in range(n_perms):
        rng.shuffle(states)                     # permute labels across tips
        score = fitch_changes(tree, dict(zip(tip_names, states.tolist())))
        null_scores[i] = score
        if score <= obs:
            as_clumped_or_more += 1

    p = (as_clumped_or_more + 1) / (n_perms + 1)   # add-one smoothing
    return obs, p, float(null_scores.mean())


def r2(x, nd: int = 2):
    """Round for CSV output; NaN becomes an empty cell rather than the string 'nan'."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    return round(float(x), nd)


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="ANI + AMR phylogenetic-signal test.")
    ap.add_argument("--genomes", type=Path, default=root / "data" / "raw")
    ap.add_argument("--results", type=Path, default=root / "results")
    ap.add_argument("--metadata", type=Path,
                    default=root / "data" / "metadata" / "isolates.csv")
    ap.add_argument("--perms", type=int, default=DEFAULT_PERMS)
    ap.add_argument("--seed", type=int, default=42, help="for reproducible permutations")
    args = ap.parse_args()

    out = args.results
    hits_path = out / "amr_hits.csv"
    tree_path = out / "core_genome.nwk"
    if not hits_path.exists():
        print(f"ERROR: {hits_path} not found -- run script 02 first.", file=sys.stderr)
        return 1
    if not tree_path.exists():
        print(f"ERROR: {tree_path} not found -- run script 05 first.", file=sys.stderr)
        return 1

    rng = np.random.default_rng(args.seed)

    # ---- 1. pairwise ANI ----------------------------------------------------
    genome_files = sorted(args.genomes.glob("*.fna"))
    print(f"Computing ANI across {len(genome_files)} genomes (sourmash, "
          f"k={ANI_K}, scaled={ANI_SCALED}) ...")
    ani = ani_matrix(genome_files)
    ani.to_csv(out / "ani_matrix.csv")
    off = ani.values[np.triu_indices(len(ani), 1)]
    print(f"  ANI range: {off.min():.2f}% - {off.max():.2f}%")
    print(f"  pairs above the 95% species boundary: {(off >= 95).sum()}"
          f"/{len(off)}")

    # ---- 2. load the tree ---------------------------------------------------
    from skbio import TreeNode
    try:
        tree = TreeNode.read(io.StringIO(tree_path.read_text()))
    except Exception:
        import skbio.io
        tree = skbio.io.read(io.StringIO(tree_path.read_text()), format="newick",
                             into=TreeNode)
    tip_names = [t.name for t in tree.tips()]
    print(f"\nTree loaded: {len(tip_names)} tips")

    # ---- 3. gene presence per isolate --------------------------------------
    hits = pd.read_csv(hits_path)
    acq = hits[(hits["database"] == "CARD") & (hits["acquired"] == True)]
    presence: dict[str, set[str]] = defaultdict(set)
    for _, r in acq.iterrows():
        presence[r["gene"]].add(r["strain"])

    # ---- 4. host groups ----------------------------------------------------
    meta = pd.read_csv(args.metadata)
    acc_col = next((c for c in ["asm_acc", "assembly_accession", "assembly",
                                "accession"] if c in meta.columns), None)
    host_col = next((c for c in ["host_group", "Host_Group", "hostgroup",
                                 "group"] if c in meta.columns), None)
    if acc_col is None or host_col is None:
        print("WARNING: could not identify accession/host columns; "
              "host breakdown will be omitted", file=sys.stderr)
        host_of = {}
    else:
        host_of = dict(zip(meta[acc_col].astype(str), meta[host_col].astype(str)))
    host_groups = sorted(set(host_of.values())) or ["Human", "Poultry", "Cattle"]

    # ---- 5. per-gene sharing + signal test ---------------------------------
    sharing_rows, signal_rows, pattern_rows = [], [], []
    print(f"\nTesting {len(presence)} acquired AMR genes "
          f"({args.perms:,} permutations each) ...")

    for gene in sorted(presence):
        carriers = presence[gene]
        # host breakdown
        counts = {h: 0 for h in host_groups}
        for strain in carriers:
            h = host_of.get(strain)
            if h in counts:
                counts[h] += 1
        n_hosts = sum(1 for h in host_groups if counts[h] > 0)

        sharing_rows.append({"gene": gene, **counts,
                             "n_hosts": n_hosts, "acquired": True})

        # Fitch + permutation null
        observed = {t: (1 if t in carriers else 0) for t in tip_names}
        changes_obs, p_clumped, expected = permutation_test(
            tree, tip_names, observed, args.perms, rng)

        signal_rows.append({"gene": gene, "acquired": True,
                            "n_isolates": len(carriers), "n_hosts": n_hosts,
                            "changes_observed": changes_obs,
                            "changes_expected": r2(expected),
                            "p_clumped": r2(p_clumped, 4)})

        if np.isnan(p_clumped):
            pattern = "invariant (no test)"
        elif p_clumped < CLUMPED_ALPHA:
            pattern = "clumped (inherited)"
        else:
            pattern = "dispersed (not explained by tree)"
        pattern_rows.append({"gene": gene, "acquired": True,
                             "n_isolates": len(carriers), "n_hosts": n_hosts,
                             "changes_observed": changes_obs,
                             "changes_expected": r2(expected),
                             "p_clumped": r2(p_clumped, 4),
                             "pattern": pattern})

    pd.DataFrame(sharing_rows).to_csv(out / "amr_sharing_matrix.csv", index=False)
    sig_df = pd.DataFrame(signal_rows)
    sig_df.to_csv(out / "amr_phylogenetic_signal.csv", index=False)
    pat_df = pd.DataFrame(pattern_rows).sort_values("p_clumped")
    pat_df.to_csv(out / "amr_inherited_vs_acquired.csv", index=False)

    # ---- 6. per-isolate one-health summary ---------------------------------
    repl_path = out / "plasmid_replicons.csv"
    repl_by_strain: dict[str, list[str]] = defaultdict(list)
    if repl_path.exists():
        rdf = pd.read_csv(repl_path)
        for _, r in rdf.iterrows():
            repl_by_strain[r["strain"]].append(r["replicon"])

    serovar_of = dict(zip(meta[acc_col].astype(str),
                          meta.get("serovar", pd.Series()).astype(str))) if acc_col else {}
    country_of = dict(zip(meta[acc_col].astype(str),
                          meta.get("geo_loc_name", meta.get("country", pd.Series())).astype(str))) if acc_col else {}

    summary = []
    for strain in sorted({s.split("|", 1)[0] for s in hits["strain"].unique()}):
        genes = sorted(g for g in presence if strain in presence[g])
        summary.append({
            "strain": strain,
            "host_group": host_of.get(strain, ""),
            "serovar": serovar_of.get(strain, ""),
            "country": country_of.get(strain, ""),
            "n_acquired_AMR": len(genes),
            "acquired_AMR": ", ".join(genes),
            "replicons": ", ".join(sorted(set(repl_by_strain.get(strain, [])))),
        })
    pd.DataFrame(summary).to_csv(out / "onehealth_summary.csv", index=False)

    # ---- report -------------------------------------------------------------
    print("-" * 78)
    if not pat_df.empty:
        print(f"Genes tested: {len(pat_df)}")
        for label, n in pat_df["pattern"].value_counts().items():
            print(f"  {label:<36} {n}")
        clumped = pat_df[pat_df["pattern"] == "clumped (inherited)"]
        if not clumped.empty:
            print("\nStrongest inheritance signal (lowest p_clumped):")
            for _, r in clumped.head(5).iterrows():
                print(f"  {r['gene']:<34} p={r['p_clumped']:.3f}  "
                      f"changes {r['changes_observed']} "
                      f"(null ~{r['changes_expected']})  hosts={r['n_hosts']}")

    print("\nFiles written:")
    for f in ["ani_matrix.csv", "amr_sharing_matrix.csv",
              "amr_phylogenetic_signal.csv", "amr_inherited_vs_acquired.csv",
              "onehealth_summary.csv"]:
        if (out / f).exists():
            print(f"  results/{f}")

    print("\nCAVEAT: 'clumped' is consistent with inheritance but is not proof;")
    print("        'dispersed' means only 'not explained by this tree'. Clonality")
    print("        among the isolates is the main confounder. See README.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
