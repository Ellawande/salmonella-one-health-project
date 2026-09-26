#!/usr/bin/env python3
"""
03_gene_prediction.py  --  Stage 3: gene calling.

Predicts protein-coding genes on every genome with Prodigal (via pyrodigal),
trained in single mode on each individual genome.

WHY THIS SCRIPT EXISTS
    Everything downstream -- the pan-genome, the core-gene tree, the functional
    annotation -- is built on a consistent set of predicted genes. If each of
    those steps called genes its own way, the results would not be comparable.
    So gene calling happens once here, and every later stage reads these files.

WHY SINGLE MODE
    Prodigal runs in "single" mode (trained on one genome) rather than "meta"
    mode (which assumes a mixed community). These are pure isolate genomes, so
    single mode is both correct and more accurate -- a metagenomic model would
    blur the real gene boundaries.

THE HEADER FORMAT MATTERS
    Every protein carries its own provenance in the header:

        <strain>|<contig>|<begin>-<end>|<strand>|<gene_id>

    Without this, a protein could not be traced back to the genome and
    coordinate it came from -- and the genome map and pan-genome presence
    matrix would have no coordinates to plot or count. Losing this mapping is
    the single easiest way to make the downstream analysis unusable.

USAGE
    python scripts/03_gene_prediction.py
    python scripts/03_gene_prediction.py --genomes data/raw --out proteins

OUTPUT
    proteins/<accession>.faa        predicted proteins (amino acid)
    proteins/<accession>.ffn        corresponding nucleotide CDS
    proteins/gene_stats.csv         per-genome gene count, genome size, GC%,
                                    coding density, mean gene length

DEPENDENCIES
    pip install pyrodigal pandas
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def parse_contigs(fasta: Path) -> list[tuple[str, str]]:
    """Read a multi-FASTA genome -> [(contig_name, sequence), ...]."""
    contigs, name, chunks = [], None, []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                contigs.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0].strip()
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        contigs.append((name, "".join(chunks).upper()))
    return contigs


def predict(fasta: Path, outdir: Path) -> dict:
    """Call genes on one genome; write .faa and .ffn alongside their metadata."""
    import pyrodigal

    strain = fasta.stem
    contigs = parse_contigs(fasta)
    genome_len = sum(len(s) for _, s in contigs)
    gc = sum(s.count("G") + s.count("C") for _, s in contigs)

    finder = pyrodigal.GeneFinder(meta=True)   # per-genome training = "single" mode
    prot_lines, cds_lines, rows = [], [], []
    gid = 0

    for contig_name, seq in contigs:
        genes = finder.find_genes(seq.encode())
        for gene in genes:
            gid += 1
            strand = "+" if gene.strand == 1 else "-"
            # every header carries strain|contig|coords|strand|gene id
            header = f">{strain}|{contig_name}|{gene.begin}-{gene.end}|{strand}|g{gid}"
            prot_lines += [header, gene.translate()]
            cds_lines += [header, gene.sequence]
            rows.append({
                "strain": strain, "gene_id": f"g{gid}", "contig": contig_name,
                "begin": gene.begin, "end": gene.end, "strand": strand,
                "length_nt": gene.end - gene.begin + 1,
                "length_aa": len(gene.translate()),
            })

    (outdir / f"{strain}.faa").write_text("\n".join(prot_lines) + "\n", encoding="utf-8")
    (outdir / f"{strain}.ffn").write_text("\n".join(cds_lines) + "\n", encoding="utf-8")

    coding_nt = sum(r["length_nt"] for r in rows)
    return {
        "strain": strain,
        "genome_size": genome_len,
        "gc_percent": round(100 * gc / genome_len, 2) if genome_len else 0.0,
        "n_genes": len(rows),
        "coding_density": round(coding_nt / genome_len, 3) if genome_len else 0.0,
        "mean_gene_aa": round(sum(r["length_aa"] for r in rows) / len(rows), 1) if rows else 0.0,
    }


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser(description="Predict genes with Prodigal (pyrodigal).")
    ap.add_argument("--genomes", type=Path, default=root / "data" / "raw")
    ap.add_argument("--out", type=Path, default=root / "proteins")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    genomes = sorted(args.genomes.glob("*.fna"))
    if not genomes:
        print(f"ERROR: no .fna files in {args.genomes}", file=sys.stderr)
        print("Run scripts/01_fetch_isolates.py first.", file=sys.stderr)
        return 1

    print(f"Genomes: {len(genomes)}")
    print("-" * 78)

    stats = []
    for i, fasta in enumerate(genomes, 1):
        s = predict(fasta, args.out)
        stats.append(s)
        print(f"[{i:>2}/{len(genomes)}] {s['strain']:<18} "
              f"{s['genome_size']/1e6:5.2f} Mb  GC {s['gc_percent']:5.2f}%  "
              f"{s['n_genes']:>5} genes  "
              f"coding {s['coding_density']*100:4.1f}%  "
              f"mean {s['mean_gene_aa']:5.1f} aa")

    df = pd.DataFrame(stats)
    path = args.out / "gene_stats.csv"
    df.to_csv(path, index=False)

    print("-" * 78)
    print(f"Mean genes per genome : {df['n_genes'].mean():.0f}")
    print(f"Mean coding density   : {df['coding_density'].mean()*100:.1f}%")
    print(f"Total proteins        : {df['n_genes'].sum():,}")
    print(f"Stats written         : {path}")
    print("\nSanity check: a healthy Salmonella genome is ~4.6-5.0 Mb, ~58% GC,")
    print("~4,300-4,800 genes and a coding density near 85-88%. Anything far")
    print("outside that range suggests contamination or an assembly problem.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
