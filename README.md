# Salmonella One Health

**Genomic epidemiology of antimicrobial resistance across the human–livestock interface**

`18 genomes` · `3 host groups` · `50 resistance genes` · `8,829 pan-genome families`

---

## Summary

Antimicrobial resistance (AMR) in *Salmonella enterica* is a textbook One Health problem:
the same resistance genes turn up in patients, in poultry, and in cattle. But *finding the same
gene in two hosts is not evidence that it spread between them.* Two very different mechanisms
produce the identical observation:

1. **Vertical inheritance** — the isolates share a recent common ancestor that already carried the gene.
2. **Horizontal gene transfer (HGT)** — the gene was acquired independently, in different lineages, in different hosts.

This project distinguishes the two by testing each resistance gene against a **core-genome
phylogeny**, using Fitch parsimony with a permutation null. The result is that resistance in
these isolates is split roughly 2:1 between the two mechanisms — and the distinction changes
the public-health interpretation completely.

---

## Research question

> Is antimicrobial resistance shared between human and livestock *Salmonella* because of
> **shared ancestry** (clonal expansion), or because of **horizontal gene transfer** across
> the animal–human interface?

Answering this requires a phylogeny, not a presence/absence table. Ancestry is the null
hypothesis; HGT is the deviation from it.

---

## Key findings

**1. All 18 isolates are the same species and extremely closely related.**
Pairwise ANI ranges from **98.48% to 99.99%** — every pair is far above the ~95% species
boundary. *Salmonella enterica* is highly clonal, so this is expected — but it means that
"shared resistance" is the *default* expectation and must be tested, not assumed.

**2. Genomic relatedness is not structured by host.**
Mean ANI between isolates from the **same host is 98.86%**, versus **98.92% between
different hosts**. Isolates from different host species are, on average, very slightly *more*
similar to each other than isolates from the same host. Host is not the organising principle
of genomic variation in this set.

**3. Both mechanisms of resistance sharing are operating simultaneously.**

| Gene | Isolates | Host groups | Gain events needed | p (clumped) | Interpretation |
| --- | --- | --- | --- | --- | --- |
| `AAC(6')-Iy` | 16 / 18 | 3 | **1** | 0.027 | **Inherited** — one ancient acquisition, spread with the lineage |
| `FosA7` | 4 | 3 | 2 | 0.027 | **Inherited** |
| `AAC(6')-Iaa` | 2 | 2 | 1 | 0.026 | **Inherited** |
| `tet(A)` | 3 | **3** | **3** | 1.000 | **Horizontally acquired** — three independent gains |
| `sul1` | 3 | 2 | 3 | 1.000 | **Horizontally acquired** |
| `APH(3'')-Ib` | 2 | 2 | 2 | 1.000 | **Horizontally acquired** |
| `APH(6)-Id` | 2 | 2 | 2 | 1.000 | **Horizontally acquired** |
| `aadA`, `aadA2`, `floR` | 2 each | 1–2 | 2 each | 1.000 | **Horizontally acquired** |

**The two headline cases:**

- **`AAC(6')-Iy` is present in 16 of 18 isolates across all three host groups but requires only
  a single evolutionary gain.** It sits on the tree exactly where ancestry places it. This is
  *not* cross-host transmission — it is one ancient acquisition that propagated with the
  dominant lineage. A naive presence/absence analysis would have reported this as
  "resistance shared between humans and livestock."

- **`tet(A)` is present in a human, a poultry, and a cattle isolate, and requires three
  independent gains** on the tree. The same gene was acquired separately, three times, in
  three different host-associated lineages. This is the genuine signature of horizontal gene
  transfer — and it is the defensible genomic evidence that tetracycline resistance moves
  across the animal–human interface.

**4. The accessory genome is large and the resistance complement is mostly mobilisable.**
Of 50 detected resistance genes, **20 are acquired/mobile** (plasmid- or transposon-associated)
rather than intrinsic chromosomal. **12 of 18 isolates (67%) carry at least one plasmid
replicon**, with `IncFII(S)` the most common (7 isolates) — the plasmid family most frequently
associated with AMR carriage in Enterobacteriaceae.

---

## Dataset

18 publicly available *Salmonella enterica* genomes, six per host group, spanning 12 serovars
and 11 geographic locations. All are real, downloaded assemblies with accessions recorded in
`data/metadata/isolates.csv`.

| Host group | n | Serovars |
| --- | --- | --- |
| Human | 6 | Enteritidis, Typhimurium, Saintpaul, Paratyphi B, Limete, (unresolved) |
| Poultry | 6 | Enteritidis, Gallinarum, Heidelberg, Infantis, Montevideo, Newport |
| Cattle | 6 | Enteritidis, Typhimurium, Dublin, Agona, Heidelberg, (unresolved) |

**Geography:** Belgium, Canada (British Columbia), China, China (Shanghai), France, Haiti,
Hungary, Mexico, USA, USA (SD), United Kingdom.
**Collection dates:** recorded for 13 of 18 isolates.

> **Honest scope:** 18 genomes is a *feasibility* set — enough to demonstrate the method and
> generate a hypothesis, not enough to support population-level inference. Scaling the cohort
> is the next phase (see *Limitations*).

---

## Pipeline

Six numbered, rerunnable steps in `scripts/`. Each is a thin wrapper around a real, tested
implementation — nothing is hand-waved.

| Script | Stage | Method |
| --- | --- | --- |
| `01_fetch_isolates.py` | Data acquisition | NCBI Pathogen Detection + Datasets API; real assembly accessions |
| `02_amr_and_plasmids.py` | Resistance & mobility | CARD (AMR), VFDB (virulence), PlasmidFinder (replicons) by sequence homology |
| `03_gene_prediction.py` | Gene calling | Prodigal via pyrodigal, trained per genome |
| `04_pangenome_annotation.py` | Pan-genome | k-mer prefilter + verified pairwise alignment (≥70% identity, ≥70% coverage); Pfam-A domains via pyhmmer |
| `05_phylogeny.py` | Core-genome tree | Single-copy core genes → FAMSA alignment → concatenation → Neighbour-Joining |
| `06_ani_and_amr_overlay.py` | Relatedness & signal test | sourmash MinHash ANI; Fitch parsimony + 1,000-permutation null per gene |

---

## Methods

**Relatedness (ANI).** sourmash MinHash sketches, k=31, scaled=1000. Reports relatedness
*between the supplied genomes* — it is not a taxonomic species call against reference type
strains (that requires GTDB-Tk or FastANI-vs-type-strains).

**Core-genome phylogeny.** 3,293 single-copy core gene families (present exactly once in all
18 genomes) were aligned individually with FAMSA, concatenated into a **1,031,143 amino-acid
column** supermatrix, and used to build a Neighbour-Joining tree on pairwise p-distances.

**Phylogenetic signal of AMR.** For each resistance gene, Fitch parsimony computes the minimum
number of independent gain/loss events required to explain its distribution on the tree. This
observed count is compared against a null distribution from **1,000 random permutations** of
the same number of positives across the tips. Genes requiring no more changes than chance
(p ≥ 0.05, ≥2 changes) are flagged as repeatedly acquired; genes requiring significantly fewer
changes (p < 0.05) are flagged as inherited with the lineage.

---

## Results

Full tables are in `results/`; figures in `results/figures/`.

- `ani_matrix.csv` — pairwise ANI across all 18 genomes
- `core_genome.nwk` — the core-genome phylogeny (Newick; openable in iTOL/FigTree)
- `core_genome_distance.csv` — the underlying core-gene distance matrix
- `pan_genome_matrix.csv` — presence/absence matrix, 18 strains × 8,829 gene families
- `amr_phylogenetic_signal.csv` — per-gene parsimony test (observed vs expected changes, p)
- `amr_inherited_vs_acquired.csv` — genes classified by acquisition mechanism
- `plasmid_replicons_all.csv` — replicon types detected per isolate, with % identity/coverage
- `virulence_genes.csv` — 206 virulence-associated genes detected

**Pan-genome partition** (8,829 gene families total):

| Class | Gene families |
| --- | --- |
| Core (18/18) | 3,341 |
| Soft-core (16–17) | 453 |
| Accessory (2–15) | 2,239 |
| Unique (1) | 2,795 |

Note the deliberate distinction between the **3,341 core** families (present in all 18) and the
**3,293 single-copy core** families (present *exactly once* in all 18) used for the tree —
single-copy filtering removes paralogs and their alignment artefacts, which is why the tree
uses the smaller set.

---

## Limitations

Stated plainly, because they determine what this analysis can and cannot claim:

- **Sample size.** 18 genomes gives limited statistical power for phylogenetic signal on 18
  tips. The permutation p-values are real but modest; the findings are hypothesis-generating.
  Population-level claims require a scaled cohort.
- **Tree inference.** The phylogeny is distance-based Neighbour-Joining without bootstrap
  support. A maximum-likelihood tree with bootstrap (IQ-TREE/RAxML) and **recombination
  detection** is required for publication — recombination in *Salmonella* can mislead
  phylogenies.
- **Parsimony infers minimum events.** It can demonstrate that multiple independent gains are
  necessary, but it cannot establish the *direction* of a transfer or name a donor–recipient
  pair.
- **Gene–plasmid linkage is association, not proof.** The co-occurrence of acquired resistance
  genes and plasmid replicons in the same isolate is consistent with plasmid carriage, but
  proving a specific gene sits on a specific plasmid needs assembly + plasmid binning.
- **Metagenomic/gene-content caveat.** Resistance detection is sequence-homology based; it
  reports genes present in the assembly and does not measure expression or phenotypic
  resistance.

---

## Reproducing this analysis

```bash
# clone / initialise
git init
git add .
git commit -m "Salmonella One Health: AMR acquisition analysis"

# environment
pip install -r requirements.txt

# run the pipeline in order
python scripts/01_fetch_isolates.py
python scripts/02_amr_and_plasmids.py
python scripts/03_gene_prediction.py
python scripts/04_pangenome_annotation.py
python scripts/05_phylogeny.py
python scripts/06_ani_and_amr_overlay.py

# explore the results interactively
streamlit run dashboard/app.py
```

> Git operations are run by you in a terminal — this analysis was produced in a notebook
> environment and the repository is prepared here as a complete, ready-to-commit project.

---

## Next steps

1. Scale the cohort to hundreds of isolates per host for population-level inference.
2. Maximum-likelihood phylogeny with bootstrap support (IQ-TREE) + recombination detection.
3. Plasmid binning to establish gene–plasmid linkage rather than co-occurrence.
4. Statistical association testing corrected for population structure (clonality).

---

## License

MIT — see `LICENSE`.

*Analyses performed with Omicsboard Lab.*
