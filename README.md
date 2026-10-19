# MagCluster <img src='docs/_static/logo_new.svg' align='right' width=220 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.0-orange?style=flat-square)](https://github.com/RunJiaJi/magcluster)

[English](README.md) | [中文说明](README_zh.md)

**MagCluster** is a tool for finding, annotating, and visualizing magnetosome gene clusters (MGCs) in bacterial genomes.

---

## What's New in v1.0

- **One command**: `magcluster run` runs everything from gene prediction to cluster visualization in a single step.
- **No Prokka needed**: gene calling is handled natively by [Pyrodigal](https://github.com/althonos/pyrodigal) (fast, no Perl or external setup).
- **92 HMM profiles**: curated protein HMMs covering 8 magnetosome families (`Mam`, `Mms`, `Mad`, `Man`, `Mae`, `Maq`, `Mca`, `Feo`), searched using [PyHMMER](https://github.com/althonos/pyhmmer).
- **Flexible inputs**: accepts GenBank (`.gbk`/`.gbff`), nucleotide FASTA contigs (`.fasta`/`.fna`/`.fa`), and protein FASTA (`.faa`), including `.gz` compressed files.
- **Per-genome outputs**: writes separate `{genome}_mgc_genes.tsv` and `{genome}_magnetosome_proteins.faa` files for each genome.


---

## Installation

### 1. Create a Conda environment (recommended)

```bash
conda create -n magcluster-1.0 -c conda-forge -c bioconda python=3.12 pip
conda activate magcluster-1.0
```

### 2. Install MagCluster

From GitHub:
```bash
pip install git+https://github.com/RunJiaJi/magcluster.git
```

Or from local source:
```bash
git clone https://github.com/RunJiaJi/magcluster.git
cd MagCluster
pip install -e .
```

Verify the installation:
```bash
magcluster --version
magcluster run --help
```

---

## Quick Start

### 1. Annotated GenBank file (`.gbff` / `.gbk`)
```bash
magcluster run genome.gbff -o results/
```

### 2. Unannotated FASTA contigs (`.fasta` / `.fna` / `.fa`)
MagCluster automatically predicts genes and uses all CPU cores by default:
```bash
magcluster run contigs.fasta -o results/ -t 8
```

### 3. Multiple genomes in a folder
```bash
magcluster run genomes_dir/ -o results/
```

### 4. Protein FASTA (`.faa`)
Note if your input contains protein sequences without coordinates, MagCluster will annotate magnetosome proteins and output TSV table and FASTA files without plotting the gene cluster:
```bash
magcluster run proteins.faa -o results/
```

### 5. Overwrite an existing output folder
```bash
magcluster run genome.gbff -o results/ --force
```

---

## Commands & Options

### `magcluster run`

```bash
magcluster run [OPTIONS] INPUTS...
```

| Option | Shorthand | Default | Description |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | `magcluster_results` | Output directory |
| `--force` | `-f` | `False` | Overwrite files if the output directory is not empty |
| `--threads` | `-t` | `0` | Number of CPU threads (`0` = use all cores) |
| `--evalue` | `-e` | `1e-5` | HMM search E-value cutoff |
| `--score` | `-s` | `25.0` | Minimum HMM bit score |
| `--max-gap` | | `5` | Maximum non-magnetosome genes allowed between hits in a cluster |
| `--threshold` | `-th` | `3` | Minimum number of distinct magnetosome families required to call an MGC |
| `--flank` | `-un` | `3` | Number of non-magnetosome flanking genes to keep on each side |
| `--meta / --single` | | Auto | Force metagenome mode or single-genome mode for gene prediction |
| `--identity` | `-i` | `0.3` | Minimum sequence identity cutoff for Clinker alignments (0.0–1.0) |
| `--clinker-opts` | | `None` | Extra arguments passed to Clinker (e.g. `"-na -ufo"`) |
| `--clinker-help` | | | Show available Clinker options guide and exit |
| `--static-plot` | | `False` | Also generate static SVG figures |
| `--no-plot` | | `False` | Skip Clinker visualization (tables only) |

### `magcluster db`

View the built-in reference database, the 10 core defining families, and all 92 HMM profiles sorted by sequence count:

```bash
magcluster db
```

---

## Output Files

Running `magcluster run <inputs> -o results/` produces:

```text
results/
├── {genome}_mgc_genes.tsv              # Table of all identified magnetosome genes
├── {genome}_magnetosome_proteins.faa   # Extracted magnetosome protein sequences
├── MGC_align.html                      # Interactive gene cluster map (open in browser)
└── clusters_gbk/                       # Sliced GenBank file for each cluster
    └── {cluster_id}.gbk
```

### 1. Gene table (`{genome}_mgc_genes.tsv`)

| Column | Meaning |
| :--- | :--- |
| `cluster_id` | Cluster ID (e.g. `MSR-1_genome_MGC`), or `unclustered` for isolated genes |
| `genome` | Genome name |
| `gene_id` | Locus tag or gene identifier |
| `contig` | Contig or scaffold name |
| `start`, `end` | Gene coordinates (1-based) |
| `strand` | Strand (`+` or `-`) |
| `profile` | Best matching HMM profile (e.g. `MamE_full`, `MamE_Cter`, `MamQ_greigite`) |
| `parent_family` | Family name (e.g. `MamE`, `MamQ`) |
| `role` | `core`, `accessory`, or `iron_context` |
| `score`, `evalue` | HMM match score and E-value |

### 2. Protein sequences (`{genome}_magnetosome_proteins.faa`)

Amino acid sequences of all magnetosome proteins found in that genome, with detailed header info:
```text
>MGMSRV2_RS11765 [MSR-1_genome] cluster=MSR-1_genome_MGC family=MamE profile=MamE_full score=1000.4 evalue=1.97e-301
MSERRPVGRTL...
```

### 3. Interactive visualization (`MGC_align.html`)

Open this file in your webbrowser:

- **Gene labels**: Arrows show specific domain annotations like `MamE_Cter`, `MamE_Nter`, and `MamQ_greigite`.
- **Family colors**:
  - Red: `mam` family
  - Yellow: `mms` family
  - Green: `mad` family
  - Blue: `man` family
  - Gold: `feo` family
  - Purple: `mae` family
  - Teal: `maq` family
  - Pink: `mca` family
  - Grey: Other neighboring non-magnetosome genes
- **Customizable**: Drag genes, change colors, re-align, and export publication-ready **SVG** figures directly from the browser.

---

## Helpful Tips

### Faster Clinker plotting on large datasets
Pairwise sequence alignments between large clusters can take time. Use `--clinker-opts` to speed it up:

```bash
# Skip pairwise sequence alignment:
magcluster run *.gbff -o results/ --clinker-opts "-na"

# Keep the original file order instead of automatic reordering:
magcluster run *.gbff -o results/ --clinker-opts "-ufo"

# Both:
magcluster run *.gbff -o results/ --clinker-opts "-na -ufo"
```

To see all Clinker visualization options, run:
```bash
magcluster run --clinker-help
```

---

## Core Magnetosome Families

MagCluster considers an MGC authentic if it contains at least 3 distinct magnetosome families. The 10 core defining families are:

$$\text{MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ}$$

---

## Citation

If you use MagCluster, please cite:

> Runjia Ji, Wensi Zhang, Yongxin Pan, and Wei Lin. (2022).  
> **MagCluster: a Tool for Identification, Annotation, and Visualization of Magnetosome Gene Clusters.**  
> *Microbiology Resource Announcements*, 11:e01031-21.  
> [https://doi.org/10.1128/mra.01031-21](https://doi.org/10.1128/mra.01031-21)

Dependencies used by MagCluster:
- **Pyrodigal**: [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **PyHMMER**: [doi:10.1093/bioinformatics/btad211](https://doi.org/10.1093/bioinformatics/btad211)
- **Clinker**: [doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007)

---

## Contact

Questions or suggestions? Feel free to reach out:
- **GitHub Issues**: [Open an issue](https://github.com/RunJiaJi/magcluster/issues)
- **Email**: Runjia Ji (`jirunjia@gmail.com`) or Wei Lin (`weilin@mail.iggcas.ac.cn`)
