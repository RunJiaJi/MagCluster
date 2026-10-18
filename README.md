# MagCluster <img src='docs/_static/logo_new.svg' align='right' width=240 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.0-orange?style=flat-square)](https://github.com/RunJiaJi/magcluster)
[![Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://bioconda.github.io/recipes/magcluster/README.html)

**MagCluster** is a high-performance, standalone bioinformatics tool for the **identification, annotation, classification, and interactive visualization** of Magnetosome Gene Clusters (MGCs) across magnetotactic bacteria (MTB) and environmental metagenomic datasets.

---

## 🌟 What's New in MagCluster v1.0

MagCluster v1.0 is a complete architectural redesign built for speed, accuracy, and ease of use:

- ⚡ **One-Step End-to-End Execution (`magcluster run`)**: Replaces the old multi-step workflow (`prokka` &rarr; `mgc_screen` &rarr; `clinker`) with a single, streamlined command.
- 🧬 **Prokka-Free Native Gene Calling**: Powered by [Pyrodigal](https://github.com/althonos/pyrodigal) (native C/Cython binding of Prodigal). No Perl dependencies, BLAST indexing, or external binaries required.
- 🔍 **Curated 92-Profile HMM Library**: Replaces heuristic keyword matching with **92 profile-HMMs** built from 817 representative magnetosome proteins across 8 major families (`Mam`, `Mms`, `Mad`, `Man`, `Mae`, `Maq`, `Mca`, `Feo`), searched via [PyHMMER](https://github.com/althonos/pyhmmer) (~20&times; faster than standard HMMER3).
- 📂 **Multi-Format Input Flexibility**: Supports pre-annotated GenBank (`.gbff`/`.gbk`), unannotated nucleotide contigs (`.fasta`/`.fna`), and protein FASTA (`.faa`), with transparent `.gz` decompression.
- 💎 **Synteny Clustering & Mineral Crystal Typing**: Groups magnetosome genes based on chromosome synteny and intergenic distances, automatically predicting the mineral type (**Magnetite** $\text{Fe}_3\text{O}_4$ vs. **Greigite** $\text{Fe}_3\text{S}_4$).
- 🎨 **Publication-Ready Interactive Visualization**: Generates an interactive web alignment (`MGC_align.html`) powered by Clinker with:
  - Specific domain and sub-profile labels (e.g., `MamE_Cter`, `MamE_Nter`, `MamQ_greigite`, `MamO_domain`).
  - Standard clean naming for full-length proteins without redundant suffixes (`MamE`, `MamO`, `Man5`).
  - Unified family color palettes (`mam` red, `mms` yellow, `mad` green, `man` blue, etc.).
  - Pass-through support for native Clinker options (`--clinker-opts`).
- 📁 **Per-Genome Clean Deliverables**: Generates `{genome}_mgc_genes.tsv` and `{genome}_magnetosome_proteins.faa` per genome without mixed files.
- 🛡️ **Overwrite Protection**: Safeguards existing non-empty output directories using `--force` (`-f`).

---

## 📋 Table of Contents

- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Input Formats](#-input-formats)
- [Commands & Options](#-commands--options)
  - [`magcluster run`](#magcluster-run)
  - [`magcluster db`](#magcluster-db)
- [Output Files & Interpretation](#-output-files--interpretation)
- [Clinker Customization Options](#-clinker-customization-options)
- [Citation](#-citation)
- [Contact & Support](#-contact--support)

---

## 📦 Installation

### Option 1: Conda Environment (Recommended)

Create a clean Conda environment with Python 3.10+:

```bash
# Create and activate environment
conda create -n magcluster -c conda-forge -c bioconda python=3.12 pip
conda activate magcluster

# Install MagCluster
pip install git+https://github.com/RunJiaJi/magcluster.git
```

Or install from source repository:
```bash
git clone https://github.com/RunJiaJi/magcluster.git
cd MagCluster
pip install -e .
```

### Verification

Check that MagCluster is installed correctly:
```bash
magcluster --version
magcluster --help
```

---

## 🚀 Quick Start

### 1. Analyze an Annotated GenBank / GBFF Genome

```bash
magcluster run MSR-1_genome.gbff -o results/
```

### 2. Analyze Unannotated FASTA Contigs (with Multi-threading)

MagCluster will automatically predict genes using Pyrodigal and analyze syntenic clusters across all CPU cores:

```bash
magcluster run contigs.fasta -o results/ -t 8
```

### 3. Batch Process Multiple Genomes in a Directory

```bash
magcluster run genomes_dir/ -o results/ --force
```

### 4. Analyze Protein FASTA Sequences (`.faa`)

When providing protein sequences without genomic coordinates, MagCluster performs HMM annotation and exports per-genome TSV and FASTA files while cleanly bypassing syntenic cluster mapping:

```bash
magcluster run proteins.faa -o results/
```

---

## 📄 Input Formats

MagCluster automatically detects input formats based on file extension and headers:

| File Type | Supported Extensions | Processing Pipeline |
| :--- | :--- | :--- |
| **GenBank / GBFF** | `.gbk`, `.gbff`, `.gbf`, `.gb` (or `.gz`) | Extracts pre-annotated CDS coordinates and translations directly. |
| **FASTA Contigs** | `.fasta`, `.fna`, `.fa`, `.fsa` (or `.gz`) | Predicts genes and translates proteins natively using Pyrodigal. |
| **Protein FASTA** | `.faa`, `.pep` (or `.gz`) | Annotates magnetosome proteins via HMM; exports tabular and sequence results. |

---

## 🛠 Commands & Options

### `magcluster run`

Run the end-to-end MagCluster analysis pipeline:

```bash
magcluster run [OPTIONS] INPUTS...
```

#### Available Options

| Option | Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | `magcluster_results` | Output directory path. |
| `--force` | `-f` | `False` | Overwrite existing files if the output directory already exists and is non-empty. |
| `--threads` | `-t` | `0` | Number of CPU threads (0 = utilize all available cores). |
| `--evalue` | `-e` | `1e-5` | HMM search E-value cutoff. |
| `--score` | `-s` | `25.0` | Minimum bit score threshold for HMM hits. |
| `--max-gap` | | `5` | Maximum intergenic non-MGC genes allowed between cluster hits. |
| `--threshold` | `-th` | `3` | Minimum number of distinct magnetosome families required to define an MGC. |
| `--flank` | `-un` | `3` | Number of flanking non-MGC genes to reserve on each side of the cluster. |
| `--meta / --single` | | Auto | Force metagenomic mode vs single-genome mode in Pyrodigal. |
| `--identity` | `-i` | `0.3` | Minimum sequence identity cutoff for Clinker alignments (0.0–1.0). |
| `--clinker-opts` | | `None` | Pass custom arguments directly to Clinker (e.g. `--clinker-opts "-na -ufo"`). |
| `--clinker-help` | | | Display available Clinker native options and usage guide, then exit. |
| `--static-plot` | | `False` | Also render static publication-quality SVG figures via pyGenomeViz. |
| `--no-plot` | | `False` | Skip Clinker alignment and visualization (faster headless execution). |

---

### `magcluster db`

Inspect the built-in reference database, version, defining core families, and all 92 HMM profiles sorted by reference sequence count:

```bash
magcluster db
```

*Example output:*
```text
=== MagCluster Reference Database Info ===
Rules version: 1.0.0
HMM database path: .../src/magcluster/data/hmms/magnetosome.hmm
Total HMM profiles: 92
Core defining families (10): MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ

         Reference Magnetosome HMM Profiles (92)          
┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┓
┃ Profile            ┃ Protein Family ┃ Ref Seqs ┃ Median length (aa) ┃
┡━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━┩
│ MamI               │ MamI           │       33 │                 76 │
│ MamK               │ MamK           │       31 │                347 │
│ MamQ               │ MamQ           │       31 │                270 │
│ MamA               │ MamA           │       30 │                217 │
│ MamB               │ MamB           │       27 │                296 │
│ MamE_full          │ MamE           │       27 │                728 │
│ MamM               │ MamM           │       25 │                307 │
│ FeoB               │ FeoB           │       24 │                699 │
│ MamP               │ MamP           │       24 │                305 │
│ ...                │ ...            │      ... │                ... │
└────────────────────┴────────────────┴──────────┴────────────────────┘
```

---

## 📊 Output Files & Interpretation

After running `magcluster run <inputs> -o results/`, the output folder contains:

```text
results/
├── {genome_name}_mgc_genes.tsv              # Comprehensive tabular annotations per genome
├── {genome_name}_magnetosome_proteins.faa   # Extracted magnetosome protein sequences
├── MGC_align.html                           # Interactive Clinker alignment & figure builder
├── clusters_gbk/                            # Sliced GenBank files for each cluster
│   └── {cluster_id}.gbk
└── figures/                                 # (Optional) SVG figures if --static-plot is set
    └── {cluster_id}.svg
```

### 1. Tabular Report (`{genome_name}_mgc_genes.tsv`)

A tab-delimited file documenting every magnetosome hit detected in the genome:

| Column | Description |
| :--- | :--- |
| `cluster_id` | MGC cluster identifier (or `unclustered` for lone magnetosome hits outside clusters). |
| `genome` | Genome prefix identifier. |
| `gene_id` | Locus tag or predicted gene ID. |
| `contig` | Contig or scaffold identifier. |
| `start` | 1-based start coordinate. |
| `end` | 1-based end coordinate. |
| `strand` | Strand orientation (`+` or `-`). |
| `profile` | Best matching HMM profile (e.g., `MamE_full`, `MamE_Cter`, `MamQ_greigite`). |
| `parent_family` | Standardized parent protein family (e.g., `MamE`, `MamQ`). |
| `role` | Functional role classification (`core`, `accessory`, or `iron_context`). |
| `score` | HMM bit score. |
| `evalue` | HMM search E-value. |

### 2. Protein FASTA (`{genome_name}_magnetosome_proteins.faa`)

Amino acid sequences of all identified magnetosome proteins with structured FASTA headers:
```text
>MGMSRV2_RS11765 [MSR-1_genome] cluster=MSR-1_genome_MGC family=MamE profile=MamE_full score=1000.4 evalue=1.97e-301
MSERRPVGRTL...
```

### 3. Interactive Web Visualization (`MGC_align.html`)

Open `MGC_align.html` in any modern web browser to view, explore, and customize your gene cluster alignments:

- **Specific Profile Labels**: Arrows display specific domain or subtype annotations (e.g., `MamE_Cter`, `MamE_Nter`, `MamE_Cter_greigite`, `MamQ_greigite`, `MamO_domain`), while standard full-length proteins are cleanly labeled as `MamE`, `MamO`, `Man5`.
- **Classic Color Palette**:
  - <span style="color:#EA4335">■</span> **`mam`** (`#EA4335`): Core/classical Mam proteins (Google Red)
  - <span style="color:#FBBC05">■</span> **`mms`** (`#FBBC05`): Mms proteins (Google Yellow)
  - <span style="color:#34A853">■</span> **`mad`** (`#34A853`): Desulfobacterota Mad proteins (Google Green)
  - <span style="color:#4385F4">■</span> **`man`** (`#4385F4`): Nitrospirota Man proteins (Google Blue)
  - <span style="color:#FFE066">■</span> **`feo`** (`#FFE066`): Iron transport proteins (Gold Yellow)
  - <span style="color:#8338EC">■</span> **`mae`** (`#8338EC`): Environmental novel Mae proteins (Purple)
  - <span style="color:#06D6A0">■</span> **`maq`** (`#06D6A0`): Novel Maq proteins (Teal)
  - <span style="color:#FF006E">■</span> **`mca`** (`#FF006E`): Novel Mca proteins (Pink)
  - <span style="color:#BDBDBD">■</span> **`Other`** (`#BDBDBD`): Flanking non-magnetosome genes (Muted grey)
- **Vector Export**: Adjust spacing, gene labels, and alignment links in real-time, then export publication-ready **SVG** figures directly from the browser.

---

## 💡 Clinker Customization Options

You can pass any native Clinker option directly to `magcluster run` using `--clinker-opts`:

```bash
# Skip pairwise sequence alignment for fast rendering:
magcluster run *.gbff -o results/ --clinker-opts "-na"

# Retain original file order without automatic clustering reordering:
magcluster run *.gbff -o results/ --clinker-opts "-ufo"

# Combine multiple options (skip alignments + keep file order):
magcluster run *.gbff -o results/ --clinker-opts "-na -ufo"
```

To view all supported Clinker options directly in your terminal, run:
```bash
magcluster run --clinker-help
```

---

## 🔬 Core Magnetosome Defining Families

MagCluster defines an authentic Magnetosome Gene Cluster based on synteny and the presence of magnetosome protein families. The 10 core defining families are:

$$\text{MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ}$$

---

## 📖 Citation

If you use MagCluster in your research, please cite:

> Runjia Ji, Wensi Zhang, Yongxin Pan, and Wei Lin. (2022).  
> **MagCluster: a Tool for Identification, Annotation, and Visualization of Magnetosome Gene Clusters.**  
> *Microbiology Resource Announcements*, 11:e01031-21.  
> DOI: [10.1128/mra.01031-21](https://doi.org/10.1128/mra.01031-21)

Please also acknowledge the core dependencies powering MagCluster v1.0:
- **Pyrodigal**: Larralde (2022). *Journal of Open Source Software*, 7(72), 4296. [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **PyHMMER**: Larralde (2023). *Bioinformatics*, 39(5), btad211. [doi:10.1093/bioinformatics/btad211](https://doi.org/10.1093/bioinformatics/btad211)
- **Clinker**: Gilchrist & Chooi (2021). *Bioinformatics*, 37(16), 2473–2475. [doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007)

---

## 📬 Contact & Support

For questions, issues, or feature requests:
- **GitHub Issues**: [Open an issue on GitHub](https://github.com/RunJiaJi/magcluster/issues)
- **Email**: Runjia Ji (`jirunjia@gmail.com`) or Wei Lin (`weilin@mail.iggcas.ac.cn`)
