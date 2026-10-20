# MagCluster <img src='https://raw.githubusercontent.com/RunJiaJi/MagCluster/main/docs/_static/logo_new.svg' align='right' width=220 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.1-orange?style=flat-square)](https://github.com/RunJiaJi/magcluster)

[English](README.md) | [中文说明](README_zh.md)

**MagCluster** is a tool for finding, annotating, and visualizing magnetosome gene clusters (MGCs) in bacterial genomes.

---

## What's New in v1.0

- **One command**: `magcluster run` runs everything from gene prediction to cluster visualization in a single step.
- **No Prokka needed**: gene calling is handled natively by [Pyrodigal](https://github.com/althonos/pyrodigal) (fast, no Perl or external setup).
- **92 HMM profiles**: curated protein HMMs covering 8 magnetosome families (`Mam`, `Mms`, `Mad`, `Man`, `Mae`, `Maq`, `Mca`, `Feo`), searched using [PyHMMER](https://github.com/althonos/pyhmmer).
- **Flexible inputs**: accepts GenBank (`.gbk`/`.gbff`), nucleotide FASTA contigs (`.fasta`/`.fna`/`.fa`), and protein FASTA (`.faa`), including `.gz` compressed files.
- **Per-genome outputs**: writes separate `{genome}_mgc_genes.tsv` and `{genome}_magnetosome_proteins.faa` files for each genome.
- **Run log** (v1.0.1): every run saves `magcluster.log` in the output folder, so you can review the parameters, software versions and per-genome results later. See [CHANGELOG.md](https://github.com/RunJiaJi/MagCluster/blob/main/CHANGELOG.md).


---

## Installation

MagCluster requires Python 3.10 or newer.

### Option A: pip (PyPI)

```bash
pip install magcluster
```

### Option B: conda / mamba (Bioconda)

```bash
conda create -n magcluster -c conda-forge -c bioconda magcluster
conda activate magcluster
```

To add it to an existing environment: `conda install -c conda-forge -c bioconda magcluster`.
List `conda-forge` before `bioconda` as shown. If your conda is configured with strict channel priority, `conda install -c bioconda magcluster` alone can fail to solve.

### Option C: development version from GitHub

```bash
pip install git+https://github.com/RunJiaJi/MagCluster.git
# or, from a local clone:
git clone https://github.com/RunJiaJi/MagCluster.git
cd MagCluster
pip install -e ".[dev]"
```

Verify the installation:
```bash
magcluster --version
magcluster run --help      # -h works too
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

### 4. Protein FASTA (`.faa` / `.pep`)
Prodigal-style headers (`>id # start # end # strand # ...`) keep the genomic coordinates, so clusters are detected and plotted. If your input contains protein sequences without coordinates, MagCluster will annotate magnetosome proteins and output TSV table and FASTA files without plotting the gene cluster:
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

`INPUTS` can be files or folders. In a folder, only files with a genome extension are used (`.gbk`, `.gbff`, `.gbf`, `.gb`, `.fa`, `.fasta`, `.fna`, `.fsa`, `.fsa_nt`, `.faa`, `.pep`, each optionally `.gz`). Other files (e.g. `README.txt`) and hidden files are skipped with a warning. A file given explicitly is used even without an extension; its format is detected from the first line. GFF3 files are not supported because they contain no sequences. Use the matching GenBank/GBFF, FASTA or Prodigal `.faa` file instead.

| Option | Aliases | Default | Description |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | `magcluster_results` | Output directory (must not be an existing file) |
| `--force` | `-f` | `False` | Allow a non-empty output directory. Outputs of a previous MagCluster run there (tables, FASTA, `MGC_align.html`, `clusters_gbk/`, `figures/`, `magcluster.log`) are deleted first; other files are kept |
| `--threads` | `-t` | `0` | Number of CPU threads (`0` = all cores; values above the available cores are capped) |
| `--evalue` | `-e` | `1e-5` | HMM search E-value cutoff (≥ 0) |
| `--score` | `-s` | per-profile `trusted_cutoff` (25.0) | Minimum HMM bit score. If omitted, each profile's `trusted_cutoff` from `rules.yaml` is used (25.0 for all 92 built-in profiles). An explicit value replaces it for every profile, including values below 25 |
| `--max-gap` | | `5` | Maximum non-magnetosome genes allowed between hits in a cluster (≥ 0). Hits must also be ≤ 25 kb apart (fixed) |
| `--threshold` | `-th`, `--min-families` | `3` | Minimum number of distinct magnetosome families (any non-Feo family) required to call an MGC (≥ 1) |
| `--flank` | `-un`, `--unmag-num` | `3` | Number of non-magnetosome flanking genes to keep on each side (≥ 0) |
| `--meta / --single` | | Auto | Force metagenome or single-genome mode for gene prediction. Auto uses meta mode for inputs < 100 kb or with > 200 contigs all < 20 kb. `--single` falls back to meta mode (with a warning) when there is too little sequence (< 20 kb) to train. Ignored for inputs that already contain CDS/proteins |
| `--identity` | `-i` | `0.3` | Minimum sequence identity cutoff for Clinker alignments (0.0–1.0) |
| `--clinker-opts` | | `None` | Extra arguments passed to Clinker (e.g. `"-na -ufo"`). Relative output paths such as `-mo matrix.tsv` are written inside the output directory |
| `--clinker-help` | | | Show available Clinker options guide and exit |
| `--static-plot` | | `False` | Also generate static SVG figures in `figures/` |
| `--no-plot` | | `False` | Skip Clinker visualization (tables and GenBank slices only) |
| `--log` | | `<outdir>/magcluster.log` | Write the run log to another file |
| `--help` | `-h` | | Show help and exit |

Global options: `magcluster --version` (`-v`) and `magcluster --help` (`-h`).

Genome IDs (used in output file names, the `genome` column and gene IDs of predicted genes) are taken from the input file name minus `.gz` and the genome extension, e.g. `GCF_000513295.1.fna` → `GCF_000513295.1` and `strainX.v1.fna` → `strainX.v1`. Spaces and special characters are replaced by `_` (e.g. `my strain.fna` → `my_strain`). If two inputs give the same ID (e.g. NCBI Datasets' `genomic.gbff`), the later one gets a suffix (`genomic_2`) and a warning is printed. The same input file listed twice is used once.

Exit status: `0` on success. `1` if any input could not be read, Clinker failed, or the output directory cannot be used. Results of the other genomes are still written. `2` for invalid command-line options.

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
├── MGC_align.html                      # Interactive gene cluster map (open in browser; not with --no-plot)
├── magcluster.log                      # Run log (see below)
├── clusters_gbk/                       # Sliced GenBank file for each cluster
│   └── {cluster_id}.gbk
└── figures/                            # Only with --static-plot
    └── {cluster_id}.svg                # Static cluster map (pyGenomeViz)
```

Genomes without any magnetosome hit produce no table or FASTA file. They are still listed in `magcluster.log`.

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

### 4. Run log (`magcluster.log`)

A plain-text log of the run, overwritten on each run into the same folder (e.g. with `--force`). It contains:

- the start and end time, MagCluster, Python and dependency versions;
- the full command line and every effective parameter;
- the genome ID assigned to each input file;
- per-genome progress, timings and results (genes, hits, MGC coordinates, subclusters with their families, mineral type);
- all warnings and errors, including skipped files.

Use `--log FILE` to write it somewhere else.

### 5. Static figures (`figures/{cluster_id}.svg`)

Written with `--static-plot`: one SVG per MGC (one track per subcluster), drawn with pyGenomeViz.

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

## How an MGC is called

1. Magnetosome HMM hits on the same contig are grouped when at most `--max-gap` (default 5) non-magnetosome genes and at most 25 kb separate consecutive hits. Feo (iron-transport) hits do not start or extend a group.
2. A group is kept if it contains at least `--threshold` (default 3) distinct magnetosome families. **Any non-Feo family counts** (core Mam families as well as accessory families such as Mad, Mms, Man, Mae, Maq, Mca). Feo families never count.
3. Each kept group is extended by `--flank` (default 3) genes on each side. All kept groups of one genome are reported together as one MGC (`{genome}_MGC`) with one subcluster (`loc1`, `loc2`, …) per group.
4. On a complete circular replicon (GenBank `circular` topology, or a FASTA header containing "complete genome", "complete sequence" or "circular"), an MGC whose subclusters lie on both sides of the origin is reported with `start > end` and marked "(spans origin)".

The 10 core defining families shown by `magcluster db` are:

$$\text{MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ}$$

They are listed for reference; MGC calling does not require them.

The keys `clustering_params` and the per-profile `evalue_cutoff` in the bundled `rules.yaml` are informational only and are not read. Use the command-line options above instead.

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
