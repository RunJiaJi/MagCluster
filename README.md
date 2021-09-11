# MagCluster
MagCluster is a tool for identification, annotation and visualization of magnetosome gene clusters (MGCs) from genomes of magnetotactic bacteria (MTB).

## Table of Contents
- [Installation](#installation)
  - [Conda](#conda)
  - [Pip](#pip)
- [Usage](#usage)
  - [Genome annotation](#genome-annotation)
  - [MGCs screening](#MGCs-screening)
  - [MGCs alignment and mapping](#MGCs-alignment-and-mapping)
- [Getting Started](https://github.com/RunJiaJi/magcluster/blob/main/Getting_started.md)
- [Citation](#Citation)
- [Contact us](#contact-us)
---

## Installation
MagCluster requires a working **[Conda](https://www.anaconda.com/products/individual)** installation.

### Conda
MagCluster can be installed through conda. We recommend creating a ***new environment*** for MagCluster.
```bash
# Create magcluster environment
conda create -n magcluster

# Activate magcluster environment
conda activate magcluster

# Install MagCluster through bioconda channel
conda install -c conda-forge -c bioconda -c defaults magcluster

# Check for the usage of MagCluster
magcluster -h
```
### Pip
Alternatively, you can install magcluster through pip in an existing environment. In this way, please make sure you have [prokka](https://github.com/tseemann/prokka) installed.

```bash
# Install MagCluster through pip
pip install magcluster

# Check for the usage of MagCluster
magcluster -h
```

## Usage


MagCluster comprises three modules for MGCs batch processing: 
(i) MTB genome annotation with **[Prokka](https://github.com/tseemann/prokka)**
(ii) MGCs screening with **MGC_Screen**
(iii) MGCs mapping with **[Clinker](https://github.com/gamcil/clinker)**


```bash
usage: magcluster [options]

Options:
  {prokka,mgc_screen,clinker}
    prokka              Genome annotation with Prokka
    mgc_screen          Magnetosome gene cluster screening with MGC_Screen
    clinker             Magnetosome gene cluster mapping with Clinker
```
#### Genome annotation
MagCluster allows users to input **multiple genome files** or **genome-containing folder(s)** in one command for batch annotation. The general usage is same as Prokka yet some parameters are set with default value for MTB genome batch annotation.

To avoid confusion, the name of each genome is used as the output folder’s name (**--outdir** GENOME_NAME), output files’ prefix (**--prefix** GENOME_NAME), and GenBank file’s locus_tag (**--locustag** GENOME_NAME) by default. The ‘**--compliant**’ parameter is also used by default to ensure the standard GenBank files. 

For MGCs annotation, we provide a **[reference MGCs file](https://github.com/RunJiaJi/magcluster/releases/download/v1.0/Magnetosome_protein_data.fasta.faa)** containing magnetosome protein sequences from representative MTB strains which is used by default. The value of '**--evalue**' is recommended to set to 1e-05.
```bash
example usage: 

# MGCs annotation with multiple MTB genomes as input
$ magcluster prokka --evalue 1e-05 --proteins Magnetosome_protein_data.fasta MTB_genome1.fasta MTB_genome2.fasta MTB_genome3.fasta

# MGCs annotation with MTB genomes containing folder as input
$ magcluster prokka --evalue 1e-05 --proteins Magnetosome_protein_data.fasta /MTB_genomes_folder
```
#### MGCs screening
MGC_Screen module retrieves MGC-containing contigs/scaffolds in GenBank files. As magnetosome genes are always physically clustered in MTB genomes, MGC_Screen identify MGC based on the number of magnetosome genes gathered. 
Three parameters involved in MGC screening: '**--contiglength**', '**--windowsize**' and '**--threshold**' (see below). Users can adjust them according to needs. 
For each genome, MGC_Screen produces two files as output: a **GenBank file of MGCs containing contigs** and a **csv file summarizing all magnetosome proteins sequences**.
```bash

usage: magcluster mgc_screen [-h] [-l CONTIGLENGTH] [-win WINDOWSIZE] [-th THRESHOLD] [-o OUTDIR] gbkfile [gbkfile ...]

positional arguments:
  gbkfile               .gbk/.gbf files to analyzed. Multiple files or files-containing folder is acceptable.

optional arguments:
  -h, --help            show this help message and exit
  -l CONTIGLENGTH, --contiglength CONTIGLENGTH
                        The minimum size of a contig for screening (default '2,000 bp')
  -w WINDOWSIZE, --windowsize WINDOWSIZE
                        The window size in the text mining of magnetosome proteins (default '10,000 bp')
  -th THRESHOLD, --threshold THRESHOLD
                        The minimum number of magnetosome genes existed in a window size (default '3')
  -o OUTDIR, --outdir OUTDIR
                        Output folder (default 'mgc_screen')
```
```bash
example usage: 

# MGCs screening with multiple GenBank files as input
$ magcluster mgc_screen --threshold 3 --contiglength 2000 --windowsize 10000 file1.gbk file2.gbk file3.gbk

# MGCs screening with GenBank files containing folder as input
$ magcluster mgc_screen --threshold 3 --contiglength 2000 --windowsize 10000 /gbkfiles_folder
```
#### MGCs alignment and mapping
We use [Clinker](https://github.com/gamcil/clinker) for MGCs alignment and visualization. Note that the '**-p**' parameter is used by default to generate an interactive HTML web page where you can modify the MGCs figure and export it as publication-quality file.

```bash
example usage: 

# MGCs screening with multiple GenBank files as input
$ magcluster clinker -p MGC_align.html /MGCs_files_folder/*.gbk
```
## Citation
The manuscript is in preparation.

## Contact us
If you have any questions or suggestions, feel free to contact us.

**jirunjia@gmail.com**
or 
**weilin@mail.iggcas.ac.cn**

