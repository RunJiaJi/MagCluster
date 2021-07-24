# magcluster
Magnetosome gene cluster annotation, screening and mapping tool
## Installation
We recommend creating a ***new environment*** for the magcluster release being installed through conda.
```bash
wget https://github.com/RunJiaJi/magcluster/releases/download/0.1.4/magcluster-0.1.4.yml
```
```bash
conda env create -n magcluster --file magcluster-0.1.4.yml
```
```bash
# OPTIONAL CLEANUP
rm magcluster-0.1.4.yml
```
Alternatively, you can install magcluster through pip in an existing environment. In this way, please make sure you have prokka installed.
```bash
#Prokka installation
conda install -c conda-forge -c bioconda -c defaults prokka=1.13.4
```
```bash
pip install magcluster
```

## General Usage
```bash
usage: magcluster [options]

Magnetosome Gene Cluster Analysis Tool

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show magcluster version number and exit

Options:
  {prokka,mgc_screen,clinker}
    prokka              Genome annotation with Prokka
    mgc_screen          Magnetosome gene cluster screening with magscreen
    clinker             Magnetosome gene cluster mapping with Clinker

General usage
-------------
Magnetosome gene annotation:
  $ magcluster prokka XXX.fa

Magnetosome gene screen:
  $ magcluster mgc_screen XXX.gbk

Magnetosome gene cluster mapping:
  $ magcluster clinker XXX_screened.gbk

Runjia, 2021
```
