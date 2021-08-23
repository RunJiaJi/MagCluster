# MagCluster
MagCluster is a tool for annotating, screening and mapping magnetosome gene clusters (MGCs) from genomes of magnetotactic bacteria (MTB).
## Installation
MagCluster requires a working **[Conda](https://www.anaconda.com/products/individual)** installation.
We recommend creating a ***new environment*** for the magcluster release being installed through conda.
```bash
wget https://github.com/RunJiaJi/magcluster/releases/download/0.1.6/magcluster-0.1.6.yml

conda env create -n magcluster --file magcluster-0.1.6.yml

# OPTIONAL CLEANUP
rm magcluster-0.1.6.yml
```
Alternatively, you can install magcluster through pip in an existing environment. In this way, please make sure you have prokka installed.
```bash
#Prokka installation
conda install -c conda-forge -c bioconda -c defaults prokka
```
```bash
pip install magcluster
```

## Usage
MagCluster comprises three modules for MGCs batch processing: 
(i) MTB genome annotation with **[Prokka](https://github.com/tseemann/prokka)**
(ii) magnetosome gene cluster screening with **Mgc_Screen**
(iii) MGCs mapping with **[Clinker](https://github.com/gamcil/clinker)**


```bash
usage: magcluster [options]

Options:
  {prokka,mgc_screen,clinker}
    prokka              Genome annotation with Prokka
    mgc_screen          Magnetosome gene cluster screening with magscreen
    clinker             Magnetosome gene cluster mapping with Clinker
```
#### MTB genome annotation
MagCluster allows users to input **multiple genome files** or **genome-containing folder(s)** in one command for batch annotation. The general usage is same as Prokka yet some parameters are set with default value for MTB genome batch annotation.

To avoid confusion, the name of each genome is used as the output folder’s name (**--outdir** GENOME_NAME), output files’ prefix (**--prefix** GENOME_NAME), and GenBank file’s locus_tag (**--locustag** GENOME_NAME) by default. The ‘**--compliant**’ parameter is also used by default to ensure the standard GenBank files. 

For MGCs annotation, we provide a **[reference MGCs file](https://github.com/RunJiaJi/magcluster/releases/download/v1.0/Magnetosome_protein_data.fasta.faa)** containing magnetosome protein sequences of novel MTB strains which is highly recommended to use with ‘**--proteins**’ parameter. The value of '**--evalue**' is recommended to set to 1e-05.
```bash
example usage: 

# MGCs annotation with multiple MTB genomes as input
$ magcluster prokka --evalue 1e-05 --proteins Magnetosome_protein_data.fasta MTB_genome1.fasta MTB_genome2.fasta MTB_genome3.fasta

# MGCs annotation with MTB genomes containing folder as input
$ magcluster prokka --evalue 1e-05 --proteins Magnetosome_protein_data.fasta /MTB_genomes_folder
```
#### MGCs screening
Mgc_Screen module retrieves MGC_containing contigs in GenBank files. As magnetosome genes are always physically clustered on genome, Mgc_Screen identify MGC based on the number of magnetosome genes gathered. 
Three parameters involved in MGC screening: '**--minlength**', '**--maxlength**',  '**--threshold**' (see below). Users can adjust them according to needs. 
Mgc_screen produces two files as output: a **GenBank file of MGCs containing contigs** and a **csv file summarizing all magnetosome proteins sequences**.
```bash
usage: magcluster mgc_screen [-h] [-th THRESHOLD] [-o OUTDIR] [-min MINLENGTH] [-max MAXLENGTH] gbkfile [gbkfile ...]

positional arguments:
  gbkfile               .gbk files to analyzed. Multiple files or files-containing folder is acceptable.

optional arguments:
  -h, --help            show this help message and exit
  -th THRESHOLD, --threshold THRESHOLD
                        The minimum number of magnetosome genes in one contig/scaffold to screen (default '2')
  -o OUTDIR, --outdir OUTDIR
                        Output folder (default 'mgc_screen')
  -min MINLENGTH, --minlength MINLENGTH
                        Minimum length of contigs to be considered (default '2000bp')
  -max MAXLENGTH, --maxlength MAXLENGTH
                        Maximum length of contigs containing magnetosome gene (default '10000bp')
```
```bash
example usage: 

# MGCs screening with multiple GenBank files as input
$ magcluster mgc_screen -th 3 -min 2000 -max 10000 file1.gbk file2.gbk file3.gbk

# MGCs screening with GenBank files containing folder as input
$ magcluster mgc_screen -th 3 -min 2000 -max 10000 /gbkfiles_folder
```
#### MGCs alignment and mapping
We use [Clinker](https://github.com/gamcil/clinker) for MGCs alignment and visualization. We recommend using the '**-p**' parameter to generate an interactive html web page where you can modify the MGCs figure and export it as publication-quality file.

```bash
example usage: 

# MGCs screening with multiple GenBank files as input
$ magcluster clinker -o MGCs_alignment_result -p /MGCs_files_folder/*.gbk
```
## Citation
The manuscript is in preparation.

