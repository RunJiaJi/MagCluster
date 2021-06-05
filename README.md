# magcluster
Magnetosome gene cluster annotation, screening and mapping tool
## Installation
We recommend creating a new environment for the magcluster release being installed through conda.
```bash
wget https://github.com/RunJiaJi/magcluster/releases/download/0.0.6/magcluster-0.0.6.yml

conda env create -n magcluster --file magcluster-0.0.6.yml

# OPTIONAL CLEANUP
rm magcluster-0.0.6.yml
```

## General Usage
Magnetosome gene annotation:  
```
magcluster maga XXX.fa
```
  
Magnetosome gene screen:  
```
magcluster magsc XXX.faa XXX.gbk
```
  
Magnetosome gene cluster mapping:  
```
magcluster magm XXX_screened.gbk
```
