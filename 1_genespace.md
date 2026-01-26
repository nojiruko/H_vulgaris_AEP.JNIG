# 1_genespace.md

This document describes a minimal workflow to run **GENESPACE** for comparative genomics / synteny analysis.

---

## Overview

1. Install required tools  
2. Prepare GENESPACE input files  
3. Run GENESPACE and generate plots  

---

## 1. Tool Installation

```bash

# 1.1 Create a conda environment and install dependencies
cd genespace_wd

conda create -n genespace
conda activate genespace

conda install -c bioconda orthofinder=2.5.5
conda install -c conda-forge mamba -y
mamba install -c bioconda agat gffread -y

#1.2 Install MCScanX
wget https://github.com/wyp1125/MCScanX/archive/refs/tags/v1.0.0.tar.gz
tar -xzf v1.0.0.tar.gz

cd MCScanX-1.0.0
make
```


## 2. Prepare GENESPACE input files  

```bash
#Keep only the longest isoform per gene
agat_sp_keep_longest_isoform.pl --gff H_vulgaris_AEP.gff3 -o H_vulgaris_AEP.longestIso.gff3
gffread H_vulgaris_AEP.longestIso.gff3 -g H_vulgaris_AEP.fa -y H_vulgaris_AEP.prot.longestIso.fa

agat_sp_keep_longest_isoform.pl --gff jnig15_withGENE_withCodons.gtf -o uc.FINAL.longestIso.gff3
gffread uc.FINAL.longestIso.gff3 -g uc.FINAL.fasta -y uc.FINAL.prot.longestIso.fa

#Generate BED files from annotation
awk -F'\t' '!/^#/ && $3=="mRNA"{
  id=$9; sub(/^.*ID=/,"",id); sub(/;.*$/,"",id);
  print $1,$4-1,$5,id
}' OFS='\t' H_vulgaris_AEP.longestIso.gff3 > H_vulgaris_AEP.prot.longestIso.bed

awk -F'\t' '!/^#/ && $3=="transcript"{
  id=$9; sub(/^.*ID=/,"",id); sub(/;.*$/,"",id);
  print $1,$4-1,$5,id
}' OFS='\t' uc.FINAL.longestIso.gff3 > uc.FINAL.prot.longestIso.bed

#Organize files into directories
mkdir peptide
mkdir bed

cp ./*.bed ./bed/
cp ./*.prot.longestIso.fa ./peptide/

#Replace dots (.) with underscores (_) in filenames/IDs
cd ./bed
sed -i 's/\./_/g' *.bed
cd ..

cd ./peptide
sed -i 's/\./_/g' *.fa
cd ..
```

## 2. Prepare GENESPACE input files  

```bash
#Install GENESPACE in R
install.packages("remotes")
remotes::install_github("jtlovell/GENESPACE", upgrade = "never")

#Run GENESPACE
library(GENESPACE)

wd <- "/data2/nojiri/genespace_wd"
path2mcscanx <- "/data2/nojiri/genespace_wd/MCScanX-1.0.0/"

gpar <- init_genespace(wd = wd, path2mcscanx = path2mcscanx)
out <- run_genespace(gpar, overwrite = TRUE)
