# Comparative genomics reveals signatures of distinct metabolic strategies and gene loss associated with _Hydra_ immortality

This repository contains step-by-step descriptions of all analyses and associated code related to the following manuscript:

> **Kaoruko Nojiri, Koryu Kin, Akimasa Someya, Tetsuo Kon, Koto Kon-Nanjo, Hiroshi Shimizu, Kazuharu Arakawa, and Etsuo A. Susaki**  
> **Comparative genomics reveals signatures of distinct metabolic strategies and gene loss associated with _Hydra_ immortality**  
> (in preparation 2026).  
> doi: ***
> PMID: ***

## Structure of this repository

This repository contains markdown documentation for the computational analyses presented in this manuscript.

Specifically, the following markdown documents contain scripts to reproduce the figures and tables reported in the manuscript:

- **1_genespace.md**: Fig. 2C
- **2_plotsr.md**: Fig. 2D, Supp Fig. 2A 
- **3_Earlgrey.md**: Table 5, Table 6, Supp Fig. 2B, Supp Fig. 2C
- **4_quarTeT.md**: Fig. 2D (annotation track), Fig. 3A (centromere track)
- **5_karyoplotR.md**: Fig. 3A
- **6_5mC.md**: Fig. 3A–F, Supp Fig. 3  
- **7_sprioustranscription.md**: Fig. 3G
- **8_Orthofinder.md**: Fig. 4, Fig. 5, Supp Table 3, Supp Table 4, Supp Table 5, Supp Table 6, Supp Fig. 4, Supp Fig. 5   


The code was originally written to rely on directory/file structures and software path configurations that are specific to the systems on which the analyses were performed.  
Therefore, cloning this repository alone does not fully reproduce the original file organization.

In addition, due to file size limitations, we are not able to host all files required for every analysis on GitHub.  
As a result, users will need to modify file paths within the scripts as needed and download additional required files from external sources (described below).

## Accessing additional required files

> Nanopore and Hi-C sequencing data have been deposited in ***** under accession numbers ***** and *****, respectively.
> The genome assembly has been deposited in NCBI with BioProject ID: PRJNA1413991 and BioSample accession: SAMN54890366.

## System requirements
```bash
### System
- OS: Linux (Ubuntu 22.04.5 LTS)
- Conda: 24.4.0
- Mamba: 2.4.0

### Core tools
- HISAT2: 2.2.1
- MCScanX: 1.0.0
- GENESPACE: 1.2.3
- OrthoFinder: 2.5.5
- SyRI: 1.7.1
- minimap2: 2.28-r1209
- samtools: 1.23
- bedtools: 2.31.1
- bedops: 2.4.42
- gtf2bed: 2.4.42
- bedGraphToBigWig: 2.10
- EarlGrey: 6.3.3
- EDTA: 2.2.2
- quarTeT: 1.2.5
- dorado: 0.7.1
- modkit: 0.4.3
- bedtools: 2.31.1
- bedops: 2.4.42
- gtf2bed: 2.4.42

### R environment
- R: 4.4.2
- Packages:
  - karyoploteR: 1.32.0
  - GenomicRanges: 1.58.0
  - rtracklayer: 1.66.0
  - GenomicFeatures: 1.58.0
  - BSgenome: 1.74.0

### Python environment
- Python: 3.9.16
- Packages:
  - pandas: 2.2.3
  - numpy: 1.26.4
  - matplotlib: 3.9.2
  - seaborn: 0.13.2
  - scipy: 1.13.1
  - AGAT: 1.5.1
```
