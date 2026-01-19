# Whole-genome comparative analysis revealing signatures of aging-related genes in immortal *Hydra vulgaris*

This repository contains step-by-step descriptions of all analyses and associated code related to the following manuscript:

> **Kaoruko Nojiri†, Koryu Kin†, Akimasa Someya, Tetsuo Kon, Koto Kon-Nanjo, Hiroshi Shimizu, Kazuharu Arakawa, and Etsuo A. Susaki**  
> **Whole-genome comparative analysis revealing signatures of aging-related genes in immortal *Hydra vulgaris***  
> (in preparation / 2026).  
> doi: [***](https://doi.org/***)  
> PMID: [***](https://pubmed.ncbi.nlm.nih.gov/***/)

---

## Structure and purpose of this repository

This repository provides **markdown documents** that describe each computational component of the analyses presented in the manuscript. 
Each markdown document includes all code used for the corresponding analysis, along with explanatory text describing the purpose and rationale of the code.

Specifically, the following markdown documents　contain scripts to reproduce the figures and tables reported in the manuscript:

- **Orthofinder.md**: Fig. 4, Fig. 5, Supplementary Table 3, Supplementary Table 4, Supplementary Table 5, Supplementary Table 6, Supplementary Fig. 4, Supplementary Fig. 5  
- **plotsr.md**: Fig. 2D, Supplementary Fig. 2A  
- **quarTeT.md**: Fig. 2D (annotation track), Fig. 3A (centromere track)  
- **5mC.md**: Fig. 3A–F, Supplementary Fig. 3  
- **Earlgrey.md**: Table 5, Table 6, Supplementary Fig. 2B, Supplementary Fig. 2C  
- **genespace.md**: Fig. 2C  
- **sprious_transcription.md**: Fig. 3G  


The goal of this repository is to document our methodology in sufficient detail to enable other researchers to reproduce the results reported in the manuscript.

However, the code was originally written to rely on directory/file structures and software path configurations that are specific to the systems on which the analyses were performed.  
Therefore, cloning this repository alone does not fully reproduce the original file organization.

In addition, due to file size limitations, we are not able to host all files required for every analysis on GitHub.  
As a result, users will need to modify file paths within the scripts as needed and download additional required files from external sources (described below).

---

## Accessing additional required files

All files required to run the analyses described in this repository can be obtained from the following link:  
[https://***]

Raw sequencing data are also available from NCBI under the BioProject accession **PRJNAxxx**.  
In addition, the strain AEP.JNIG *H. vulgaris* genome assembly is available from GenBank under accession **JALDPZxxx**.
