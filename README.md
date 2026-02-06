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
- **2_plotsr.md**: Fig. 2D, Supplementary Fig. 2A 
- **3_Earlgrey.md**: Table 5, Table 6, Supplementary Fig. 2B, Supplementary Fig. 2C
- **4_quarTeT.md**: Fig. 2D (annotation track), Fig. 3A (centromere track)
- **5_karyoplotR.md**: Fig. 3A
- **6_5mC.md**: Fig. 3A–F, Supplementary Fig. 3  
- **7_sprioustranscription.md**: Fig. 3G
- **8_Orthofinder.md**: Fig. 4, Fig. 5, Supplementary Table 3, Supplementary Table 4, Supplementary Table 5, Supplementary Table 6, Supplementary Fig. 4, Supplementary Fig. 5   


The code was originally written to rely on directory/file structures and software path configurations that are specific to the systems on which the analyses were performed.  
Therefore, cloning this repository alone does not fully reproduce the original file organization.

In addition, due to file size limitations, we are not able to host all files required for every analysis on GitHub.  
As a result, users will need to modify file paths within the scripts as needed and download additional required files from external sources (described below).

## Accessing additional required files

Nanopore and Hi-C sequencing data have been deposited in ***** under accession numbers ***** and *****, respectively. The genome assembly has been deposited in NCBI with BioProject ID: PRJNA1413991 and BioSample accession: SAMN54890366.
