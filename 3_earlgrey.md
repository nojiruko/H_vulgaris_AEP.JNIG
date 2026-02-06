# 3_earlgrey

This document describes the script used to generate Table 5, Table 6, Supp Fig. 2B and Supp Fig. 2C.

## Overview
1. Install EarlGrey
2. Run EarlGrey

## 1. Install EarlGrey
```bash
cd /data2/nojiri/earlgrey_wd

conda create -n earlgrey -c conda-forge -c bioconda earlgrey=6.3.1

conda activate earlgrey

# When installing EarlGrey, it is necessary to specify Dfam partitions.
# Select partitions 0 (root) and 9 (metazoa).
./configure_dfam39.sh
```

## 2. Run EarlGrey
```bash
earlGrey -g ./H_vulgaris_AEP.fa -s AEP -o ./AEPearlGreyOutputs -t 40
earlGrey -g ./uc.FINAL.fasta -s JNIG -o ./JNIGearlGreyOutputs -t 40
```
