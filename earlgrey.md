# earlgrey.md

This document describes a minimal workflow to run **earlgrey** for detecting TEs.

---

## Overview
1. Install EarlGrey
2. Run EarlGrey

---

## 1. Install EarlGrey

cd /data2/nojiri/earlgrey_wd

conda create -n earlgrey -c conda-forge -c bioconda earlgrey=6.3.1

conda activate earlgrey

#earlgreyインストール時にDfamのpartitionを設定する必要がある
#0(root)と9(metazoa)を選択
./configure_dfam39.sh


## 2. Run EarlGrey
earlGrey -g ./H_vulgaris_AEP.fa -s AEP -o ./AEPearlGreyOutputs -t 40
earlGrey -g ./uc.FINAL.fasta -s JNIG -o ./JNIGearlGreyOutputs -t 40
