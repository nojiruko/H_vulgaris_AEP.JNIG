# 2_plotsr

This document describes a minimal workflow to run **plotsr** for stractual variation.

---

## Overview
1. Detect variants using SyRI
2. Visualize using plotsr

## 1. Detect variants using SyRI
```bash
cd plotsr_wd

conda activate ntLink
conda install -c bioconda syri

minimap2 -ax asm5 -t 40 --eqx H_vulgaris_AEP.fa uc.FINAL.fasta | samtools sort -O BAM - > AEP_JNIG.bam

samtools index AEP_JNIG.bam

syri -c AEP_JNIG.bam -q uc.FINAL.fasta -r H_vulgaris_AEP.fa -F B --prefix AEP_JNIG
```

## 2. Visualize using plotsr
```bash
conda install -c bioconda plotsr 

#geneの位置
awk -F'\t' '$3=="gene" {print $1, $4-1, $5}' OFS='\t' H_vulgaris_AEP.gff3 > H_vulgaris_AEP_gene.bed

#repeatの位置
WS2の/data2/nojiri/earlgrey_wd/AEPearlGreyOutputs/AEP_EarlGrey/AEP_summaryFiles/AEP.filteredRepeats.bedをコピー

#centromere bedの作成
#/data2/nojiri/quarTeT_wd/aep_quaTeT/Candidatesを見てcentromere.bedを作成

#genomes.txtの作成(nanoコマンドで以下を記載)
#file	name	tags
H_vulgaris_AEP.fa	H_vulgaris_AEP	lw:1.5
uc.FINAL.fasta	H_vulgaris_JNIG	lw:1.5

#tracks.txtの作成(nanoコマンドで以下を記載)
H_vulgaris_AEP_gene.bed	Genes	bw:5000;ft:bed;lc:crimson;lw:0.01;bc:mistyrose;ti:1;ns:1
AEP.filteredRepeats.bed	Repeat	bw:5000;ft:bed;lc:dodgerblue;lw:0.01;bc:aliceblue;ti:2;ns:1
centromere.bed	Centromere	bw:5000;ft:bed;lc:forestgreen;lw:0.01;bc:honeydew;ti:3;ns:1

plotsr --sr AEP_JNIGsyri.out --genomes genomes.txt --tracks tracks.txt -o plotsr.tracks.pdf
```
