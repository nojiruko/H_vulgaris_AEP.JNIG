# Spurious Transcription / Exon coverage analysis

This repository contains a pipeline to:

1. Create exon BED files annotated with exon number and total exon count (strand-aware)
2. Compute mean coverage per exon using RNA-seq mapping + `bedtools map`
3. Visualize exon coverage ratios grouped by methylation status and gene conservation

---

## Requirements

### Tools
- conda (recommended)
- awk, sort, comm
- bedtools
- samtools
- hisat2
- python >= 3.8
- python packages: pandas, numpy, matplotlib

---

## 1) Create strand-aware exon BED with exon number and total exon count

```bash
cd SpriousTranscription_wd

conda activate ntLink
conda install -c bioconda bedtools
