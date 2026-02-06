# 7_SpuriousTranscription
This document describes the script used to generate Fig. 3A–F and Supp Fig. 3G.

## Overview
1. Create exon BED files annotated with exon number and total exon count
2. Compute mean coverage per exon using RNA-seq
3. Visualize exon coverage ratios grouped by methylation status and gene conservation


## 1. Create strand-aware exon BED with exon number and total exon count

```bash
cd SpriousTranscription_wd

conda activate ntLink
conda install -c bioconda bedtools

#Extract transcript_id and gene_id as new columns
awk -F '\t' '{
    tid=""; gid="";
    if (match($9, /transcript_id "([^"]+)"/, t)) tid=t[1];
    if (match($9, /gene_id "([^"]+)"/, g))       gid=g[1];
    print $1, $2, $3, $4, $5, $6, $7, $8, tid, gid;
}' OFS='\t' jnig15_withGENE_withCodons.gtf > jnig_15chrom.txt

# Identify genes with multiple isoforms (.t2–.t9)
# The maximum number of isoforms per JNIG gene was 9
awk -F'\t' '$9 ~ /\.t[2-9]$/ { print $10 }' jnig_15chrom.txt | sort -u > gene_list_t2_to_t9.txt
wc -l gene_list_t2_to_t9.txt
# 1713

#Identify all genes
cut -f10 jnig_15chrom.txt | sort -u > gene_list.txt
wc -l gene_list.txt
# 28516

#Extract single-isoform genes
comm -23 gene_list.txt gene_list_t2_to_t9.txt > single_isoform_gene_list.txt
wc -l single_isoform_gene_list.txt
# 26803

#Filter annotation table by single-isoform genes
awk 'NR==FNR{keep[$1]; next} ($10 in keep)' single_isoform_gene_list.txt jnig_15chrom.txt \
    > braker_single_isoform.txt

#+ strand: exon numbering and total exon count
awk -F'\t' '$3=="exon" && $7=="+"' braker_single_isoform.txt \
  | sort -t$'\t' -k10,10 -k4,4n \
  > plus_exon.txt

awk -F'\t' '{count[$10]++; print $0, count[$10];}' OFS='\t' plus_exon.txt \
  > plus_exon_number.txt

awk -F'\t' '
    NR==FNR {
        gene_id = $10
        total[gene_id]++
        next
    }
    {
        gene_id = $10
        print $0, total[gene_id]
    }
' OFS='\t' plus_exon_number.txt plus_exon_number.txt \
  > plus_exon_number_with_total.txt

#- strand: exon numbering and total exon count
awk -F'\t' '$3=="exon" && $7=="-"' braker_single_isoform.txt \
  | sort -t$'\t' -k10,10 -k5,5nr \
  > minus_exon.txt

awk -F'\t' '{count[$10]++; print $0, count[$10];}' OFS='\t' minus_exon.txt \
  > minus_exon_number.txt

awk -F'\t' '
    NR==FNR {
        gene_id = $10
        total[gene_id]++
        next
    }
    {
        gene_id = $10
        print $0, total[gene_id]
    }
' OFS='\t' minus_exon_number.txt minus_exon_number.txt \
  > minus_exon_number_with_total.txt

#Merge strands and build BED (genes with >= 6 exons)
cat plus_exon_number_with_total.txt minus_exon_number_with_total.txt > exon_number.txt

awk -F'\t' '
    $12 >= 6 {
        $4 = $4 - 1
        print $1, $4, $5, $6, $7, $8, $9, $10, $11, $12
    }
' OFS='\t' exon_number.txt > exon_number_over6.bed

sort -k1,1 -k2,2n -k3,3n exon_number_over6.bed > exon_number_over6.sorted.bed
```


## 2. Compute exon coverage
```bash
#Mapping RNA-seq reads (HISAT2)
hisat2-build uc.FINAL.fasta hydra_index

hisat2 -p 90 -x hydra_index \
  -1 S2_RNA_Cont_1.fq.gz \
  -2 S2_RNA_Cont_2.fq.gz \
  > mapped.sam

#Convert SAM -> sorted BAM
samtools view -@ 90 -bS mapped.sam > mapped.bam
samtools sort -@ 90 -o mapped.sorted.bam mapped.bam

#Per-base depth
bedtools genomecov -ibam mapped.sorted.bam -d > cont.depth.txt

#Convert depth to BED4 (chr, start0, end1, depth)
awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' cont.depth.txt > cont.depth.bed
sort -k1,1 -k2,2n -k3,3n cont.depth.bed > cont.depth.sorted.bed

#Mean coverage per exon
bedtools map \
  -a exon_number_over6.sorted.bed \
  -b cont.depth.sorted.bed \
  -c 4 \
  -o mean \
  > single_isoform_exon_over6_count_cov.bed
```

## 3. Plot exon coverage ratios
```bash

#python
import pandas as pd
import numpy as np
import re

df = pd.read_csv("single_isoform_exon_over6_count_cov.bed", sep="\t", header=None)

# Rename columns
df = df.rename(columns={df.columns[0]: "chr", df.columns[1]: "start", df.columns[2]: "end", df.columns[4]: "strand", df.columns[7]: "gene_id", df.columns[6]: "transcript_id", df.columns[8]: "exon_number", df.columns[10]: "cov"})
df["exon_number"] = pd.to_numeric(df["exon_number"], errors="coerce")
df["cov"] = pd.to_numeric(df["cov"], errors="coerce")

df_plus  = df[df["strand"] == "+"].copy()
df_plus_gene_cov = df_plus[df_plus["exon_number"] == 1][["gene_id", "cov"]].copy()
df_minus = df[df["strand"] == "-"].copy()
df_minus_gene_cov = df_minus[df_minus["exon_number"] == 1][["gene_id", "cov"]].copy()
df_plus_gene_cov = df_plus_gene_cov.rename(columns={"cov": "cov_exon1"})
df_minus_gene_cov = df_minus_gene_cov.rename(columns={"cov": "cov_exon1"})


df_plus = df_plus.merge(df_plus_gene_cov, how="left", left_on="gene_id", right_on="gene_id")
df_minus = df_minus.merge(df_minus_gene_cov, how="left", left_on="gene_id", right_on="gene_id")
df_plus["ratio"] = df_plus["cov"] / df_plus["cov_exon1"].replace(0, np.nan)
df_minus["ratio"] = df_minus["cov"] / df_minus["cov_exon1"].replace(0, np.nan)

df_all = pd.concat([df_plus, df_minus], axis=0, ignore_index=True)

methylated = pd.read_csv("/data2/nojiri/5mC_wd/methylated_gene.txt", header=0)
methylated = df["name"].str.extract(r"gene_id=([^;]+)", expand=False)
methylated_set = set(methylated)
df_all["methylation_status"] = df_all["gene_id"].apply(lambda x: "methylated" if x in methylated_set else "nonmethylated")

conserved = [line.strip() for line in open("/data2/nojiri/5mC_wd/all_common_og_lists.txt")]
conserved = [re.sub(r"^file_\d+_file_\d+_", "", x) for x in conserved]
conserved_set = set(conserved)
df_all["conserved_gene"] = df_all["gene_id"].apply(lambda x: "conserved" if x in conserved_set else "nonconserved")

# Create labels for conserved_gene × methylation_status
df_all['group'] = df_all['conserved_gene'] + "_" + df_all['methylation_status']
df_all["group"].unique() 

df_all = df_all[df_all['ratio'] > 0]
df_all['log_ratio'] = np.log(df_all['ratio'])

groups = ['nonconserved_nonmethylated', 'nonconserved_methylated', 'conserved_nonmethylated', 'conserved_methylated']
results = []


for exon in [1, 2, 3, 4, 5, 6]:
    g1 = df_all.loc[(df_all["exon_number"] == exon) & (df_all["group"] == groups[0]), "log_ratio"].dropna()
    g2 = df_all.loc[(df_all["exon_number"] == exon) & (df_all["group"] == groups[1]), "log_ratio"].dropna()
    g3 = df_all.loc[(df_all["exon_number"] == exon) & (df_all["group"] == groups[2]), "log_ratio"].dropna()
    g4 = df_all.loc[(df_all["exon_number"] == exon) & (df_all["group"] == groups[3]), "log_ratio"].dropna()
    # means
    mean_g1 = g1.mean()
    mean_g2 = g2.mean()
    mean_g3 = g3.mean()
    mean_g4 = g4.mean()
    results.append([exon, len(g1), len(g2), len(g3), len(g4), mean_g1, mean_g2, mean_g3, mean_g4])

df_results = pd.DataFrame(results, columns=["exon", "len_g1", "len_g2", "len_g3", "len_g4", "mean_g1", "mean_g2", "mean_g3", "mean_g4"])

#Note: df_results
>>> df_results
   exon  len_g1  len_g2  len_g3  len_g4   mean_g1   mean_g2   mean_g3   mean_g4
0     1    2960    3230     356    1090  0.000000  0.000000  0.000000  0.000000
1     2    2917    3213     354    1089  0.348294  0.265944  0.532704  0.294763
2     3    2875    3211     350    1089  0.480233  0.316983  0.759117  0.331617
3     4    2864    3213     354    1090  0.546983  0.340703  0.803252  0.333202
4     5    2851    3214     351    1090  0.576894  0.344364  0.812908  0.331614
5     6    2820    3212     352    1090  0.578766  0.343732  0.819413  0.313879

# Plot a line graph
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
groups = ["mean_g1", "mean_g2", "mean_g3", "mean_g4"]
for g in groups:
    plt.plot(df_results["exon"], df_results[g], marker="o", label=g)

plt.xlabel("Exon number")
plt.ylabel("Mean value")
plt.title("Mean values per exon for 4 groups")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("mean_exon_cov_plot.svg", format="svg")
```
