# quarTeT.md

This document describes a minimal workflow to run **quarTeT** for detecting centromerer.

---

## Overview
1. Predict repeat sequences using EDTA
2. Predict centromeres using quarTeT


## 1. Predict repeat sequences using EDTA

```bash
cd quarTeT_wd

conda create -n edta_env
conda activate edta_env
mamba install -c conda-forge -c bioconda edta
EDTA.pl -v #2.2.2

EDTA.pl --genome ./uc.FINAL.fasta --anno 1 --threads 48
EDTA.pl --genome ./H_vulgaris_AEP.fa --anno 1 --threads 48
```

## 2. Predict centromeres using quarTeT
```bash
mkdir jnig_quarTeT
cd jnig_quarTeT
python3 /data2/nojiri/quarTeT_wd/quarTeT/quartet.py CentroMiner -i /data2/nojiri/quarTeT_wd/uc.FINAL.fasta --TE /data2/nojiri/quarTeT_wd/uc.FINAL.fasta.mod.EDTA.TEanno.gff3 --gene /data2/nojiri/quarTeT_wd/jnig15_withGENE_withCodons.gtf -t 48
cd ..

mkdir aep_quarTeT
cd aep_quarTeT
python3 /data2/nojiri/quarTeT_wd/quarTeT/quartet.py CentroMiner -i /data2/nojiri/quarTeT_wd/H_vulgaris_AEP.fa --TE /data2/nojiri/quarTeT_wd/H_vulgaris_AEP.fa.mod.EDTA.TEanno.gff3 --gene /data2/nojiri/quarTeT_wd/H_vulgaris_AEP.gff3 -t 48

#メモ：AEPのセントロメア#
#/data2/nojiri/quarTeT_wd/aep_quaTeT/Candidatesから抜粋
HVAEP1	47693748	48293151
HVAEP2	4307655	4418075
HVAEP3	49220529	49331983
HVAEP4	11322108	11422229
HVAEP5	36636307	36766067
HVAEP6	33988638	34308723
HVAEP7	14970473	15265460
HVAEP8	33806121	33952685
HVAEP9	27753077	28136273
HVAEP10	24642363	24808711
HVAEP11	29123407	29506872
HVAEP12	44289010	45050735
HVAEP13	44559691	44709459
HVAEP14	7742636	7989881
HVAEP15	23696240	23856110

#メモ：JNIGのセントロメア#
#/data2/nojiri/quarTeT_wd/jnig_quaTeT/Candidatesから抜粋
HiC_scaffold_1	44862403	45129867
HiC_scaffold_2	46070841	46936703
HiC_scaffold_3	48816250	48933293
HiC_scaffold_4	31395376	32345466
HiC_scaffold_5	16262765	16825010
HiC_scaffold_6	36742516	37030628
HiC_scaffold_7	13526727	14357103
HiC_scaffold_8	38305618	38953643
HiC_scaffold_9	35083759	35476109
HiC_scaffold_10	22917584	23826895
HiC_scaffold_11	1733186	1853474
HiC_scaffold_12	57656989	57985783
HiC_scaffold_13	28514404	29741176
HiC_scaffold_14	38441342	38546612
HiC_scaffold_15	23609324	23903717
```

