# 8_orthofinder

This document describes a minimal workflow to run **orthofinder** for comparative genomics analysis.

---

## Overview
1. Download the sequences of 33 species as input for OrthoFinder
2. Generate proteomes from genome sequences
3. Run OrthoFinder
4. Run BLAST
5. Annotate orthogroups
6. Perform enrichment analysis
7. Summarize the number of orthogroups per species and per pathway
8. Identify Ho/Hv-specific orthogroups

## 1. Download the sequences of 33 species as input for OrthoFinder
```bash
#Acropora millepora
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna.gz 
gunzip GCF_013753865.1_Amil_v2.1_genomic.fna.gz
mv GCF_013753865.1_Amil_v2.1_genomic.fna A_millepora.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.gff.gz
gunzip GCF_013753865.1_Amil_v2.1_genomic.gff.gz
mv GCF_013753865.1_Amil_v2.1_genomic.gff A_millepora.gff

#Adineta vaga
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/adineta_vaga/dna/Adineta_vaga.AMS_PRJEB1171_v1.dna.toplevel.fa.gz
gunzip Adineta_vaga.AMS_PRJEB1171_v1.dna.toplevel.fa.gz
mv Adineta_vaga.AMS_PRJEB1171_v1.dna.toplevel.fa A_vaga.fa

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/adineta_vaga/Adineta_vaga.AMS_PRJEB1171_v1.60.gff3.gz
gunzip Adineta_vaga.AMS_PRJEB1171_v1.60.gff3.gz
mv Adineta_vaga.AMS_PRJEB1171_v1.60.gff3 A_vaga.gff3

#Amphimedon queenslandica
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/amphimedon_queenslandica_gca000090795v2rs/dna/Amphimedon_queenslandica_gca000090795v2rs.v1.1.dna.toplevel.fa.gz && gunzip Amphimedon_queenslandica_gca000090795v2rs.v1.1.dna.toplevel.fa.gz && mv Amphimedon_queenslandica_gca000090795v2rs.v1.1.dna.toplevel.fa A_queenslandica.fa && wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/amphimedon_queenslandica_gca000090795v2rs/Amphimedon_queenslandica_gca000090795v2rs.v1.1.60.gff3.gz && gunzip Amphimedon_queenslandica_gca000090795v2rs.v1.1.60.gff3.gz && mv Amphimedon_queenslandica_gca000090795v2rs.v1.1.60.gff3 A_queenslandica.gff3

#Bombyx mori
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/bombyx_mori/dna/Bombyx_mori.Bmori_2016v1.0.dna.toplevel.fa.gz && gunzip Bombyx_mori.Bmori_2016v1.0.dna.toplevel.fa.gz && mv Bombyx_mori.Bmori_2016v1.0.dna.toplevel.fa B_mori.fa

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/bombyx_mori/Bombyx_mori.Bmori_2016v1.0.60.gff3.gz && gunzip Bombyx_mori.Bmori_2016v1.0.60.gff3.gz && mv Bombyx_mori.Bmori_2016v1.0.60.gff3 B_mori.gff3

#Caenorhabditis elegans
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz && gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz && mv Caenorhabditis_elegans.WBcel235.dna.toplevel.fa C_elegans.fa

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.60.gff3.gz && gunzip Caenorhabditis_elegans.WBcel235.60.gff3.gz && mv Caenorhabditis_elegans.WBcel235.60.gff3 C_elegans.gff3 

#Capsaspora owczarzaki
https://juntendo.box.com/s/6234gjolpyjk0anaslgdvb8k4vq4emg9
mv Cowc_genome_v5.fasta C_owczarzaki.fa
mv Cowc_genome_v5.gtf C_owczarzaki.gtf

#Callorhinchus milii
wget https://ftp.ensembl.org/pub/release-113/fasta/callorhinchus_milii/dna/Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna.toplevel.fa.gz && gunzip Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna.toplevel.fa.gz && mv Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna.toplevel.fa C_milii.fa

wget https://ftp.ensembl.org/pub/release-113/gtf/callorhinchus_milii/Callorhinchus_milii.Callorhinchus_milii-6.1.3.113.gtf.gz && gunzip Callorhinchus_milii.Callorhinchus_milii-6.1.3.113.gtf.gz && mv Callorhinchus_milii.Callorhinchus_milii-6.1.3.113.gtf C_milii.gtf

#Clytia hemisphaerica
wget http://marimba.obs-vlfr.fr/sites/default/files/download/clytia_hm2.fasta && mv clytia_hm2.fasta C_hemisphaerica.fa

wget http://marimba.obs-vlfr.fr/sites/default/files/download/merged_transcript_models.gff3 && mv merged_transcript_models.gff3 C_hemisphaerica.gff3

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/728/285/GCF_902728285.1_Clytia_hemisphaerica_genome_assembly/GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.fna.gz && gunzip GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.fna.gz && mv GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.fna C_hemisphaerica.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/728/285/GCF_902728285.1_Clytia_hemisphaerica_genome_assembly/GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.gff.gz && gunzip GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.gff.gz && mv GCF_902728285.1_Clytia_hemisphaerica_genome_assembly_genomic.gff C_hemisphaerica.gff

#Daphnia pulex
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/daphnia_pulex_gca021134715v1rs/dna/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa.gz && gunzip Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa.gz && mv Daphnia_pulex_gca021134715v1rs.ASM2113471v1.dna.toplevel.fa D_pulex.fa

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/daphnia_pulex_gca021134715v1rs/Daphnia_pulex_gca021134715v1rs.ASM2113471v1.60.gff3.gz && gunzip Daphnia_pulex_gca021134715v1rs.ASM2113471v1.60.gff3.gz && mv Daphnia_pulex_gca021134715v1rs.ASM2113471v1.60.gff3 D_pulex.gff

#Drosophila melanogaster
wget https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz && gunzip Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz && mv Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa D_melanogaster.fa

wget https://ftp.ensembl.org/pub/release-113/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.113.gff3.gz && gunzip Drosophila_melanogaster.BDGP6.46.113.gff3.gz && mv Drosophila_melanogaster.BDGP6.46.113.gff3 D_melanogaster.gff3

#Exaiptasia diaphana
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz && gunzip GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz && mv GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna E_diaphana.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_genomic.gff.gz && gunzip GCF_001417965.1_Aiptasia_genome_1.1_genomic.gff.gz && mv GCF_001417965.1_Aiptasia_genome_1.1_genomic.gff E_diaphana.gff

#Gallus gallus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz && gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz && mv GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna G_gallus.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz && gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz && mv GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff G_gallus.gff

#Homo sapiens
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz && gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz && mv Homo_sapiens.GRCh38.dna.toplevel.fa H_sapiens.fa

wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz && gunzip Homo_sapiens.GRCh38.113.gff3.gz && mv Homo_sapiens.GRCh38.113.gff3 H_sapiens.gff3

#Hydra oligactis
wget https://research.nhgri.nih.gov/HydraAEP/download/sequences/h_oli/HOLI.genome.fa.gz && gunzip HOLI.genome.fa.gz && mv HOLI.genome.fa H_oligactis.fa && wget https://research.nhgri.nih.gov/HydraAEP/download/coordinates/h_oli/HOLI.GeneModels.gff3.gz && gunzip HOLI.GeneModels.gff3.gz && mv HOLI.GeneModels.gff3 H_oligactis.gff3

#Hydra viridissima
wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_genome_hm2_250116_renamed.fa.gz && gunzip hvir_genome_hm2_250116_renamed.fa.gz && mv hvir_genome_hm2_250116_renamed.fa H_viridissima.fa && wget https://marinegenomics.oist.jp/hydra_viridissima_a99/download/hvir_r06.all.recounted.gff.gz && gunzip hvir_r06.all.recounted.gff.gz && mv hvir_r06.all.recounted.gff H_viridissima.gff

#Hydra vulgaris 105
wget https://arusha.nhgri.nih.gov/hydra/download/assembly/Hm105_Dovetail_Assembly_1.0.fa.gz && gunzip Hm105_Dovetail_Assembly_1.0.fa.gz && mv Hm105_Dovetail_Assembly_1.0.fa H_vulgaris_105.fa && wget https://arusha.nhgri.nih.gov/hydra/download/genemodels_gff3/hydra2.0_genemodels.gff3.gz && gunzip hydra2.0_genemodels.gff3.gz && mv hydra2.0_genemodels.gff3 Hydra_vulgaris_105.gff3
cat ./Hydra_vulgaris_105.gff3 > H_vulgaris_105.gff3

#Hydractinia echinata
wget https://research.nhgri.nih.gov/hydractinia/download/assembly/echinata/Hech_primary_v1.0.fa.gz && gunzip Hech_primary_v1.0.fa.gz && mv Hech_primary_v1.0.fa H_echinata.fa && wget https://research.nhgri.nih.gov/hydractinia/download/genemodels_gff3/echinata/Hech_primary_v1.0.gff3.gz && gunzip Hech_primary_v1.0.gff3.gz && mv Hech_primary_v1.0.gff3 H_echinata.gff3

#Latimeria chalumnae
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_genomic.fna.gz && gunzip GCF_037176945.1_fLatCha1.pri_genomic.fna.gz && mv GCF_037176945.1_fLatCha1.pri_genomic.fna L_chalumnae.fa && wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_genomic.gff.gz && gunzip GCF_037176945.1_fLatCha1.pri_genomic.gff.gz && mv GCF_037176945.1_fLatCha1.pri_genomic.gff L_chalumnae.gff

#Lepisosteus oculatus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/954/835/GCF_040954835.1_fLepOcu1.hap2/GCF_040954835.1_fLepOcu1.hap2_genomic.fna.gz && gunzip GCF_040954835.1_fLepOcu1.hap2_genomic.fna.gz && mv GCF_040954835.1_fLepOcu1.hap2_genomic.fna L_oculatus.fa && wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/954/835/GCF_040954835.1_fLepOcu1.hap2/GCF_040954835.1_fLepOcu1.hap2_genomic.gff.gz && gunzip GCF_040954835.1_fLepOcu1.hap2_genomic.gff.gz && mv GCF_040954835.1_fLepOcu1.hap2_genomic.gff L_oculatus.gff

#Mnemiopsis leidyi
wget https://research.nhgri.nih.gov/mnemiopsis/download/genome/MlScaffold09.nt.gz && gunzip MlScaffold09.nt.gz && mv MlScaffold09.nt M_leidyi.fa && wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.gff3.gz && gunzip ML2.2.gff3.gz && mv ML2.2.gff3 M_leidyi.gff3

#Monosiga brevicollis
wget http://ftp.ensemblgenomes.org/pub/protists/release-60/fasta/protists_choanoflagellida1_collection/monosiga_brevicollis_mx1_gca_000002865/dna/Monosiga_brevicollis_mx1_gca_000002865.V1.0.dna.toplevel.fa.gz && gunzip Monosiga_brevicollis_mx1_gca_000002865.V1.0.dna.toplevel.fa.gz && mv Monosiga_brevicollis_mx1_gca_000002865.V1.0.dna.toplevel.fa M_brevicollis.fa && wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-60/gff3/protists_choanoflagellida1_collection/monosiga_brevicollis_mx1_gca_000002865/Monosiga_brevicollis_mx1_gca_000002865.V1.0.60.gff3.gz && gunzip Monosiga_brevicollis_mx1_gca_000002865.V1.0.60.gff3.gz && mv Monosiga_brevicollis_mx1_gca_000002865.V1.0.60.gff3 M_brevicollis.gff3

#Morbakka virulenta
wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_genome.fa.gz && gunzip MOR05_r06_genome.fa.gz && mv MOR05_r06_genome.fa M_virulenta.fa && wget https://marinegenomics.oist.jp/morbakka_virulenta/download/MOR05_r06_wref.gff3.gz && gunzip MOR05_r06_wref.gff3.gz && mv MOR05_r06_wref.gff3 M_virulenta.gff3

#Mus Musculus
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz && gunzip Mus_musculus.GRCm39.dna.toplevel.fa.gz && mv Mus_musculus.GRCm39.dna.toplevel.fa M_musculus.fa && wget https://ftp.ensembl.org/pub/release-113/gff3/mus_musculus/Mus_musculus.GRCm39.113.gff3.gz && gunzip Mus_musculus.GRCm39.113.gff3.gz && mv Mus_musculus.GRCm39.113.gff3 M_musculus.gff3

#Nematostella vectensis
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/Nvec200.fasta && mv Nvec200.fasta N_vectensis.fa && wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/aligned/tcs_v2/20240221/NV2g.20240221.gff && mv NV2g.20240221.gff N_vectensis.gff

#Patiria miniata
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/706/575/GCF_015706575.1_ASM1570657v1/GCF_015706575.1_ASM1570657v1_genomic.fna.gz && gunzip GCF_015706575.1_ASM1570657v1_genomic.fna.gz && mv GCF_015706575.1_ASM1570657v1_genomic.fna P_miniata.fa && wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/706/575/GCF_015706575.1_ASM1570657v1/GCF_015706575.1_ASM1570657v1_genomic.gff.gz && gunzip GCF_015706575.1_ASM1570657v1_genomic.gff.gz && mv GCF_015706575.1_ASM1570657v1_genomic.gff P_miniata.gff

#Porites lutea
wget http://plut.reefgenomics.org/download/plut_final_2.1.fasta.gz && gunzip plut_final_2.1.fasta.gz && mv plut_final_2.1.fasta P_lutea.fa && wget http://plut.reefgenomics.org/download/plut2v1.1.genes.gff3.gz && gunzip plut2v1.1.genes.gff3.gz && mv plut2v1.1.genes.gff3 P_lutea.gff3

#Rhopilema esculentum
wget https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100720/genome.fasta && mv genome.fasta R_esculentum.fa && wget https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100720/genome.gff && mv genome.gff R_esculentum.gff

#Scolanthus callimorphus
wget --no-check-certificate https://simrbase.stowers.org/files/pub/nematostella/Scal/genomes/Scal100/Scal100.fasta && mv Scal100.fasta S_callimorphus.fa && wget --no-check-certificate https://simrbase.stowers.org/files/pub/seaanemone/Scal/genomes/Scal100/aligned/NY_Scal100_v1/NY_Scal100_v1.20200813.gff && mv NY_Scal100_v1.20200813.gff S_callimorphus.gff

#Strongylocentrotus purpuratus
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/strongylocentrotus_purpuratus/dna/Strongylocentrotus_purpuratus.Spur_5.0.dna.toplevel.fa.gz && gunzip Strongylocentrotus_purpuratus.Spur_5.0.dna.toplevel.fa.gz && mv Strongylocentrotus_purpuratus.Spur_5.0.dna.toplevel.fa S_purpuratus.fa && wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/strongylocentrotus_purpuratus/Strongylocentrotus_purpuratus.Spur_5.0.60.gff3.gz && gunzip Strongylocentrotus_purpuratus.Spur_5.0.60.gff3.gz && mv Strongylocentrotus_purpuratus.Spur_5.0.60.gff3 S_purpuratus.gff3

#Trichoplax adhaerens
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-60/fasta/trichoplax_adhaerens/dna/Trichoplax_adhaerens.ASM15027v1.dna.toplevel.fa.gz && gunzip Trichoplax_adhaerens.ASM15027v1.dna.toplevel.fa.gz && mv Trichoplax_adhaerens.ASM15027v1.dna.toplevel.fa T_adhaerens.fa && wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/gff3/trichoplax_adhaerens/Trichoplax_adhaerens.ASM15027v1.60.gff3.gz && gunzip Trichoplax_adhaerens.ASM15027v1.60.gff3.gz && mv Trichoplax_adhaerens.ASM15027v1.60.gff3 T_adhaerens.gff3

#Xenopus tropicalis
wget https://ftp.ensembl.org/pub/release-113/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz && gunzip Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz && mv Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa X_tropicalis.fa && wget https://ftp.ensembl.org/pub/release-113/gff3/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.113.gff3.gz && gunzip Xenopus_tropicalis.UCB_Xtro_10.0.113.gff3.gz && mv Xenopus_tropicalis.UCB_Xtro_10.0.113.gff3 X_tropicalis.gff3

#Hydra vulgaris AEP 
wget https://research.nhgri.nih.gov/HydraAEP/download/sequences/hv_aep/HVAEP.genome.fa.gz && gunzip HVAEP.genome.fa.gz && mv HVAEP.genome.fa H_vulgaris_AEP.fa && wget https://research.nhgri.nih.gov/HydraAEP/download/coordinates/hv_aep/HVAEP.GeneModels.gff3.gz && gunzip HVAEP.GeneModels.gff3.gz && mv HVAEP.GeneModels.gff3 H_vulgaris_AEP.gff3

#Hydra vulgaris JNIG
cp /mnt/4e68f331-9abf-401e-87a7-745930b0c5ce/Hydra-hic/uc.fasta ./H_vulgaris_JNIG.fa
cp /mnt/4e68f331-9abf-401e-87a7-745930b0c5ce/Hydra-hic/braker.gtf ./H_vulgaris_JNIG.gtf

#calvadosia_cd-hit
wget https://github.com/josephryan/Ohdera_et_al_2018/blob/master/AA_Files/Ccrux.Trinity.faa.gz && gunzip Ccrux.Trinity.faa.gz 
mamba install -c bioconda -y cd-hit
cd-hit -i ./Ccrux.Trinity.faa -o C_cruxmelitensis.aa -c 0.95
#calvadosiaはinputから削除した
```

## 2. Generate proteomes from genome sequences
```bash
##最長アイソフォームからタンパク質配列を抽出
mamba install -c bioconda agat -y
for f in *.gtf ; do sed agat_convert_sp_gxf2gxf.pl --gff ./braker.gtf $f ; done
for f in *.gff ; do out="${f%.gff}.longestIso.gff"; agat_sp_keep_longest_isoform.pl --gff "$f" -o "$out"; done
for f in *.gff3; do out="${f%.gff3}.longestIso.gff"; agat_sp_keep_longest_isoform.pl --gff "$f" -o "$out"; done
for f in *.longestIso.gff; do base="${f%.longestIso.gff}"; gffread "$f" -g "${base}.fa" -y "${base}.prot.longestIso.fa"; done
##gffに文字化けがあったため
agat_sp_extract_sequences.pl --gff C_owczarzaki.longestIso.gff --fasta C_owczarzaki.fa -p --output C_owczarzaki.prot.longestIso.fa

##orthofinderでのエラー防止
mkdir ./prot.longestIso.fa
cp ./*.prot.longestIso.fa ./prot.longestIso.fa/
cd prot.longestIso.fa/

for f in *fa ; do sed -i 's/\(^>[^ \|]\+\).*/\1/g' $f ; done
for f in *.fa; do sed -i 's/\*//g' $f ; done
for f in *fa ; do sed -i '/^[^>]/ s/\.//g' "$f"; done
cd ..

#phylogenetic_tree.txt
((((((Hydractinia_echinata,Clytia_hemisphaerica),(Hydra_viridissima,(((Hydra_vulgaris_AEP,Hydra_vulgaris_JNIG),Hydra_vulgaris_105),Hydra_oligactis))),(Rhopilema_esculentum,Morbakka_virulenta)),((Acropora_millepora,Porites_lutea),(Exaiptasia_diaphana,(Nematostella_vectensis,Scolanthus_callimorphus)))),Trichoplax_adhaerens),((Adineta_vaga,((Daphnia_pulex,(Bombyx_mori,Drosophila_melanogaster)),Caenorhabditis_elegans)),((Patiria_miniata,Strongylocentrotus_purpuratus),(Callorhinchus_milii,((Latimeria_chalumnae,(Xenopus_tropicalis,((Homo_sapiens,Mus_musculus),Gallus_gallus))),Lepisosteus_oculatus))))),Amphimedon_queenslandica),Mnemiopsis_leidyi),Monosiga_brevicollis),Capsaspora_owczarzaki);

#/prot.longestIso.fa/以下のファイル名を以下に変更
Acropora_millepora.fa
Adineta_vaga.fa
Amphimedon_queenslandica.fa
Bombyx_mori.fa
Caenorhabditis_elegans.fa
Callorhinchus_milii.fa
Capsaspora_owczarzaki.fa
Clytia_hemisphaerica.fa
Daphnia_pulex.fa
Drosophila_melanogaster.fa
Exaiptasia_diaphana.fa
Gallus_gallus.fa
Homo_sapiens.fa
Hydra_oligactis.fa
Hydra_viridissima.fa
Hydra_vulgaris_105.fa
Hydra_vulgaris_AEP.fa
Hydra_vulgaris_JNIG.fa
Hydractinia_echinata.fa
Latimeria_chalumnae.fa
Lepisosteus_oculatus.fa
Mnemiopsis_leidyi.fa
Monosiga_brevicollis.fa
Morbakka_virulenta.fa
Mus_musculus.fa
Nematostella_vectensis.fa
Patiria_miniata.fa
Porites_lutea.fa
Rhopilema_esculentum.fa
Scolanthus_callimorphus.fa
Strongylocentrotus_purpuratus.fa
Trichoplax_adhaerens.fa
Xenopus_tropicalis.fa
```

## 3. Run OrthoFinder
```bash
conda install orthofinder -c bioconda
orthofinder -S diamond_ultra_sens -M msa -s ./phylogenetic_tree.txt -t 48 -f ./prot.longestIso.fa/
```

## 4. Run BLAST
```bash
#下記annotation.pyに必要なファイルを準備
cd /Orthofinder_wd

makeblastdb -in ./orthofinder/prot.longestIso.fa/Hydra_oligactis.fa -dbtype prot -out ./blast_db/H_oligactis_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_105.fa -dbtype prot -out ./blast_db/Hv_105_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_AEP.fa -dbtype prot -out ./blast_db/Hv_AEP_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_JNIG.fa -dbtype prot -out ./blast_db/Hv_JNIG_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -dbtype prot -out ./blast_db/Homo_sapiens_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Caenorhabditis_elegans.fa -dbtype prot -out ./blast_db/Caenorhabditis_elegans_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Drosophila_melanogaster.fa -dbtype prot -out ./blast_db/Drosophila_melanogaster_db
makeblastdb -in ./orthofinder/prot.longestIso.fa/Mus_musculus.fa -dbtype prot -out ./blast_db/Mus_musculus_db

#blastpを13通り分行う
#1 q:Ho s:Hv105
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_oligactis.fa -db ./blast_db/Hv_105_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_Holi_s_Hv105_blast.out
awk '$6 < 1e-5' ./q_Holi_s_Hv105_blast.out > ./q_Holi_s_Hv105_blast_significantHits.out

#2 q:Ho s:HvAEP
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_oligactis.fa -db ./blast_db/Hv_AEP_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_Holi_s_HvAEP_blast.out
awk '$6 < 1e-5' ./q_Holi_s_HvAEP_blast.out > ./q_Holi_s_HvAEP_blast_significantHits.out

#3 q:Ho s:HvJNIG
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_oligactis.fa -db ./blast_db/Hv_JNIG_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_Holi_s_HvJNIG_blast.out
awk '$6 < 1e-5' ./q_Holi_s_HvJNIG_blast.out > ./q_Holi_s_HvJNIG_blast_significantHits.out

#4 q:Hv105 s:Ho
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_105.fa -db ./blast_db/H_oligactis_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_Hv105_s_Holi_blast.out
awk '$6 < 1e-5' ./q_Hv105_s_Holi_blast.out > ./q_Hv105_s_Holi_blast_significantHits.out

#5 q:HvAEP s:Ho
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_AEP.fa -db ./blast_db/H_oligactis_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HvAEP_s_Holi_blast.out
awk '$6 < 1e-5' ./q_HvAEP_s_Holi_blast.out > ./q_HvAEP_s_Holi_blast_significantHits.out

#6 q:HvJNIG s:Ho
blastp -query ./orthofinder/prot.longestIso.fa/Hydra_vulgaris_JNIG.fa -db ./blast_db/H_oligactis_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HvJNIG_s_Holi_blast.out
awk '$6 < 1e-5' ./q_HvJNIG_s_Holi_blast.out > ./q_HvJNIG_s_Holi_blast_significantHits.out

#7 q:Hsapi s:HvJNIG
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Hv_JNIG_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Hv_JNIG_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Hv_JNIG_blast.out > ./q_HomoSapi_s_Hv_JNIG_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Hv_JNIG_blast_significantHits.out

#8 q:Hsapi s:HvAEP
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Hv_AEP_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Hv_AEP_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Hv_AEP_blast.out > ./q_HomoSapi_s_Hv_AEP_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Hv_AEP_blast_significantHits.out

#9 q:Hsapi s:Hv105
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Hv_105_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Hv_105_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Hv_105_blast.out > ./q_HomoSapi_s_Hv_105_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Hv_105_blast_significantHits.out

#10 q:Hsapi s:Ho
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/H_oligactis_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_H_oligactis_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_H_oligactis_blast.out > ./q_HomoSapi_s_H_oligactis_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_H_oligactis_blast_significantHits.out

#11 q:Hsapi s:Cle
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Caenorhabditis_elegans_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Cele_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Cele_blast.out > ./q_HomoSapi_s_Cele_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Cele_blast_significantHits.out

#12 q:Hsapi s:Dmelano
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Drosophila_melanogaster_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Dmelano_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Dmelano_blast.out > ./q_HomoSapi_s_Dmelano_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Dmelano_blast_significantHits.out

#13 q:Hsapi s:Mmus
blastp -query ./orthofinder/prot.longestIso.fa/Homo_sapiens.fa -db ./blast_db/Mus_musculus_db -num_threads 48 -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" -out ./q_HomoSapi_s_Mmus_blast.out
awk '$6 < 1e-5' ./q_HomoSapi_s_Mmus_blast.out > ./q_HomoSapi_s_Mmus_blast_significantHits.out
sed -i 's/^transcript://' q_HomoSapi_s_Mmus_blast_significantHits.out
```

## 5. Annotate orthogroups
```bash
python3 annotation.py
```

## 6. Perform enrichment analysis
```bash
#oligactisとvulgarisで保有に差が見られたオルソログのエンリッチメント解析
#annotation.tsvからHoli_o_Hv_x.txt, Holi_x_Hv_o.txtを作成する
cd ../Orthofinder_wd
Rscript enrichment.R
```


## 7. Summarize the number of orthogroups per species and per pathway
```python
#作成したannotation.tsvを使って集計する
#種ごとのGOやKEGGタームのOG数を数える
import pandas as pd
df = pd.read_csv("annotation.tsv", sep="\t")
cols = ["Orthogroup", "h_sapiens_gene", "database", "Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]
df = df[cols]

#メモ
#(df["h_sapiens_gene"] == "").sum()
#0→遺伝子のアノテーションが無い行は全てNaN

df["Hsapi_present"] = df["h_sapiens_gene"].notna().astype(int)

cols = ["Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]
df[cols] = (df[cols] != 0).astype(int)

tmp = df.copy()
tmp["database_term"] = tmp["database"].fillna("").astype(str).str.split(",")
tmp = tmp.explode("database_term")
tmp["database_term"] = tmp["database_term"].str.strip()
tmp = tmp[tmp["database_term"] != ""]

#tmp["database_term"].unique()から得た
terms=['GO:0031012_extracellular_matrix',
	'chronic_inflammation',
	'hsa04014_Ras_signaling_pathway',
	'hsa04152_AMPK_signaling_pathway',
	'hsa04340_Hedgehog_signaling_pathway',
    	'hsa04550_Signaling_pathways_regulating_pluripotency_of_stem_cells',
    	'SASP',
	'hsa04350_TGF-beta_signaling_pathway',
    	'hsa04151_PI3K-Akt_signaling_pathway',
    	'hsa04010_MAPK_signaling_pathway',
    	'hsa04120_Ubiquitin_mediated_proteolysis',
   	'GO:0031577_spindle_checkpoint_signaling',
    	'GO:0006338_chromatin_remodeling',
    	'GO:0040029_epigenetic_regulation_of_gene_expression',
    	'GO:0045814_negative_regulation_of_gene_expression',
    	'cellular_communication',
	'hsa04137_mitophagy',
    	'hsa04620_Toll-like_receptor_signaling_pathway',
    	'hsa04140_autophagy',
   	'GO:1900225_regulation_of_NLRP3_inflammasome_complex_assembly',
   	'GO:0030968_endoplasmic_reticulum_unfolded_protein_response',
   	'hsa03083_polycomb_repressive_complex',
    	'hsa04218_Cellular_senescence',
	'stem_cell_exhaustion',
    	'GO:0038084_vascular_endothelial_growth_factor_signaling_pathway',
   	'GO:0008286_insulin_receptor_signaling_pathway',
    	'hsa04142_lysosome',
    	'GO:0031445_regulation_of_heterochromatin_formation',
    	'GO:0141137_positive_regulation_of_gene_expression',
    	'Cellular_senescence',
	'hsa04310_Wnt_signaling_pathway',
    	'ETS_domain',
	'GO:0032200_telomere_organization',
   	'hsa04330_Notch_signaling_pathway',
	'GO:0005125_cytokine_activity',
    	'GO:0140896_cGAS_STING_signaling_pathway',
	'HSP70',
    	'GO:0046902_regulation_of_mitochondrial_membrane_permeability',
    	'GO:0005840_ribosome',
   	'GO:0140807_NAD_protein_glutamate_ADP_ribosyltransferase_activity',
    	'GO:0019985_translesion_synthesis',
    	'GO:0042054_histone_methyltransferase_activity',
    	'hsa04211_longevity_regulating_pathway',
	'disabled_macroautophagy',
   	'epigenetic_alteration',
	'genomic_instability',
    	'GO:0008459_chondroitin_6-sulfotransferase_activity',
    	'GO:0004402_histone_acetyltransferase_activity',
    	'hsa03410_base_excision_repair',
	'UPR_mt',
    	'GO:0004407_histone_deacetylase_activity',
    	'hsa03420_nucleotide_excision_repair',
    	'GO:0097472_cyclin_dependent_protein_kinase_activity',
    	'GO:0032452_histone_demethylase_activity',
    	'hsa03440_Homologous_recombination',
    	'hsa03450_Non-homologous_end-joining',
    	'deregelated_nutrient_sensing',
	'loss_of_proteostasis',
   	'mitochondrial_dysfunction',
	'hsa03430_Mismatch_repair',
    	'GO:0031931_TORC1_complex',
	'hsa03050_proteasome',
    	'GO:0009008_DNA-methyltransferase_activity',
	'dysbiosis']

presence_cols = ["Hsapi_present", "Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]


#各database_termの和を出す
result = (
    tmp[tmp["database_term"].isin(terms)]
    .groupby("database_term", as_index=False)[presence_cols]
    .sum()
)

result.to_csv("database_term_presence_sums.tsv", sep="\t", index=False)

#種ごとの老化関連遺伝子数を集計
df_sub = df[df["database"].notna()].copy()
cols = ["Hsapi_present", "Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]
species_sum = df_sub[cols].sum(axis=0)

#メモ：>>> species_sum#
Hsapi_present      2057
Mmus_present       2055
Dmelano_present    1665
Cele_present       1577
HvAEP_present      1724
HvJNIG_present     1725
Hv105_present      1720
Hv_present         1738
Ho_present         1745
```

## 8. Identify Ho/Hv-specific orthogroups
```python
import pandas as pd
df = pd.read_csv("annotation.tsv", sep="\t")
cols = ["Orthogroup", "h_sapiens_gene", "database", "Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]
df = df[cols]
df["Hsapi_present"] = df["h_sapiens_gene"].notna().astype(int)
cols = ["Mmus_present", "Dmelano_present", "Cele_present", "HvAEP_present", "HvJNIG_present", "Hv105_present", "Hv_present", "Ho_present"]
df[cols] = (df[cols] != 0).astype(int)

Hv_specific = df[ (df["Hsapi_present"] == 1) & (df["Hv_present"] == 1) & (df["Ho_present"] == 0) ] 
Ho_specific = df[ (df["Hsapi_present"] == 1) & (df["Hv_present"] == 0) & (df["Ho_present"] == 1) ]
Hv_specific.to_csv("Hv_specific.tsv", sep="\t", index=False)
Ho_specific.to_csv("Ho_specific.tsv", sep="\t", index=False)
```
