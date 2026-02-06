# This Python script adds annotations to the OrthoFinder output.

## Overview
#1. Name genes based on protein IDs
#2. Merge blast_significantHits files 

#===== 1. Name genes based on protein IDs =====
# Download BioMart data for each organism.
# For human and mouse, download using protein stable ID version.
# For fly and worm, download using protein stable ID.
# Download a table of NCBI Gene (formerly Entrez Gene) IDs and gene names.

import pandas as pd
orthogroup = pd.read_csv('./Orthogroups.tsv', sep='\t', dtype=str)

H_sapiens = pd.read_csv('./mart_h_sapiens.txt', sep='\t', dtype=str, header=0)
M_musculus = pd.read_csv('./mart_m_musculus.txt', sep='\t', dtype=str, header=0)
D_melanogaster = pd.read_csv('./mart_d_melanogaster.txt', sep='\t', dtype=str, header=0)
C_elegans = pd.read_csv('./mart_c_elegans.txt', sep='\t', dtype=str, header=0)

# Check for NaN values
orthogroup.isna().sum()
H_sapiens.isna().sum()

# Check for empty strings
(orthogroup == "").sum() #0
(H_sapiens == "").sum() #0

# Remove rows with missing values from mart.txt
H_sapiens = H_sapiens[H_sapiens['Transcript stable ID'].notna()]
H_sapiens = H_sapiens[H_sapiens['Gene name'].notna()]
M_musculus = M_musculus[M_musculus['Transcript stable ID'].notna()]
M_musculus = M_musculus[M_musculus['Gene name'].notna()]
D_melanogaster = D_melanogaster[D_melanogaster['Transcript stable ID'].notna()]
D_melanogaster = D_melanogaster[D_melanogaster['Gene name'].notna()]
C_elegans = C_elegans[C_elegans['Transcript stable ID'].notna()]
C_elegans = C_elegans[C_elegans['Gene name'].notna()]

# Remove duplicates to make protein IDs unique
H_sapiens = H_sapiens.drop_duplicates(subset=['Transcript stable ID', 'Gene name'])
H_sapiens = H_sapiens.groupby('Transcript stable ID').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()
M_musculus = M_musculus.drop_duplicates(subset=['Transcript stable ID', 'Gene name'])
M_musculus = M_musculus.groupby('Transcript stable ID').agg(lambda x: ','.join(x.dropna().unique())).reset_index()
D_melanogaster = D_melanogaster.drop_duplicates(subset=['Transcript stable ID', 'Gene name'])
D_melanogaster = D_melanogaster.groupby('Transcript stable ID').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()
C_elegans = C_elegans.drop_duplicates(subset=['Transcript stable ID', 'Gene name'])
C_elegans = C_elegans.groupby('Transcript stable ID').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()

orthogroup = orthogroup[['Orthogroup', 'Homo_sapiens', 'Mus_musculus', 'Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Hydra_oligactis', 'Hydra_vulgaris_105', 'Hydra_vulgaris_AEP', 'Hydra_vulgaris_JNIG']]
orthogroup = orthogroup.assign(Homo_sapiens=orthogroup["Homo_sapiens"].str.split(",")).explode("Homo_sapiens")
orthogroup["Homo_sapiens"] = orthogroup["Homo_sapiens"].str.strip()

annotation = pd.merge(orthogroup, H_sapiens, left_on='Homo_sapiens', right_on='Transcript stable ID', how='left')
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index() #NaNだった部分は''に置換され、annotation内はNaNが０、''だけになる

annotation = annotation.rename(columns={'Transcript stable ID': 'h_sapiens_id', 'Gene name': 'h_sapiens_gene'})

annotation = annotation.assign(Mus_musculus=annotation["Mus_musculus"].str.split(",")).explode("Mus_musculus")
annotation["Mus_musculus"] = annotation["Mus_musculus"].str.strip()
annotation = pd.merge(annotation, M_musculus, left_on='Mus_musculus', right_on='Transcript stable ID', how='left')
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().astype(str).unique())).reset_index()

annotation = annotation.rename(columns={'Transcript stable ID': 'm_musculus_id', 'Gene name': 'm_musculus_gene'})

annotation = annotation.assign(Drosophila_melanogaster=annotation["Drosophila_melanogaster"].str.split(",")).explode("Drosophila_melanogaster")
annotation["Drosophila_melanogaster"] = annotation["Drosophila_melanogaster"].str.strip()
annotation = pd.merge(annotation, D_melanogaster, left_on='Drosophila_melanogaster', right_on='Transcript stable ID', how='left')
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()

annotation = annotation.rename(columns={'Transcript stable ID': 'd_melanogaster_id', 'Gene name': 'd_melanogaster_gene'})

annotation = annotation.assign(Caenorhabditis_elegans=annotation["Caenorhabditis_elegans"].str.split(",")).explode("Caenorhabditis_elegans")
annotation["Caenorhabditis_elegans"] = annotation["Caenorhabditis_elegans"].str.strip()
annotation = pd.merge(annotation, C_elegans, left_on='Caenorhabditis_elegans', right_on='Transcript stable ID', how='left')
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()

annotation = annotation.rename(columns={'Transcript stable ID': 'c_elegans_id', 'Gene name': 'c_elegans_gene'})

annotation = annotation[['Orthogroup', 'Homo_sapiens', 'h_sapiens_gene', 'Mus_musculus', 'm_musculus_gene', 'Drosophila_melanogaster', 'd_melanogaster_gene', 'Caenorhabditis_elegans', 'c_elegans_gene', 'Hydra_oligactis', 'Hydra_vulgaris_105', 'Hydra_vulgaris_AEP', 'Hydra_vulgaris_JNIG']]

annotation.to_csv('./annotation_tmp1.csv', sep='\t', index=False)


# Annotation of aging-related genes
# Gene list aging_gene_lists_in_this_study.tsv was created based on the Hallmarks of Aging
aging_gene_lists = pd.read_csv('./aging_gene_lists_in_this_study.tsv', sep='\t', dtype=str, header=None)
annotation = annotation.assign(h_sapiens_gene=annotation["h_sapiens_gene"].str.split(",")).explode("h_sapiens_gene")
annotation["h_sapiens_gene"] = annotation["h_sapiens_gene"].str.strip()
annotation = pd.merge(annotation, aging_gene_lists, left_on='h_sapiens_gene', right_on=0, how='left') 
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()
annotation = annotation.rename(columns={0:'aging_gene_in_this_study', 1:'database', 2:'hallmarks'})
annotation = annotation[['Orthogroup', 'Homo_sapiens', 'h_sapiens_gene', 'Mus_musculus', 'm_musculus_gene', 'Drosophila_melanogaster', 'd_melanogaster_gene', 'Caenorhabditis_elegans', 'c_elegans_gene', 'Hydra_oligactis', 'Hydra_vulgaris_105', 'Hydra_vulgaris_AEP', 'Hydra_vulgaris_JNIG', 'aging_gene_in_this_study', 'database', 'hallmarks']]

# Add a column to check whether the same gene name is assigned to multiple orthogroups
annotation = annotation.assign(h_sapiens_gene=annotation["h_sapiens_gene"].str.split(',')).explode("h_sapiens_gene")
annotation["h_sapiens_gene"] = annotation["h_sapiens_gene"].str.strip()

annotation['h_sapiens_gene_is_duplicated'] = annotation['h_sapiens_gene'].duplicated(keep=False)

annotation = annotation.astype(str)
annotation = annotation.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()

annotation.to_csv('./annotation_tmp2.csv', sep='\t', index=False)


#===== 2. Merge blast_significantHits files =====
# Repeat explode and merge operations
## Reduce blast_significantHits.out to the required columns and make query IDs unique

def merge_blast_hits(
    annotation_df,
    species_col,          # ex) 'Hydra_oligactis'
    query_label,          # ex) 'Holi'
    seq_label,            # ex) 'Hv105'
    file_path,   # ex) './q_Holi_s_Hv105_blast_significantHits.out'
    added_colname        # ex) 'qHoli_sHv105_BlastSignificantHits'
):
    # Load and format BLAST files
    blast_df = (
        pd.read_csv(file_path, sep='\t', dtype=str, header=None)[[0, 1]]
        .rename(columns={0: query_label, 1: seq_label})
        .drop_duplicates(subset=[query_label, seq_label])
        .groupby(query_label, as_index=False)[seq_label]
        .agg(lambda x: ', '.join(x.dropna().unique()))
        .rename(columns={seq_label: added_colname})
    )

    # Explode and strip the target columns
    df = annotation_df
    df[species_col] = df[species_col].str.split(',')
    df = df.explode(species_col)
    df[species_col] = df[species_col].str.strip()

    # Merg
    df = pd.merge(df, blast_df, left_on=species_col, right_on=query_label, how='left').drop(columns=[query_label])

    # Group by orthogroup
    return df.groupby('Orthogroup').agg(lambda x: ', '.join(x.dropna().unique())).reset_index()

#for 13 combinations
#1
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_oligactis',
    query_label='Holi',
    seq_label='Hv105',
    file_path='./q_Holi_s_Hv105_blast_significantHits.out',
    added_colname='qHoli_sHv105_BlastSignificantHits'
)

#2
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_oligactis',
    query_label='Holi',
    seq_label='HvAEP',
    file_path='./q_Holi_s_HvAEP_blast_significantHits.out',
    added_colname='qHoli_sHvAEP_BlastSignificantHits'
)

#3
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_oligactis',
    query_label='Holi',
    seq_label='HvJNIG',
    file_path='./q_Holi_s_HvJNIG_blast_significantHits.out',
    added_colname='qHoli_sHvJNIG_BlastSignificantHits'
)

#4
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_vulgaris_105',
    query_label='Hv105',
    seq_label='Holi',
    file_path='./q_Hv105_s_Holi_blast_significantHits.out',
    added_colname='qHv105_sHoli_BlastSignificantHits'
)

#5
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_vulgaris_AEP',
    query_label='HvAEP',
    seq_label='Holi',
    file_path='./q_HvAEP_s_Holi_blast_significantHits.out',
    added_colname='qHvAEP_sHoli_BlastSignificantHits'
)

#6
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Hydra_vulgaris_JNIG',
    query_label='HvJNIG',
    seq_label='Holi',
    file_path='./q_HvJNIG_s_Holi_blast_significantHits.out',
    added_colname='qHvJNIG_sHoli_BlastSignificantHits'
)

#7
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='HvJNIG',
    file_path='./q_HomoSapi_s_Hv_JNIG_blast_significantHits.out',
    added_colname='qHomoSapi_sHvJNIG_BlastSignificantHits'
)

#8
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='HvAEP',
    file_path='./q_HomoSapi_s_Hv_AEP_blast_significantHits.out',
    added_colname='qHomoSapi_sHvAEP_BlastSignificantHits'
)

#9
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='Hv105',
    file_path='./q_HomoSapi_s_Hv_105_blast_significantHits.out',
    added_colname='qHomoSapi_sHv105_BlastSignificantHits'
)

#10
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='Ho',
    file_path='./q_HomoSapi_s_H_oligactis_blast_significantHits.out',
    added_colname='qHomoSapi_sHo_BlastSignificantHits'
)

#11
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='Cele',
    file_path='./q_HomoSapi_s_Cele_blast_significantHits.out',
    added_colname='qHomoSapi_sCele_BlastSignificantHits'
)

#12
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='Dmelano',
    file_path='./q_HomoSapi_s_Dmelano_blast_significantHits.out',
    added_colname='qHomoSapi_sDmelano_BlastSignificantHits'
)

#13
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='Homo_sapiens',
    query_label='HomoSapi',
    seq_label='Mmus',
    file_path='./q_HomoSapi_s_Mmus_blast_significantHits.out',
    added_colname='qHomoSapi_sMmus_BlastSignificantHits'
)

#14
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHvJNIG_BlastSignificantHits',
    query_label='HvJNIG',
    seq_label='Holi',
    file_path='./q_HvJNIG_s_Holi_blast_significantHits.out',
    added_colname='qHsHvJNIG_sHoli_BlastSignificantHits'
)

#15
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHvAEP_BlastSignificantHits',
    query_label='HvAEP',
    seq_label='Holi',
    file_path='./q_HvAEP_s_Holi_blast_significantHits.out',
    added_colname='qHsHvAEP_sHoli_BlastSignificantHits'
)

#16
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHv105_BlastSignificantHits',
    query_label='Hv105',
    seq_label='Holi',
    file_path='./q_Hv105_s_Holi_blast_significantHits.out',
    added_colname='qHsHv105_sHoli_BlastSignificantHits'
)

#17
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHo_BlastSignificantHits',
    query_label='Holi',
    seq_label='HvJNIG',
    file_path='./q_Holi_s_HvJNIG_blast_significantHits.out',
    added_colname='qHsHoli_sHvJNIG_BlastSignificantHits'
)

#18
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHo_BlastSignificantHits',
    query_label='Holi',
    seq_label='HvAEP',
    file_path='./q_Holi_s_HvAEP_blast_significantHits.out',
    added_colname='qHsHoli_sHvAEP_BlastSignificantHits'
)

#19
annotation = merge_blast_hits(
    annotation_df=annotation,
    species_col='qHomoSapi_sHo_BlastSignificantHits',
    query_label='Holi',
    seq_label='Hv105',
    file_path='./q_Holi_s_Hv105_blast_significantHits.out',
    added_colname='qHsHoli_sHv105_BlastSignificantHits'
)

# Create _binary columns for each column
blast_cols = [col for col in annotation.columns if col.endswith('_BlastSignificantHits')]
for col in blast_cols:
    new_col = f'{col}_binary'
    annotation[new_col] = annotation[col].apply(lambda x: 1 if pd.notna(x) and x.strip() != '' else 0)

species_cols = ['Homo_sapiens', 'Mus_musculus', 'Drosophila_melanogaster', 'Caenorhabditis_elegans', 'Hydra_oligactis', 'Hydra_vulgaris_105', 'Hydra_vulgaris_AEP', 'Hydra_vulgaris_JNIG']

for col in species_cols:
    new_col = f'{col}_binary'
    annotation[new_col] = annotation[col].apply(lambda x: 1 if pd.notna(x) and str(x).strip() != '' else 0)

binary_cols = [col for col in annotation.columns if col.endswith('_binary')]

annotation[binary_cols] = annotation[binary_cols].astype(int)

# Summarize whether sequences from each species are present in each orthogroup
# Check the presence of Mmus
annotation["Mmus_present"] = annotation["Mus_musculus_binary"] + annotation["qHomoSapi_sMmus_BlastSignificantHits_binary"]

# Check the presence of Dmelano
annotation["Dmelano_present"] = annotation["Drosophila_melanogaster_binary"] + annotation["qHomoSapi_sDmelano_BlastSignificantHits_binary"]
        
# Check the presence of Cele
annotation["Cele_present"] = annotation["Caenorhabditis_elegans_binary"] + annotation["qHomoSapi_sCele_BlastSignificantHits_binary"]
        
# Check the presence of HvAEP
annotation["HvAEP_present"] = annotation["Hydra_vulgaris_AEP_binary"] + annotation["qHomoSapi_sHvAEP_BlastSignificantHits_binary"] + annotation["qHoli_sHvAEP_BlastSignificantHits_binary"] + annotation["qHsHoli_sHvAEP_BlastSignificantHits_binary"]

# Check the presence of HvJNIG
annotation["HvJNIG_present"] = annotation["Hydra_vulgaris_JNIG_binary"] + annotation["qHomoSapi_sHvJNIG_BlastSignificantHits_binary"] + annotation["qHoli_sHvJNIG_BlastSignificantHits_binary"] + annotation["qHsHoli_sHvJNIG_BlastSignificantHits_binary"]

# Check the presence of Hv105
annotation["Hv105_present"] = annotation["Hydra_vulgaris_105_binary"] + annotation["qHomoSapi_sHv105_BlastSignificantHits_binary"] + annotation["qHoli_sHv105_BlastSignificantHits_binary"] + annotation["qHsHoli_sHv105_BlastSignificantHits_binary"]

# Check the presence of Hv
annotation["Hv_present"] = annotation["Hv105_present"] + annotation["HvAEP_present"] + annotation["HvJNIG_present"] 

# Check the presence of Ho
annotation["Ho_present"] = annotation["Hydra_oligactis_binary"] + annotation["qHomoSapi_sHo_BlastSignificantHits_binary"] + annotation["qHv105_sHoli_BlastSignificantHits_binary"] + annotation["qHvAEP_sHoli_BlastSignificantHits_binary"] + annotation["qHvJNIG_sHoli_BlastSignificantHits_binary"] + annotation["qHsHvAEP_sHoli_BlastSignificantHits_binary"] + annotation["qHsHvJNIG_sHoli_BlastSignificantHits_binary"] + annotation["qHsHv105_sHoli_BlastSignificantHits_binary"]

annotation.to_csv('annotation.tsv', sep='\t', index=False)


