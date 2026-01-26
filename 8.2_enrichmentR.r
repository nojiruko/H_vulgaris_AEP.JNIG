#cd enrichment/
#conda activateしない
# ClusterProfiler, ReactomePA, org.Hs.eg.db, STRINGdb, RCy3 などのインストール
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
BiocManager::install("org.Hs.eg.db") # ヒト遺伝子アノテーションDB

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)

gene_raw <- readLines("Holi_o_Hv_x.txt")

gene_vec <- unlist(strsplit(gene_raw, ",\\s*"))  

gene.df <- bitr(gene_vec, 
                fromType = "SYMBOL",
                toType   = "ENTREZID", 
                OrgDb    = org.Hs.eg.db)
              
entrez_genes <- gene.df$ENTREZID

## GO BP (biological process) 解析
ego_bp <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",      # BP, MF, CC から選択
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

bp_plot <- barplot(ego_bp, showCategory = 10, title = "GO-BP Enrichment")
ggsave("Holi_o_Hv_x_GO_BP_Enrichment.pdf", plot = bp_plot, width = 8, height = 6)

ego_bp_df <- as.data.frame(ego_bp)
write.csv(ego_bp_df, "Holi_o_Hv_x_GO_BP_Enrichment.csv", row.names = FALSE)


## GO MF (molecular function) 解析
ego_mf <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
mf_plot <- barplot(ego_mf, showCategory = 10, title = "GO-MF Enrichment") +
  theme(axis.text.y = element_text(size = 8)) 
ggsave("Holi_o_Hv_x_GO_MF_Enrichment.pdf", plot = mf_plot, width = 8, height = 6)

ego_mf_df <- as.data.frame(ego_mf)
write.csv(ego_mf_df, "Holi_o_Hv_x_GO_MF_Enrichment.csv", row.names = FALSE)

## GO CC (cellular aomponent) 解析
ego_cc <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "CC",      # BP, MF, CC から選択
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
cc_plot <- barplot(ego_cc, showCategory = 10, title = "GO-CC Enrichment")
ggsave("Holi_o_Hv_GO_CC_Enrichment.pdf", plot = cc_plot, width = 8, height = 6)

ego_cc_df <- as.data.frame(ego_cc)
write.csv(ego_cc_df, "Holi_o_Hv_x_GO_CC_Enrichment.csv", row.names = FALSE)


# KEGG解析
ekegg <- enrichKEGG(gene         = entrez_genes,
                    organism     = "hsa",
                    keyType      = "kegg",
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

kegg_plot <- barplot(ekegg, showCategory = 10, title="KEGG Enrichment")
ggsave("Holi_o_Hv_x_KEGG_Enrichment.pdf", plot = kegg_plot, width = 8, height = 6)

ekegg_df <- as.data.frame(ekegg)
write.csv(ekegg_df, "Holi_o_Hv_x_KEGG_Enrichment.csv", row.names = FALSE)


# Reactome解析
ereact <- enrichPathway(gene          = entrez_genes,
                        organism      = "human",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

reactome_plot <- barplot(ereact, showCategory = 10, title="Reactome Pathway Enrichment")
ggsave("Holi_o_Hv_x_Reactome_Pathway_Enrichment.pdf", plot = reactome_plot, width = 8, height = 6)

ereact_df <- as.data.frame(ereact)
write.csv(ereact_df, "Holi_o_Hv_x_Reactome_Pathway_Enrichment.csv", row.names = FALSE)

#Holi_x_Hv_o
gene_raw <- readLines("Holi_x_Hv_o.txt")

gene_vec <- unlist(strsplit(gene_raw, ",\\s*"))  

gene.df <- bitr(gene_vec, 
                fromType = "SYMBOL",
                toType   = "ENTREZID", 
                OrgDb    = org.Hs.eg.db)
              
entrez_genes <- gene.df$ENTREZID

## GO BP (biological process) 解析
ego_bp <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",      # BP, MF, CC から選択
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

bp_plot <- barplot(ego_bp, showCategory = 10, title = "GO-BP Enrichment")
ggsave("Holi_x_Hv_o_GO_BP_Enrichment.pdf", plot = bp_plot, width = 8, height = 6)

ego_bp_df <- as.data.frame(ego_bp)
write.csv(ego_bp_df, "Holi_x_Hv_o_GO_BP_Enrichment.csv", row.names = FALSE)


## GO MF (molecular function) 解析
ego_mf <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
mf_plot <- barplot(ego_mf, showCategory = 10, title = "GO-MF Enrichment") +
  theme(axis.text.y = element_text(size = 8)) 
ggsave("Holi_x_Hv_o_GO_MF_Enrichment.pdf", plot = mf_plot, width = 8, height = 6)

ego_mf_df <- as.data.frame(ego_mf)
write.csv(ego_mf_df, "Holi_x_Hv_o_GO_MF_Enrichment.csv", row.names = FALSE)

## GO CC (cellular aomponent) 解析
ego_cc <- enrichGO(gene          = entrez_genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "CC",      # BP, MF, CC から選択
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
cc_plot <- barplot(ego_cc, showCategory = 10, title = "GO-CC Enrichment")
ggsave("Holi_x_Hv_o_GO_CC_Enrichment.pdf", plot = cc_plot, width = 8, height = 6)

ego_cc_df <- as.data.frame(ego_cc)
write.csv(ego_cc_df, "Holi_x_Hv_o_GO_CC_Enrichment.csv", row.names = FALSE)


# KEGG解析
ekegg <- enrichKEGG(gene         = entrez_genes,
                    organism     = "hsa",
                    keyType      = "kegg",
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

kegg_plot <- barplot(ekegg, showCategory = 10, title="KEGG Enrichment")
ggsave("Holi_x_Hv_o_KEGG_Enrichment.pdf", plot = kegg_plot, width = 8, height = 6)

ekegg_df <- as.data.frame(ekegg)
write.csv(ekegg_df, "Holi_x_Hv_o_KEGG_Enrichment.csv", row.names = FALSE)


# Reactome解析
ereact <- enrichPathway(gene          = entrez_genes,
                        organism      = "human",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

reactome_plot <- barplot(ereact, showCategory = 10, title="Reactome Pathway Enrichment")
ggsave("Holi_x_Hv_o_Reactome_Pathway_Enrichment.pdf", plot = reactome_plot, width = 8, height = 6)

ereact_df <- as.data.frame(ereact)
write.csv(ereact_df, "Holi_x_Hv_o_Reactome_Pathway_Enrichment.csv", row.names = FALSE)
