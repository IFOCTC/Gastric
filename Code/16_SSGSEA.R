## Code useful to perform:
## Single Sample Geneset Enrichment Analysis (ssGSEA)

## Pre-filtering 
columns_to_mantain_corr_gsea <- c("ID", "OS_time", "OS_Event", "PFS_time", "PFS_event",
                                  "Coorte", "CT1L_1SI", "CT1L_FLOT_1SI", "CD44", "LGALS1",
                                  "CTSZ", "PDLIM1", "GNAI2", "CD68", "IFITM3", "SH3BGRL3",
                                  "RSU1", "RAB8B", "TAF10", "SMAD3", "Condition_Q1",
                                  "Condition_Q3", "Condition_Mean")
genes_col <- c("CD44", "LGALS1", "CTSZ", "PDLIM1",
               "GNAI2", "CD68", "IFITM3", "SH3BGRL3",
               "RSU1", "RAB8B", "TAF10", "SMAD3")
df_gsea <- df_final_survival_signature_corr %>% 
  dplyr::select(all_of(columns_to_mantain_corr_gsea))

## Load Useful
library(GSEABase)
library(GSVA)
library(limma)
library(GSEABase)
library(edgeR)

## Hallmark
hallmark_pathways <- getGmt(paste0(data_path, "/h.all.v2023.2.Hs.symbols.gmt"))
expression_data <- as.matrix(df_gsea[,genes_col])
rownames(expression_data) <- df_gsea$ID
## Create the parameter object specific for ssGSEA
ssgsea_params <- ssgseaParam(t(expression_data), hallmark_pathways)
ssGSEA_scores <- gsva(ssgsea_params)

## Scaling
ssGSEA_scores <- as.data.frame(t(ssGSEA_scores))
ssGSEA_scores_scale <- scale(ssGSEA_scores)

## Heatmap
pheatmap(t(ssGSEA_scores), show_colnames = FALSE, main = "ssGSEA not scaled")
pheatmap(t(ssGSEA_scores_scale), show_colnames = FALSE, main = "ssGSEA scaled")

## Cut
info <- data.frame(Condition = ifelse(rowMeans(ssGSEA_scores) >= 
                                        median(rowMeans(ssGSEA_scores)), "High", "Low"))
rownames(info) <- rownames(ssGSEA_scores)
## Check
table(info)

## Plot
pheatmap(t(df_gsea[,genes_col]), show_colnames = FALSE, show_rownames = TRUE,
         main = "ssGSEA not scaled", annotation_col = info)
pheatmap(t(ssGSEA_scores_scale), show_colnames = FALSE, show_rownames = TRUE,
         main = "ssGSEA scaled", annotation_col = info)

