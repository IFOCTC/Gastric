## Code useful to perform wrangling operations

rm(list = ls())
## Load scripts
source("00_Support.R")

## Load data
df_tpm_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.star_tpm.tsv"), show_col_types = FALSE)
df_clinical_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
df_survival_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.survival.tsv"), show_col_types = FALSE)
df_tpm_own  <- read_tsv(paste0(data_path_own, "/salmon.merged.gene_tpm.tsv"), show_col_types = FALSE)

## Select only common genes 
top_genes_common <- c("CD44", "CTSZ", "PDLIM1", "GNAI2", "IFITM3",
                      "RSU1", "RAB8B", "TAF10", "SMAD3")
## DEGs = 0.8, p-val = 0.1 (30 genes)
top_genes <- c("CD44", "PFKP", "FAM50A", "CD82", "LGALS1",
               "HIF1A", "CTSZ", "CTSH", "HSPB1", "PDLIM1",
               "OAS2", "DBN1", "GNAI2", "ECE1", "CTSD",
               "MORF4L2", "CD68", "PHF11", "IFITM3", "NONO",
               "RSU1", "CDC123", "MCU", "INPPL1", "RAB8B",
               "TAF10", "SMAD3", "SART1", "FAM3C", "PSMB9")
## DEGs = 0.8, p-val = 0.05 (16 genes)
top_genes_02 <- c("CD44", "PFKP", "FAM50A", "CD82", "CTSZ",
                  "CTSH", "HSPB1", "PDLIM1", "DBN1", "GNAI2",
                  "IFITM3", "RSU1", "CDC123", "RAB8B", "TAF10", "SMAD3")
## Top genes corr
top_genes_corr <- c("CD44", "LGALS1", "CTSZ", "PDLIM1", "GNAI2",
                    "CD68", "IFITM3", "RSU1", "RAB8B", "TAF10",
                    "SMAD3")
## Top genes corr 02
top_genes_corr_02 <- c("CD44", "CTSZ", "GNAI2", "RSU1", "RAB8B", "TAF10", "SMAD3")

## Filtering
top_genes_id <- df_tpm_own %>% 
  dplyr::select(gene_id, gene_name) %>% 
  filter(gene_name %in% top_genes_02)
## Wrangling ensembl id gene name
Ensembl_ID <- df_tpm_tcga$Ensembl_ID
Ensembl_ID_cleaned <- sub("\\..*", "", Ensembl_ID) 
df_tpm_tcga_cleaned <- df_tpm_tcga
df_tpm_tcga_cleaned$Ensembl_ID <- Ensembl_ID_cleaned
## Check
top_genes_id$gene_id %in% Ensembl_ID_cleaned

## Filtering
df_tpm_tcga_cleaned <- df_tpm_tcga_cleaned %>% 
  filter(Ensembl_ID %in% top_genes_id$gene_id)
## Add gene name
df_tpm_tcga_cleaned$Ensembl_ID <- NULL
rownames(df_tpm_tcga_cleaned) <- top_genes_id$gene_name
## Transpose
df_tpm_tcga_cleaned <- t(df_tpm_tcga_cleaned)
df_tpm_tcga_cleaned <- data.frame(df_tpm_tcga_cleaned)

## TPM
## Select only tumor samples [01-09]
selected_rownames <- rownames(df_tpm_tcga_cleaned)[grep("01A$", rownames(df_tpm_tcga_cleaned))]
## Filtering
df_tpm_tcga_cleaned_tumor <- df_tpm_tcga_cleaned %>% 
  filter(rownames(df_tpm_tcga_cleaned) %in% selected_rownames)
## Remove the last part of barcode
selected_rownames_cleaned <- sub("-[^-]+$", "", rownames(df_tpm_tcga_cleaned_tumor))
rownames(df_tpm_tcga_cleaned_tumor) <- selected_rownames_cleaned

## CLINICAL
## Select only tumor samples [01-09]
rownames(df_clinical_tcga) <- df_clinical_tcga$sample
selected_rownames_clinical <- rownames(df_clinical_tcga)[grep("01A$", rownames(df_clinical_tcga))]
## Filtering
df_clinical_tcga_cleaned_tumor <- df_clinical_tcga %>% 
  filter(rownames(df_clinical_tcga) %in% selected_rownames_clinical)
## Remove the last part of barcode
selected_rownames_clinical_cleaned <- sub("-[^-]+$", "", rownames(df_clinical_tcga_cleaned_tumor))
rownames(df_clinical_tcga_cleaned_tumor) <- selected_rownames_clinical_cleaned

## SURVIVAL
## Select only tumor samples [01-09]
rownames(df_survival_tcga) <- df_survival_tcga$sample
selected_rownames_survival <- rownames(df_survival_tcga)[grep("01A$", rownames(df_survival_tcga))]
## Filtering
df_survival_tcga_cleaned_tumor <- df_survival_tcga %>% 
  filter(rownames(df_survival_tcga) %in% selected_rownames_survival)
## Remove the last part of barcode
selected_rownames_survival_cleaned <- sub("-[^-]+$", "", rownames(df_survival_tcga_cleaned_tumor))
rownames(df_survival_tcga_cleaned_tumor) <- selected_rownames_survival_cleaned

## Merging
length(intersect(intersect(rownames(df_survival_tcga_cleaned_tumor),
                           rownames(df_tpm_tcga_cleaned_tumor)),
                 rownames(df_clinical_tcga_cleaned_tumor)))

## Useulf operation for merge in the next script
df_survival_tcga_cleaned_tumor <-  df_survival_tcga_cleaned_tumor %>%
  rename(ID = `_PATIENT`)
df_tpm_tcga_cleaned_tumor$ID      <- rownames(df_tpm_tcga_cleaned_tumor)
df_clinical_tcga_cleaned_tumor$ID <- rownames(df_clinical_tcga_cleaned_tumor)

## Save tsv files
write.table(df_survival_tcga_cleaned_tumor,
            file = paste0(output_path, "/survival_cleaned.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df_tpm_tcga_cleaned_tumor, 
            file = paste0(output_path, "/tpm_cleaned.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df_clinical_tcga_cleaned_tumor, 
            file = paste0(output_path, "/clinical_cleaned.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
