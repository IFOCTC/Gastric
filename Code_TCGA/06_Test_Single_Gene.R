source("00_Support.R")

## Survival curve for each gene
top_genes <- c("CD44", "PFKP", "FAM50A", "CD82", "LGALS1",
               "HIF1A", "CTSZ", "CTSH", "HSPB1", "PDLIM1",
               "OAS2", "DBN1", "GNAI2", "ECE1", "CTSD",
               "MORF4L2", "CD68", "PHF11", "IFITM3", "NONO",
               "RSU1", "CDC123", "MCU", "INPPL1", "RAB8B",
               "TAF10", "SMAD3", "SART1", "FAM3C", "PSMB9")
## Load data
df_tpm_own <- read_tsv(paste0(data_path_own, "/salmon.merged.gene_tpm.tsv"),
                       show_col_types = FALSE)
## Load survival data
df_survival_own <- read_excel(paste0(data_path_own, "/stomaci_wes_all.xlsx"))

## Filter 
df_tpm_own_filtered <- df_tpm_own %>% 
  filter(df_tpm_own$gene_name %in% top_genes)
df_tpm_own_filtered <- df_tpm_own_filtered %>% 
  dplyr::select(matches("T$"))
patients   <- colnames(df_tpm_own_filtered)
genes_name <- top_genes
df_tpm_own_filtered <- data.frame(t(df_tpm_own_filtered))
rownames(df_tpm_own_filtered) <- patients
colnames(df_tpm_own_filtered) <- genes_name
df_tpm_own_filtered$ID <- rownames(df_tpm_own_filtered)
df_tpm_own_filtered$ID <- gsub("\\.1T$", "", df_tpm_own_filtered$ID)
df_tpm_own_filtered$ID <- gsub("\\.", "-", df_tpm_own_filtered$ID)
## Log-Transformation
exclude_cols <- c("ID")
df_support <- df_tpm_own_filtered %>% 
  dplyr::select(-all_of(exclude_cols)) %>%
  mutate(across(everything(), ~ log2(. + 1)))
df_tpm_own_filtered <- df_tpm_own_filtered %>%
  dplyr::select(all_of(exclude_cols)) %>%
  bind_cols(df_support)
df_tpm_own_filtered_merged <- merge(df_tpm_own_filtered,
                                    df_survival_own, by = "ID")
df_tpm_own_filtered_merged <- df_tpm_own_filtered_merged %>% 
  filter(CT1L_FLOT_1SI == 1) ## 155
df_tpm_own_filtered_merged$ID <- NULL

## Survival analysis
max_months <- 48
output_folder_os <- "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Res/Images/Survival_Test_Single_Gene/OS"
## OS
g <- colnames(df_tpm_own_filtered_merged[,1:30])  ## Get the column names
for(i in seq_along(g)){  ## Iterate over the column names by index
  gene_col <- g[i]  ## Get the current gene column name
  ## Create Group based on the median of the Mean column
  df_tpm_own_filtered_merged$Group <- cut(df_tpm_own_filtered_merged[[gene_col]],
                                          breaks = quantile(df_tpm_own_filtered_merged[[gene_col]],
                                                            probs = c(0, 1/3, 2/3, 1)),
                                          labels = c("Low", "Medium", "High"),
                                          include.lowest = TRUE)
  ## Perform survival analysis
  fit <- survfit(Surv(OS_time, OS_Event) ~ Group, data = df_tpm_own_filtered_merged)
  pairwise_pvalues_os_own <- pairwise_survdiff(Surv(OS_time, as.integer(OS_Event)) ~ Group,
                                               data = df_tpm_own_filtered_merged,
                                               p.adjust.method = "none")
  pairwise_pvalues_os_own_pvals <- pairwise_pvalues_os_own$p.value
  pairwise_pvalues_os_own_pvals_text <- paste0("Pairwise p-values:\n",
                                               "High vs Medium: ", signif(pairwise_pvalues_os_own_pvals[2, 2], 3), "\n",
                                               "High vs Low: ", signif(pairwise_pvalues_os_own_pvals[2, 1], 3), "\n",
                                               "Medium vs Low: ", signif(pairwise_pvalues_os_own_pvals[1, 1], 3))
  ## Plot survival curve with dynamic title
  p <- ggsurvplot(fit, title = gene_col,  pval = TRUE, conf.int = FALSE,
                  risk.table = TRUE, risk.table.y.text.col = TRUE, 
                  palette = c("#3498DB", "#F39C12", "#E74C3C"),
                  xlab = "Time in Months",  ylab = "OS Probability",
                  break.time.by = 4, ggtheme = theme_light(), 
                  risk.table.height = 0.25,
                  surv.median.line = "hv", xlim = c(0, 48))
  p$plot <- p$plot +
    annotate("text", x = max_months - 20, y = 0.8,
             label = pairwise_pvalues_os_own_pvals_text, size = 3, hjust = 0)
  ## Save
  output_file <- paste0(output_folder_os, "/", gene_col, "_OS_survival_plot.png")
  ggsave(output_file, plot = p$plot, width = 8, height = 6)
  ## Print the plot
  print(p)
}

output_folder_pfs <- "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Res/Images/Survival_Test_Single_Gene/PFS"
## PFS
g <- colnames(df_tpm_own_filtered_merged[,1:30])  ## Get the column names
for(i in seq_along(g)){  ## Iterate over the column names by index
  gene_col <- g[i]  ## Get the current gene column name
  ## Create Group based on the median of the Mean column
  df_tpm_own_filtered_merged$Group <- cut(df_tpm_own_filtered_merged[[gene_col]],
                                          breaks = quantile(df_tpm_own_filtered_merged[[gene_col]],
                                                            probs = c(0, 1/3, 2/3, 1)),
                                          labels = c("Low", "Medium", "High"),
                                          include.lowest = TRUE)
  ## Perform survival analysis
  fit <- survfit(Surv(PFS_time, as.numeric(PFS_event)) ~ Group, data = df_tpm_own_filtered_merged)
  pairwise_pvalues_pfs_own <- pairwise_survdiff(Surv(PFS_time, as.integer(PFS_event)) ~ Group,
                                                data = df_tpm_own_filtered_merged,
                                                p.adjust.method = "none")
  pairwise_pvalues_pfs_own_pvals <- pairwise_pvalues_pfs_own$p.value
  pairwise_pvalues_pfs_own_pvals_text <- paste0("Pairwise p-values:\n",
                                               "High vs Medium: ", signif(pairwise_pvalues_pfs_own_pvals[2, 2], 3), "\n",
                                               "High vs Low: ", signif(pairwise_pvalues_pfs_own_pvals[2, 1], 3), "\n",
                                               "Medium vs Low: ", signif(pairwise_pvalues_pfs_own_pvals[1, 1], 3))
  ## Plot survival curve with dynamic title
  p <- ggsurvplot(fit, title = gene_col,  pval = TRUE, conf.int = FALSE,
                  risk.table = TRUE, risk.table.y.text.col = TRUE, 
                  palette = c("#3498DB", "#F39C12", "#E74C3C"),
                  xlab = "Time in Months",  ylab = "PFS Probability",
                  break.time.by = 4, ggtheme = theme_light(), 
                  risk.table.height = 0.25,
                  surv.median.line = "hv", xlim = c(0, 48))
  p$plot <- p$plot +
    annotate("text", x = max_months - 20, y = 0.8,
             label = pairwise_pvalues_pfs_own_pvals_text, size = 3, hjust = 0)
  ## Save
  output_file <- paste0(output_folder_pfs, "/", gene_col, "_PFS_survival_plot.png")
  ggsave(output_file, plot = p$plot, width = 8, height = 6)
  ## Print the plot
  print(p)
}



## *******************************************************************************+++
## Testing single gene also for TCGA
df_clinical_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
df_survival_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.survival.tsv"), show_col_types = FALSE)
df_tpm_tcga      <- read_tsv(paste0(data_path, "/TCGA-STAD.star_tpm.tsv"), show_col_types = FALSE)

## Filtering stage
df_clinical_tcga <- df_clinical_tcga %>%
  filter(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA",
                                                "Stage IIIB", "Stage IIIC", "Stage IV"))
ids_df_clinical_tcga <- df_clinical_tcga$sample
df_survival_tcga     <- df_survival_tcga %>%
  filter(sample %in% ids_df_clinical_tcga)
## Filtering
top_genes_id <- df_tpm_own %>% 
  dplyr::select(gene_id, gene_name) %>% 
  filter(gene_name %in% top_genes)
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
selected_rownames_tcga <- rownames(df_tpm_tcga_cleaned)[grep("01A$", rownames(df_tpm_tcga_cleaned))]
## Filtering
df_tpm_tcga_cleaned_tumor <- df_tpm_tcga_cleaned %>% 
  dplyr::filter(rownames(df_tpm_tcga_cleaned) %in% selected_rownames_tcga)
## Remove the last part of barcode
selected_rownames_cleaned_tcga <- sub("-[^-]+$", "", rownames(df_tpm_tcga_cleaned_tumor))
rownames(df_tpm_tcga_cleaned_tumor) <- selected_rownames_cleaned_tcga
## Merge
df_tpm_tcga_cleaned_tumor$`_PATIENT` <- rownames(df_tpm_tcga_cleaned_tumor)
df_tpm_tcga_cleaned_filtered_merged  <- merge(df_tpm_tcga_cleaned_tumor,
                                              df_survival_tcga, by = "_PATIENT")


## Survival analysis
output_folder_os_tcga <- "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Res/Images/Survival_Test_Single_Gene/OS_TCGA"
## OS
g <- colnames(df_tpm_tcga_cleaned_filtered_merged[,2:31])  ## Get the column names
for(i in seq_along(g)){  ## Iterate over the column names by index
  gene_col <- g[i]  ## Get the current gene column name
  ## Create Group based on the median of the Mean column
  df_tpm_tcga_cleaned_filtered_merged$Group <- cut(df_tpm_tcga_cleaned_filtered_merged[[gene_col]],
                                          breaks = quantile(df_tpm_tcga_cleaned_filtered_merged[[gene_col]],
                                                            probs = c(0, 1/3, 2/3, 1)),
                                          labels = c("Low", "Medium", "High"),
                                          include.lowest = TRUE)
  
  ## Perform survival analysis
  fit <- survfit(Surv(OS.time/30, OS) ~ Group, data = df_tpm_tcga_cleaned_filtered_merged)
  pairwise_pvalues_os_tcga <- pairwise_survdiff(Surv(OS.time/30, OS) ~ Group,
                                                data = df_tpm_tcga_cleaned_filtered_merged,
                                                p.adjust.method = "none")
  pairwise_pvalues_os_tcga_pvals <- pairwise_pvalues_os_tcga$p.value
  pairwise_pvalues_os_tcga_pvals_text <- paste0("Pairwise p-values:\n",
                                                "High vs Medium: ", signif(pairwise_pvalues_os_tcga_pvals[2, 2], 3), "\n",
                                                "High vs Low: ", signif(pairwise_pvalues_os_tcga_pvals[2, 1], 3), "\n",
                                                "Medium vs Low: ", signif(pairwise_pvalues_os_tcga_pvals[1, 1], 3))
  ## Plot survival curve with dynamic title
  p <- ggsurvplot(fit, title = gene_col,  pval = TRUE, conf.int = FALSE,
                  risk.table = TRUE, risk.table.y.text.col = TRUE, 
                  palette = c("#3498DB", "#F39C12", "#E74C3C"),
                  xlab = "Time in Months",  ylab = "OS Probability",
                  break.time.by = 4, ggtheme = theme_light(), 
                  risk.table.height = 0.25,
                  surv.median.line = "hv", xlim = c(0, 48))
  p$plot <- p$plot +
    annotate("text", x = max_months - 20, y = 0.8,
             label = pairwise_pvalues_os_tcga_pvals_text, size = 3, hjust = 0)
  ## Save
  output_file <- paste0(output_folder_os_tcga, "/", gene_col, "_OS_survival_plot.png")  
  ggsave(output_file, plot = p$plot, width = 8, height = 6)  
  ## Print the plot
  print(p)
}
