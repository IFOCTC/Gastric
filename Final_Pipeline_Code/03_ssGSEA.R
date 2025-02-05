## Code useful to perform GSVA
source("00_Support_ssGSEA.R")

## ***************************************
## LOAD DATA
## ***************************************
## Load tcga tpm
df_tpm_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.star_tpm.tsv"),
                        show_col_types = FALSE)
## Load own tpm
df_tpm_own  <- read_tsv(paste0(data_path_own, "/salmon.merged.gene_tpm.tsv"),
                        show_col_types = FALSE)

## ***************************************
## WRANGLING
## ***************************************
## Log transformation
exclude_cols <- c("gene_id", "gene_name")
df_tpm_own_filtered <- df_tpm_own %>% 
  dplyr::select(-all_of(exclude_cols)) %>%
  mutate(across(everything(), ~ log2(. + 1)))
df_tpm_own <- df_tpm_own %>%
  dplyr::select(all_of(exclude_cols)) %>%
  bind_cols(df_tpm_own_filtered)
## Load survival data tcga
df_survival_tcga <- read_tsv(paste0(output_path, "/survival_cleaned.tsv"), show_col_types = FALSE)
df_clinical_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
df_clinical_tcga_filtered <- df_clinical_tcga %>%
  filter(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA",
                                                "Stage IIIB", "Stage IIIC", "Stage IV"))
ids_df_clinical_tcga <- df_clinical_tcga_filtered$sample
df_survival_tcga_filtered     <- df_survival_tcga %>%
  filter(sample %in% ids_df_clinical_tcga)
## Load survival data own
df_survival_own <- read_excel(paste0(data_path_own, "/stomaci_wes_all.xlsx"))
df_survival_own_filtered_os  <- df_survival_own %>%
  filter(CT1L_FLOT_1SI == 1)
df_survival_own_filtered_pfs <- df_survival_own %>%
  filter(CT1L_1SI == 1)
## Ensembl id gene name
Ensembl_ID                     <- df_tpm_tcga$Ensembl_ID
Ensembl_ID_cleaned             <- sub("\\..*", "", Ensembl_ID) 
df_tpm_tcga_cleaned            <- df_tpm_tcga
df_tpm_tcga_cleaned$Ensembl_ID <- Ensembl_ID_cleaned
## Attach the gene name of own genes to tpm tcga
df_support_gene_name <- df_tpm_own %>% 
  dplyr::select(gene_id, gene_name)
## Change name column for merge
df_support_gene_name <- df_support_gene_name %>% 
  rename(Ensembl_ID = gene_id)
## Merge
df_tpm_tcga_final <- df_tpm_tcga_cleaned %>% 
  left_join(df_support_gene_name, by = "Ensembl_ID")
## Move the column
df_tpm_tcga_final <- df_tpm_tcga_final %>% 
  dplyr::select(1, gene_name, everything())
## Remove gene id
df_tpm_tcga_final$Ensembl_ID <- NULL
df_tpm_own$gene_id <- NULL
## Export row names
genes_names_tcga <- df_tpm_tcga_final$gene_name
genes_names_tpm  <- df_tpm_own$gene_name
## Remove column
df_tpm_tcga_final$gene_name <- NULL
df_tpm_own$gene_name        <- NULL

## TCGA TPM
## Consider only tumor samples
selected_colnames_tcga <- colnames(df_tpm_tcga_final)[grep("01A$", colnames(df_tpm_tcga_final))]
## Filtering
df_tpm_tcga_cleaned_tumor <- df_tpm_tcga_final %>% 
  dplyr::select(any_of(selected_colnames_tcga))
## Remove the last part of barcode
selected_colnames_cleaned_tcga <- sub("-[^-]+$", "", colnames(df_tpm_tcga_cleaned_tumor))
colnames(df_tpm_tcga_cleaned_tumor) <- selected_colnames_cleaned_tcga
## Filter Stage III-IV
df_survival_tcga_filtered$sample <-  sub("-[^-]+$", "", df_survival_tcga_filtered$sample)
common_columns <- intersect(colnames(df_tpm_tcga_cleaned_tumor), df_survival_tcga_filtered$sample)
df_tpm_tcga_cleaned_tumor_filtered <- df_tpm_tcga_cleaned_tumor %>% 
  dplyr::select(all_of(common_columns))

## OWN TPM
## Consider only tumor patients
df_tpm_own_filtered_tumor <- df_tpm_own %>% 
  dplyr::select(matches("T$"))

## Setting row names
rownames(df_tpm_tcga_cleaned_tumor) <- make.names(genes_names_tcga, unique = TRUE) 
rownames(df_tpm_own_filtered_tumor) <- make.names(genes_names_tpm, unique = TRUE) 
rownames(df_tpm_tcga_cleaned_tumor_filtered) <- make.names(genes_names_tcga, unique = TRUE)

## Define the list of genes to use
genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3",
           "MMP14", "FYN", "TUBA1A", "TNFRSF1A", "CTSZ",
           "CLIC4", "TNIP1", "C1S", "HIF1A", "PGK1",
           "PIP4K2A", "LGALS1", "SPARC", "SERPINH1", "ECE1",
           "FBN1", "TMSB10", "DCTN2", "PDLIM1", "GNAS",
           "DBN1", "IGFBP4", "CTSD", "C3", "NONO",
           "ITGB1", "LTBP3", "MYADM", "COL18A1", "CAV1",
           "CD68", "CUL4B", "C1R", "SPON2", "ABCA1",
           "TAF10", "EMP3", "CSGALNACT2", "CD44", "GRN",
           "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
           "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2",
           "SMAD3", "FAM3C", "VIM", "MORF4L2", "NRBP1",
           "PRNP", "TUBB6", "RAP1B", "CD82", "ARFGAP1",
           "NRP2", "VPS26A", "PML", "TAGLN", "SART1",
           "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
           "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
           "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
           "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
           "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP",
           "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8",
           "CYR61", "LGALS9")
## Combine dataframe
genes_set <- data.frame(Genes = genes)
genes_set <- as.list(genes_set)


## ***************************************
## SSGSEA
## ***************************************
## Compute own tpm
system.time(assign("res_own_tpm", ssgsea(as.matrix(df_tpm_own_filtered_tumor), genes_set,
                                         scale = TRUE, norm = TRUE)))
## Compute tcga tpm
system.time(assign("res_tcga_tpm", ssgsea(as.matrix(df_tpm_tcga_cleaned_tumor), genes_set,
                                          scale = TRUE, norm = TRUE)))
## Transpose
res_own_tpm_transposed  <- t(res_own_tpm)
res_tcga_tpm_transposed <- t(res_tcga_tpm)
## Z-Score the ssgsea output for comparative analysis
mat_own_tpm  <- (res_own_tpm - rowMeans(res_own_tpm))/
  (rowSds(as.matrix(res_own_tpm)))[row(res_own_tpm)]
mat_tcga_tpm <- (res_tcga_tpm - rowMeans(res_tcga_tpm))/
  (rowSds(as.matrix(res_tcga_tpm)))[row(res_tcga_tpm)]
## Transpose
mat_own_tpm_transposed  <- t(mat_own_tpm)
mat_tcga_tpm_transposed <- t(mat_tcga_tpm)
## Df
mat_own_tpm_transposed  <- data.frame(mat_own_tpm_transposed)
mat_tcga_tpm_transposed <- data.frame(mat_tcga_tpm_transposed)

## ***************************************
## SURVIVAL ANALYSIS
## ***************************************
## Set t
max_months <- 48
t <- 0
mat_own_tpm_transposed <- mat_own_tpm_transposed %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))
mat_tcga_tpm_transposed <- mat_tcga_tpm_transposed %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))

# ## Export
# write.csv(mat_own_tpm_transposed, paste(data_path, "OWN_Signature.csv"),
#           row.names = TRUE)
# write.csv(mat_tcga_tpm_transposed, paste(data_path, "TCGA_Signature.csv"),
#           row.names = TRUE)

## Boxplot molecular subtypes by signature

## SCD
## CLDN18 Distribution in High & Low
df_tpm_own_filtered_tumor_transposed <- t(df_tpm_own_filtered_tumor)
df_tpm_own_filtered_tumor_transposed <- data.frame(df_tpm_own_filtered_tumor_transposed)
df_tpm_own_filtered_tumor_transposed$Signature <- mat_own_tpm_transposed$Signature
## OWN
ggplot(df_tpm_own_filtered_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
                    y = SCD, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_own_filtered_tumor_transposed$CLDN18)) +
  labs(title = "SCD", subtitle = "IRE Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "wilcox.test", label.x = 0.6)
## TCGA
df_tpm_tcga_cleaned_tumor_transposed <- t(df_tpm_tcga_cleaned_tumor)
df_tpm_tcga_cleaned_tumor_transposed <- data.frame(df_tpm_tcga_cleaned_tumor_transposed)
df_tpm_tcga_cleaned_tumor_transposed$Signature <- mat_tcga_tpm_transposed$Signature
## Stage III-IV
df_tpm_tcga_cleaned_tumor_transposed <- df_tpm_tcga_cleaned_tumor_transposed %>% 
  filter(rownames(df_tpm_tcga_cleaned_tumor_transposed) %in% df_clinical_tcga_filtered$submitter_id)
## Plot
ggplot(df_tpm_tcga_cleaned_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = SCD, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_tcga_cleaned_tumor_transposed$SCD)) +
  labs(title = "SCD", subtitle = "TCGA Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "wilcox.test")

## TACSTD2 Distribution in High & Low
## OWN
ggplot(df_tpm_own_filtered_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = TACSTD2, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_own_filtered_tumor_transposed$TACSTD2)) +
  labs(title = "TACSTD2", subtitle = "IRE Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6)
## TCGA
mat_tcga_tpm_transposed$ID <- rownames(mat_tcga_tpm_transposed)
df_tpm_tcga_cleaned_tumor_transposed$submitter_id <- rownames(df_tpm_tcga_cleaned_tumor_transposed)
df_stage <- merge(df_tpm_tcga_cleaned_tumor_transposed,
                  mat_tcga_tpm_transposed, by = "submitter_id", how = "inner")
df_final_stage <- merge(df_stage, df_clinical_tcga, by = "submitter_id",
                        how = "outer")
## Plot
ggplot(df_final_stage,
       aes(x = ajcc_pathologic_stage.diagnoses,
           y = TACSTD2, fill = Signature.x)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_final_stage$TACSTD2)) +
  labs(title = "TACSTD2", subtitle = "TCGA Cohort Analysis", 
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6)

## ERBB2 Distribution in High & Low
## OWN
ggplot(df_tpm_own_filtered_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = ERBB2, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_own_filtered_tumor_transposed$ERBB2)) +
  labs(title = "ERBB2", subtitle = "IRE Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6)

## TCGA
ggplot(df_tpm_tcga_cleaned_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = ERBB2, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_tcga_cleaned_tumor_transposed$ERBB2)) +
  labs(title = "ERBB2", subtitle = "TCGA Cohort Analysis", 
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6)

## HIF1A Distribution in High & Low
## OWN
ggplot(df_tpm_own_filtered_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = HIF1A, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) +  
  ylim(0, max(df_tpm_own_filtered_tumor_transposed$HIF1A)) +
  labs(title = "HIF1A", subtitle = "IRE Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), 
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6, label.y = 10.1)
## TCGA
ggplot(df_tpm_tcga_cleaned_tumor_transposed,
       aes(x = factor(Signature, levels = c("Low", "High")),
           y = HIF1A, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", alpha = 0.8) + 
  ylim(0, max(df_tpm_tcga_cleaned_tumor_transposed$HIF1A)) +
  labs(title = "HIF1A", subtitle = "TCGA Cohort Analysis",
       x = "Signature",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "",  
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 14, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_line(size = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")) +
  stat_compare_means(method = "kruskal.test", label.x = 0.6, label.y = 9.1)

## Check
table(mat_own_tpm_transposed$Signature)
table(mat_tcga_tpm_transposed$Signature)

## Wrangling for survival analysis
## TCGA Tpm
mat_tcga_tpm_transposed$ID <- rownames(mat_tcga_tpm_transposed)
mat_tcga_tpm_transposed_merged <- merge(mat_tcga_tpm_transposed, df_survival_tcga, by = "ID")
mat_tcga_tpm_transposed_merged_filtered <- merge(mat_tcga_tpm_transposed,
                                                 df_survival_tcga_filtered, by = "ID")
mat_tcga_tpm_transposed_merged$OS.time.month <- mat_tcga_tpm_transposed_merged$OS.time / 30
mat_tcga_tpm_transposed_merged_filtered$OS.time.month <- mat_tcga_tpm_transposed_merged_filtered$OS.time / 30

## Own TPM
mat_own_tpm_transposed$ID <- rownames(mat_own_tpm_transposed)
mat_own_tpm_transposed$ID <- gsub("\\.1T$", "", mat_own_tpm_transposed$ID)
mat_own_tpm_transposed$ID <- gsub("\\.", "-", mat_own_tpm_transposed$ID)
mat_own_tpm_transposed_merged_os  <- merge(mat_own_tpm_transposed, df_survival_own_filtered_os, by = "ID")
mat_own_tpm_transposed_merged_pfs <- merge(mat_own_tpm_transposed, df_survival_own_filtered_pfs, by = "ID")

## Signature Density Distribution
## OS IRE 
p_signature_os <- ggplot(mat_own_tpm_transposed_merged_os, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature OS IRE") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## PFS IRE 
p_signature_pfs <- ggplot(mat_own_tpm_transposed_merged_pfs, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature PFS IRE") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## TCGA
p_signature_tcga <- ggplot(mat_tcga_tpm_transposed_merged, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature TCGA") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## TCGA Filtered
p_signature_filtered_tcga <- ggplot(mat_tcga_tpm_transposed_merged_filtered,
                                    aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature Filetered TCGA") +
  xlab("") +
  ylab("Count") +
  theme_minimal()

## Signature Barplot Distribution
## OS Status distribution
p_os_status_distribution_ire <- ggplot(mat_own_tpm_transposed_merged_os, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),
           fill = c("#E74C3C", "#3498DB"), color = c("#E74C3C", "#3498DB"), alpha = 0.8) + 
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -0.5, size = 5, fontface = "bold") + 
  scale_y_continuous(labels = percent_format()) + 
  labs(title = "OS Status Distribution - IRE", 
       subtitle = "CT1L_FLOT_1SI",
       x = "Signature",
       y = "Percentage") +
  theme_classic(base_size = 14, base_family = "Arial") + 
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),  
        axis.title.x = element_text(vjust = -0.5),  
        axis.title.y = element_text(vjust = 1.5),   
        axis.text.x = element_text(face = "bold"),  
        axis.text.y = element_text(face = "bold"), 
        panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),  
        legend.position = "none")
## PFS Status Distribution
p_pfs_status_distribution_ire <- ggplot(mat_own_tpm_transposed_merged_pfs, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),
           fill = c("#E74C3C", "#3498DB"), color = c("#E74C3C", "#3498DB"), alpha = 0.8) + 
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -0.5, size = 5, fontface = "bold") + 
  scale_y_continuous(labels = percent_format()) + 
  labs(title = "PFS Status Distribution - IRE", 
       subtitle = "CT1L_1SI",
       x = "Signature",
       y = "Percentage") +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),  
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 1.5),  
        axis.text.x = element_text(face = "bold"),  
        axis.text.y = element_text(face = "bold"), 
        panel.grid.major = element_line(color = "gray90"),  
        panel.grid.minor = element_blank(),  
        legend.position = "none")
## TCGA TPM Plot with Percentage on the Y-axis
p_os_status_distribution_tcga <- ggplot(mat_tcga_tpm_transposed_merged, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  
           fill = c("#E74C3C", "#3498DB"), color = c("#E74C3C", "#3498DB"), alpha = 0.8) +  
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -0.5, size = 5, fontface = "bold") +  
  scale_y_continuous(labels = scales::percent_format()) +  
  labs(title = "OS Status Distribution - TCGA", 
       subtitle = "Whole Population",
       x = "Signature", 
       y = "Percentage") +
  theme_classic(base_size = 14, base_family = "Arial") +  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5), 
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 1.5),  
        axis.text.x = element_text(face = "bold"),  
        axis.text.y = element_text(face = "bold"),  
        panel.grid.major = element_line(color = "gray90"),  
        panel.grid.minor = element_blank(),  
        legend.position = "none")
## OS Status Distribution Filtered Plot for TCGA
p_os_status_distribution_filtered_tcga <- ggplot(mat_tcga_tpm_transposed_merged_filtered, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  
           fill = c("#E74C3C", "#3498DB"), color = c("#E74C3C", "#3498DB"), alpha = 0.8) +  
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -0.5, size = 5, fontface = "bold") +  
  scale_y_continuous(labels = scales::percent_format()) + 
  labs(title = "OS Status Distribution - TCGA", 
       subtitle = "Stage III/IV",
       x = "Signature", 
       y = "Percentage") +
  theme_classic(base_size = 14, base_family = "Arial") + 
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),  
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1.5), 
        axis.text.x = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90"),  
        panel.grid.minor = element_blank(),  
        legend.position = "none")


## ***************************************
## Survival Analysis - IRE
## ***************************************
fit_os_own <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Signature,
                      data = mat_own_tpm_transposed_merged_os)
## Plot
os_plot_own <- ggsurvplot(fit_os_own, data = mat_own_tpm_transposed_merged_os,
                          risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                          palette = c("#E74C3C", "#3498DB"),  
                          xlim = c(0, max_months),  
                          title = "OS - IRE",  
                          subtitle = paste0("N. Genes: ", length(genes_set$Gene), ", t: ", t, " - CT1L_FLOT_1SI"), 
                          xlab = "Time (Months)",  ylab = "Survival Probability", 
                          break.time.by = 4,  
                          ggtheme = theme_classic(base_size = 16), 
                          risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                          risk.table.y.text = TRUE, conf.int.style = "step",  
                          surv.median.line = "hv",  
                          pval.size = 6,  
                          font.title = c(18, "bold"), 
                          font.subtitle = c(14, "italic"), 
                          font.x = c(14, "plain"),  
                          font.y = c(14, "plain"),  
                          font.tickslab = c(12, "plain"),  
                          legend.title = "Group",  
                          legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                          legend = c(0.85, 0.85),  
                          risk.table.fontsize = 3.5)


## ***************************************
## Progression Free Survival - IRE
## ***************************************
fit_pfs_own <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Signature,
                       data = mat_own_tpm_transposed_merged_pfs)
pfs_plot_own <- ggsurvplot(fit_pfs_own, data = mat_own_tpm_transposed_merged_pfs,
                           risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                           palette = c("#E74C3C", "#3498DB"),  
                           xlim = c(0, max_months),  
                           title = "PFS - IRE",  
                           subtitle = paste0("N. Genes: ", length(genes_set$Gene), ", t: ", t, " - CT1L_1SI"), 
                           xlab = "Time (Months)",  ylab = "Survival Probability", 
                           break.time.by = 4,  
                           ggtheme = theme_classic(base_size = 16), 
                           risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                           risk.table.y.text = TRUE, conf.int.style = "step",  
                           surv.median.line = "hv",  
                           pval.size = 6,  
                           font.title = c(18, "bold"), 
                           font.subtitle = c(14, "italic"), 
                           font.x = c(14, "plain"),  
                           font.y = c(14, "plain"),  
                           font.tickslab = c(12, "plain"),  
                           legend.title = "Group",  
                           legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                           legend = c(0.85, 0.85),  
                           risk.table.fontsize = 3.5)

## ***************************************
## Survival Analysis - TCGA 
## ***************************************
fit_os_tcga <- survfit(Surv(OS.time.month, as.integer(OS)) ~ Signature,
                       data = mat_tcga_tpm_transposed_merged)
os_plot_tcga <- ggsurvplot(fit_os_tcga, data = mat_tcga_tpm_transposed_merged,
                           risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                           palette = c("#E74C3C", "#3498DB"),  
                           xlim = c(0, max_months),  
                           title = "OS - TCGA",  
                           subtitle = paste0("N. Genes: ", length(genes_set$Gene), ", t: ", t, " - Whole Population"), 
                           xlab = "Time (Months)",  ylab = "Survival Probability", 
                           break.time.by = 4,  
                           ggtheme = theme_classic(base_size = 16), 
                           risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                           risk.table.y.text = TRUE, conf.int.style = "step",  
                           surv.median.line = "hv",  
                           pval.size = 6,  
                           font.title = c(18, "bold"), 
                           font.subtitle = c(14, "italic"), 
                           font.x = c(14, "plain"),  
                           font.y = c(14, "plain"),  
                           font.tickslab = c(12, "plain"),  
                           legend.title = "Group",  
                           legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                           legend = c(0.85, 0.85),  
                           risk.table.fontsize = 3.5)

## ***************************************
## Survival Analysis - TCGA Filtered
## ***************************************
fit_os_tcga_filtered <- survfit(Surv(OS.time.month, as.integer(OS)) ~ Signature,
                                data = mat_tcga_tpm_transposed_merged_filtered)
os_plot_tcga_filtered <- ggsurvplot(fit_os_tcga_filtered, data = mat_tcga_tpm_transposed_merged_filtered,
                                    risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                                    palette = c("#E74C3C", "#3498DB"),  
                                    xlim = c(0, max_months),  
                                    title = "OS - TCGA",  
                                    subtitle = paste0("N. Genes: ", length(genes_set$Gene), ", t: ", t, " - Stage III/IV"), 
                                    xlab = "Time (Months)",  ylab = "Survival Probability", 
                                    break.time.by = 4,  
                                    ggtheme = theme_classic(base_size = 16), 
                                    risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                                    risk.table.y.text = TRUE, conf.int.style = "step",  
                                    surv.median.line = "hv",  
                                    pval.size = 6,  
                                    font.title = c(18, "bold"), 
                                    font.subtitle = c(14, "italic"), 
                                    font.x = c(14, "plain"),  
                                    font.y = c(14, "plain"),  
                                    font.tickslab = c(12, "plain"),  
                                    legend.title = "Group",  
                                    legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                                    legend = c(0.85, 0.85),  
                                    risk.table.fontsize = 3.5)

## Plot everything together
blank <- grid::nullGrob()
grid.arrange(p_os_status_distribution_ire, os_plot_own$plot, p_pfs_status_distribution_ire, pfs_plot_own$plot,
             p_os_status_distribution_tcga, os_plot_tcga$plot, p_os_status_distribution_filtered_tcga, os_plot_tcga_filtered$plot, 
             ncol = 4)


## ***************************************
## EDA
## ***************************************
## Boxplot of expression conditioning on all stages
names(df_clinical_tcga)[names(df_clinical_tcga) == "sample"] <- "ID"
df_clinical_tcga_plot <- df_clinical_tcga %>%
  filter(ID %in% selected_colnames_tcga)
df_clinical_tcga_plot$ID <- sub("-[^-]+$", "", df_clinical_tcga_plot$ID)
## Merge
df_tcga_boxplot_stages <- merge(mat_tcga_tpm_transposed_merged, df_clinical_tcga_plot, by = "ID")
## Group by
df_plot_stage <- df_tcga_boxplot_stages %>%
  filter(!is.na(ajcc_pathologic_stage.diagnoses)) %>%
  mutate(ajcc_pathologic_stage.diagnoses = recode(ajcc_pathologic_stage.diagnoses,
                                                  "Stage I"    = "Stage I-II",
                                                  "Stage IA"   = "Stage I-II",
                                                  "Stage IB"   = "Stage I-II",
                                                  "Stage II"   = "Stage I-II",
                                                  "Stage IIA"  = "Stage I-II",
                                                  "Stage IIB"  = "Stage I-II",
                                                  "Stage III"  = "Stage III-IV",
                                                  "Stage IIIA" = "Stage III-IV",
                                                  "Stage IIIB" = "Stage III-IV",
                                                  "Stage IIIC" = "Stage III-IV",
                                                  "Stage IV"   = "Stage III-IV"))

## Boxplot Stage I-II Vs. Stage III-Iv
custom_colors <- c("Stage I-II"   = "steelblue",
                   "Stage III-IV" = "red")
ggplot(df_plot_stage, aes(x = ajcc_pathologic_stage.diagnoses, y = Genes, fill = ajcc_pathologic_stage.diagnoses)) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.7, outlier.shape = NA) +
  labs(title = "Gene Expression by AJCC Pathologic Stage",
       x = "AJCC Pathologic Stage", y = "Gene Expression",
       fill = "Stage", subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = custom_colors) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank()) +
  # stat_compare_means(comparison = my_comparisons,
  #                    label = "p.signif", method = "wilcox.test", label.y = 2.2, size = 5,
  #                    ref.group = "Stage IV") + 
  stat_compare_means(method = "kruskal.test", label.y = 2.2) +
  geom_jitter(shape = 16, color = "black", size = 1, width = 0.15, alpha = 0.5)

## Molecular Subtype TCGA
## Load clinical from bioportal where we consider the subtype variable
df_clinical_bioportal <- read_tsv(paste0(data_path, "/stad_tcga_pan_can_atlas_2018_clinical_data.tsv"),
                                  show_col_types = FALSE)

df_clinical_bioportal_filtered <- df_clinical_bioportal %>%
  filter(`Patient ID` %in% mat_tcga_tpm_transposed_merged$ID)
names(df_clinical_bioportal_filtered)[names(df_clinical_bioportal_filtered) == "Patient ID"] <- "ID"
df_clinical_bioportal_filtered <- merge(mat_tcga_tpm_transposed_merged, df_clinical_bioportal_filtered, by = "ID")
names(df_clinical_bioportal_filtered)[names(df_clinical_bioportal_filtered) == "Neoplasm Disease Stage American Joint Committee on Cancer Code"] <- "Stage"

## Check
table(df_clinical_bioportal_filtered$Subtype)

## Filter NA
df_clinical_bioportal_filtered_plot <- df_clinical_bioportal_filtered %>%
  filter(!is.na(Subtype))

## Tests and Plots
## Calculate percentage within each Signature category
df_percent <- df_clinical_bioportal_filtered_plot %>%
  group_by(Signature) %>%
  count(Subtype) %>%
  mutate(Percentage = n / sum(n))

## Barplot Molecular Subtypes with p-value and position fill
ggplot(df_percent, aes(x = Signature, y = Percentage, fill = Subtype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7, alpha = 0.85) +
  labs(title = "Molecular Subtype Distribution by Signature",
       x = "Signature",
       y = "Percentage",
       fill = "Molecular Subtype", subtitle = "TCGA Cohort Analysis") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.7)) +
  scale_fill_manual(values = c("STAD_CIN" = "steelblue", "STAD_EBV" = "firebrick", "STAD_GS" = "forestgreen",
                               "STAD_MSI" = "purple", "STAD_POLE" = "darkorange"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank())   

## Conditioning
## Remove NA values for relevant columns
df_clinical_bioportal_filtered_plot_no_na <- df_clinical_bioportal_filtered_plot %>%
  filter(!is.na(Signature), !is.na(Subtype), !is.na(Genes), !is.na(Stage))

## Define pairwise comparisons
## Plot
ggplot(df_clinical_bioportal_filtered_plot_no_na, aes(x = Signature, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Signature",
       x = "Signature", y = "Genes Expression", fill = "Subtype", subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_EBV" = "#d62728",   
                               "STAD_GS" = "#2ca02c",   
                               "STAD_MSI" = "#9467bd",   
                               "STAD_POLE" = "#ff7f0e"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(comparison = list(c("High", "Low")),
                     method = "wilcox.test", label.y = 2.5) +
  stat_compare_means(label = "p.signif")
## Molecular Subtype by Signature
ggplot(df_clinical_bioportal_filtered_plot_no_na, aes(x = Subtype, y = Genes, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Molecular Subtype",
       x = "Subtype", y = "Genes Expression", fill = "Signature", subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "#d62728")) +
  scale_x_discrete(labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(label = "p.signif", method = "kruskal.test")
## Conditioning for stage
## Group by for stage
df_clinical_bioportal_filtered_plot_no_na_stage <- df_clinical_bioportal_filtered_plot_no_na %>%
  mutate(Stage = recode(Stage,
                        "STAGE I"    = "STAGE I-II",
                        "STAGE IA"   = "STAGE I-II",
                        "STAGE IB"   = "STAGE I-II",
                        "STAGE II"   = "STAGE I-II",
                        "STAGE IIA"  = "STAGE I-II",
                        "STAGE IIB"  = "STAGE I-II",
                        "STAGE III"  = "STAGE III-IV",
                        "STAGE IIIA" = "STAGE III-IV",
                        "STAGE IIIB" = "STAGE III-IV",
                        "STAGE IIIC" = "STAGE III-IV",
                        "STAGE IV"   = "STAGE III-IV"))
## Plot
stage_separated <- ggplot(df_clinical_bioportal_filtered_plot_no_na_stage, aes(x = Stage, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Stage",
       x = "Stage", y = "Genes Expression", fill = "Subtype", subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_EBV" = "#d62728",   
                               "STAD_GS" = "#2ca02c",   
                               "STAD_MSI" = "#9467bd",   
                               "STAD_POLE" = "#ff7f0e"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(label = "p.format")
## Filter for only Stage I-II
df_clinical_bioportal_filtered_plot_no_na_stage_stage_I_II <- df_clinical_bioportal_filtered_plot_no_na_stage %>% 
  filter(Stage %in% "STAGE I-II")
df_clinical_bioportal_filtered_plot_no_na_stage_stage_III_IV <- df_clinical_bioportal_filtered_plot_no_na_stage %>% 
  filter(Stage %in% "STAGE III-IV")
## Plot STAGE I-II
stage_I_II <- ggplot(df_clinical_bioportal_filtered_plot_no_na_stage_stage_I_II, aes(x = Subtype, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Stage",
       x = "Molecular Subtype", y = "Genes Expression", fill = "Subtype", subtitle = "STAGE I-II - TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_EBV" = "#d62728",   
                               "STAD_GS" = "#2ca02c",   
                               "STAD_MSI" = "#9467bd",   
                               "STAD_POLE" = "#ff7f0e"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  scale_x_discrete(labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(ref.group = "STAD_GS", label = "p.signif", method = "wilcox.test") +
  stat_compare_means(label.y = 2.2)
## Plot STAGE III-IV
stage_III_IV <- ggplot(df_clinical_bioportal_filtered_plot_no_na_stage_stage_III_IV, aes(x = Subtype, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Stage",
       x = "Molecular Subtype", y = "Genes Expression", fill = "Subtype", subtitle = "STAGE III-IV - TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_EBV" = "#d62728",   
                               "STAD_GS" = "#2ca02c",   
                               "STAD_MSI" = "#9467bd",   
                               "STAD_POLE" = "#ff7f0e"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  scale_x_discrete(labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(ref.group = "STAD_GS", label = "p.signif", method = "wilcox.test") +
  stat_compare_means(label.y = 2.2)
layout <- rbind(c(1, 1), c(2, 3))
grid.arrange(stage_separated, stage_I_II, stage_III_IV, layout_matrix = layout)
    

## Filter
df_clinical_bioportal_filtered_plot_no_na_high <- df_clinical_bioportal_filtered_plot_no_na %>% 
  filter(Signature == "High")
## Plot
ggplot(df_clinical_bioportal_filtered_plot_no_na_high, aes(x = Subtype, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Molecular Subtype",
       x = "Molecular Subtype", y = "Genes Expression", fill = "Subtype",
       subtitle = "High Signature - TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_EBV" = "#d62728",   
                               "STAD_GS" = "#2ca02c",   
                               "STAD_MSI" = "#9467bd",   
                               "STAD_POLE" = "#ff7f0e"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  scale_x_discrete(labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(ref.group = "STAD_GS", label = "p.signif", method = "wilcox.test") +
  stat_compare_means(label.y = 2.2)

## Now we can compare the CIN and GS group in High-Low
df_cin_gs <- df_clinical_bioportal_filtered_plot_no_na %>% 
  filter(Subtype %in% c("STAD_CIN", "STAD_GS"))
## Plot
ggplot(df_cin_gs, aes(x = Signature, y = Genes, fill = Subtype)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Molecular Subtype",
       x = "Signature", y = "Genes Expression", fill = "Subtype", subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = c("STAD_CIN" = "#1f77b4",  
                               "STAD_GS" = "#2ca02c"),
                    labels = c("CIN", "GS")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(label.y = 2.2)
## Plot
ggplot(df_cin_gs, aes(x = Subtype, y = Genes, fill = Signature)) +
  geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
               colour = "black", 
               alpha = 0.8) +  
  labs(title = "Genes Expression by Molecular Subtype",
       x = "Molecular Subtype", y = "Genes Expression", fill = "Signature",
       subtitle = "TCGA Cohort Analysis") +
  scale_fill_manual(values = c("Low" = "#1f77b4",  
                               "High" = "red")) +
  scale_x_discrete(labels = c("CIN", "GS")) +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"), 
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 13),  
        axis.title = element_text(size = 14, face = "bold"),  
        panel.grid.major = element_line(size = 0.2, color = "gray90"), 
        panel.grid.minor = element_blank()) +
  stat_compare_means(label.y = 2.2)

