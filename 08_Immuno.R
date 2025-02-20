## Code useful to perform Immuno Signature
source("00_Support_ssGSEA.R")
library(GSEABase)
library(GSVA)
library(limma)
library(GSEABase)
library(edgeR)
library(pheatmap)
library(decoupleR)
library(OmnipathR)

## ***************************************
## LOAD DATA
## ***************************************
## Load tcga tpm
df_tpm_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.star_tpm.tsv"),
                        show_col_types = FALSE)
## Load own tpm
df_tpm_own  <- read_tsv(paste0(data_path_own, "/salmon.merged.gene_tpm.tsv"),
                        show_col_types = FALSE)
## Signature
df_signature_own  <- read_csv(paste0(data_path_own, "/OWN_Signature.csv"))
df_signature_tcga <- read_csv(paste0(data_path_own, "/TCGA_Signature.csv"))
df_signature_own <- df_signature_own %>% 
  rename(ID = ...1)
df_signature_tcga <- df_signature_tcga %>% 
  rename(ID = ...1)
## Load survival data tcga
df_survival_tcga        <- read_tsv(paste0(output_path, "/survival_cleaned.tsv"), show_col_types = FALSE)
df_clinical_tcga        <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
## Load survival data own
df_survival_own <- read_excel(paste0(data_path_own, "/stomaci_wes_all.xlsx"))

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
## Wrangling survival data tcga
df_clinical_tcga_filtered <- df_clinical_tcga %>%
  filter(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA",
                                                "Stage IIIB", "Stage IIIC", "Stage IV"))
ids_df_clinical_tcga <- df_clinical_tcga_filtered$sample
df_survival_tcga_filtered     <- df_survival_tcga %>%
  filter(sample %in% ids_df_clinical_tcga)
## Wrangling survival data own
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

## OWN TPM
## Consider only tumor patients
df_tpm_own_filtered_tumor <- df_tpm_own %>% 
  dplyr::select(matches("T$"))

## Setting row names
rownames(df_tpm_tcga_cleaned_tumor) <- make.names(genes_names_tcga, unique = TRUE) 
rownames(df_tpm_own_filtered_tumor) <- make.names(genes_names_tpm, unique = TRUE) 



## SIGNATURES
## OWN
## Load Immuno
signature_gmt <- readr::read_delim("signature_gmt.csv", 
                                   delim = ";", escape_double = FALSE, trim_ws = TRUE)
expression_data     <- as.matrix(df_tpm_own_filtered_tumor)
ssgsea_params       <- ssgseaParam(expression_data, signature_gmt)
ssGSEA_scores       <- gsva(ssgsea_params)
ssGSEA_scores       <- as.data.frame(t(ssGSEA_scores))
ssGSEA_scores_scale <- scale(ssGSEA_scores)
ssGSEA_scores_scale <- data.frame(ssGSEA_scores_scale)
## Wrangling
ssGSEA_scores_scale$ID <- rownames(ssGSEA_scores_scale)
## Merge
ssGSEA_scores_scale_final <- merge(ssGSEA_scores_scale, df_signature_own, by = "ID")
rownames(ssGSEA_scores_scale_final) <- ssGSEA_scores_scale_final$ID
ssGSEA_scores_scale_final$Genes <- NULL 

## TCGA
## Load Immuno
signature_gmt <- readr::read_delim("signature_gmt.csv", 
                                   delim = ";", escape_double = FALSE, trim_ws = TRUE)
expression_data_tcga     <- as.matrix(df_tpm_tcga_cleaned_tumor)
ssgsea_params_tcga       <- ssgseaParam(expression_data_tcga, signature_gmt)
ssGSEA_scores_tcga       <- gsva(ssgsea_params_tcga)
ssGSEA_scores_tcga       <- as.data.frame(t(ssGSEA_scores_tcga))
ssGSEA_scores_scale_tcga <- scale(ssGSEA_scores_tcga)
ssGSEA_scores_scale_tcga <- data.frame(ssGSEA_scores_scale_tcga)
ssGSEA_scores_scale_tcga$Class <- rep("TCGA", nrow(ssGSEA_scores_scale_tcga))
## Wrangling
ssGSEA_scores_scale_tcga$ID <- rownames(ssGSEA_scores_scale_tcga)
## Merge
ssGSEA_scores_scale_final_tcga <- merge(ssGSEA_scores_scale_tcga, df_signature_tcga, by = "ID")
rownames(ssGSEA_scores_scale_final_tcga) <- ssGSEA_scores_scale_final_tcga$ID
ssGSEA_scores_scale_final_tcga$Genes <- NULL 

## Merge
df_final_immuno <- rbind(ssGSEA_scores_scale_final, ssGSEA_scores_scale_final_tcga)

## Reshape
df_final_immuno_long <- tidyr::pivot_longer(df_final_immuno, cols = colnames(df_final_immuno)[1:29], 
                                              names_to = "Variable", values_to = "Value")

## Boxplot
immuno_unique <- unique(df_final_immuno_long$Variable)
# for (immuno in immuno_unique){
#   df_subset <- subset(df_final_immuno_long, Variable == immuno)
#   p <- ggplot(df_subset,
#               aes(x = factor(Class, levels = c("IRE", "TCGA")),
#                   y = Value, fill = factor(Signature, levels = c("Low", "High")))) +
#     geom_boxplot(outlier.size = 1, outlier.colour = "black",
#                  width = 0.7, colour = "black", alpha = 0.8) +
#     labs(title = paste0(immuno),
#          x = "",
#          y = "Value Expression",
#          fill = "Signature") +
#     scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#d62728")) +
#     theme_minimal(base_size = 14) +
#     theme(legend.position = "right",
#           axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
#           axis.text.y = element_text(size = 12, face = "plain"),
#           plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
#           axis.title = element_text(size = 14, face = "bold"),
#           panel.grid.major = element_line(size = 0.2, color = "gray90"),
#           panel.grid.minor = element_blank(),
#           plot.background = element_rect(fill = "white", color = "white")) +
#     stat_compare_means(method = "kruskal.test", label = "p.signif", label.x = 1.5, size = 5)
#   ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Immuno_Signatures/",
#                 immuno, "_boxplot.png"), plot = p, width = 8, height = 6)}

## Plot - IRE
color_scale <- colorRampPalette(c("blue", "white", "red"))(100)
my_sample_col <- ssGSEA_scores_scale_final %>% 
  dplyr::select("Signature")
ssGSEA_scores_scale_final_plot <- ssGSEA_scores_scale_final %>% 
  dplyr::select(-c("ID", "Signature"))
ordered_columns <- order(my_sample_col$Signature, decreasing = FALSE)
ordered_matrix  <- ssGSEA_scores_scale_final_plot[ordered_columns, ]
ordered_annotation_col <- my_sample_col[ordered_columns, , drop = FALSE]
annotation_colors <- list(Signature = c("Low" = "#1f77b4", "High" = "#d62728"))
pheatmap(as.matrix(t(ordered_matrix)),
         annotation_col = ordered_annotation_col,
         annotation_colors = annotation_colors,
         label_row = NA,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         cluster_row = TRUE,
         main = "Immuno Signatures - IRE Cohort",
         color = color_scale,
         cutree_cols = 2)

## Plot - TCGA
my_sample_col_tcga <- ssGSEA_scores_scale_final_tcga %>% 
  dplyr::select("Signature")
ssGSEA_scores_scale_final_tcga_plot <- ssGSEA_scores_scale_final_tcga %>% 
  dplyr::select(-c("ID", "Signature"))
ordered_columns_tcga <- order(my_sample_col_tcga$Signature, decreasing = FALSE)
ordered_matrix_tcga  <- ssGSEA_scores_scale_final_tcga_plot[ordered_columns_tcga, ]
ordered_matrix_tcga$Class <- NULL
ordered_annotation_col_tcga <- my_sample_col_tcga[ordered_columns_tcga, , drop = FALSE]
pheatmap(as.matrix(t(ordered_matrix_tcga)),
         annotation_col = ordered_annotation_col_tcga,
         annotation_colors = annotation_colors,
         label_row = NA,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         cluster_row = TRUE,
         main = "Immuno Signatures - TCGA Cohort",
         color = color_scale)


## Pecentage of mutated genes





































