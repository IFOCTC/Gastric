## Code useful to perform Hallmark
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
df_survival_tcga <- read_tsv(paste0(output_path, "/survival_cleaned.tsv"), show_col_types = FALSE)
df_clinical_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
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

## Select only signature genes columns
df_tpm_own_filtered_tumor_signature <- df_tpm_own_filtered_tumor %>% 
  filter(rownames(df_tpm_own_filtered_tumor) %in% genes)

## Load Hallmark
hallmark_pathways   <- getGmt(paste0(data_path_own, "/h.all.v2023.2.Hs.symbols.gmt"))
expression_data     <- as.matrix(df_tpm_own_filtered_tumor_signature)
ssgsea_params       <- ssgseaParam(expression_data, hallmark_pathways)
ssGSEA_scores       <- gsva(ssgsea_params)
ssGSEA_scores       <- as.data.frame(t(ssGSEA_scores))
ssGSEA_scores_scale <- scale(ssGSEA_scores)
ssGSEA_scores_scale <- data.frame(ssGSEA_scores_scale)

## Wrangling
ssGSEA_scores_scale$ID <- rownames(ssGSEA_scores_scale)
## Merge
ssGSEA_scores_scale_final <- merge(ssGSEA_scores_scale, df_signature_own, by = "ID")
rownames(ssGSEA_scores_scale_final) <- ssGSEA_scores_scale_final$ID
ssGSEA_scores_scale_final$ID <- NULL
ssGSEA_scores_scale_final$Genes <- NULL 

## Reshape
df_long_own <- tidyr::pivot_longer(ssGSEA_scores_scale_final, cols = starts_with("HALLMARK"), 
                                   names_to = "Variable", values_to = "Value")
## Plot
hallmark_unique <- unique(df_long_own$Variable)
for (hallmark in hallmark_unique){
  df_subset <- subset(df_long_own, Variable == hallmark)
  p <- ggplot(df_subset,
              aes(x = factor(Signature, levels = c("Low", "High")),
                  y = Value, fill = Signature)) +
    geom_boxplot(outlier.size = 1, outlier.colour = "black",
                 width = 0.7, colour = "black", alpha = 0.8) +
    labs(title = paste0(hallmark, "- IRE"),
         x = "Signature",
         y = "Value Expression",
         fill = "") +
    scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#d62728")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
          axis.text.y = element_text(size = 12, face = "plain"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          panel.grid.major = element_line(size = 0.2, color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    stat_compare_means(method = "kruskal.test", label = "p.signif", label.x = 1.5, size = 5)
  ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Hallmark/IRE/",
                hallmark, "_boxplot.png"), plot = p, width = 8, height = 6)}




## TCGA
## Select only signature genes columns
df_tpm_tcga_cleaned_tumor_signature <- df_tpm_tcga_cleaned_tumor %>% 
  filter(rownames(df_tpm_tcga_cleaned_tumor) %in% genes)

## Load Hallmark
hallmark_pathways   <- getGmt(paste0(data_path_own, "/h.all.v2023.2.Hs.symbols.gmt"))
expression_data_tcga     <- as.matrix(df_tpm_tcga_cleaned_tumor_signature)
ssgsea_params_tcga       <- ssgseaParam(expression_data_tcga, hallmark_pathways)
ssGSEA_scores_tcga       <- gsva(ssgsea_params_tcga)
ssGSEA_scores_tcga       <- as.data.frame(t(ssGSEA_scores_tcga))
ssGSEA_scores_scale_tcga <- scale(ssGSEA_scores_tcga)
ssGSEA_scores_scale_tcga <- data.frame(ssGSEA_scores_scale_tcga)

## Wrangling
ssGSEA_scores_scale_tcga$ID <- rownames(ssGSEA_scores_scale_tcga)
## Merge
ssGSEA_scores_scale_final_tcga <- merge(ssGSEA_scores_scale_tcga, df_signature_tcga, by = "ID")
rownames(ssGSEA_scores_scale_final_tcga) <- ssGSEA_scores_scale_final_tcga$ID
ssGSEA_scores_scale_final_tcga$ID <- NULL
ssGSEA_scores_scale_final_tcga$Genes <- NULL 

## Reshape
df_long_tcga <- tidyr::pivot_longer(ssGSEA_scores_scale_final_tcga, cols = starts_with("HALLMARK"), 
                                    names_to = "Variable", values_to = "Value")
## Plot
hallmark_unique <- unique(df_long_tcga$Variable)
for (hallmark in hallmark_unique){
  df_subset <- subset(df_long_tcga, Variable == hallmark)
  p <- ggplot(df_subset,
              aes(x = factor(Signature, levels = c("Low", "High")),
                  y = Value, fill = Signature)) +
    geom_boxplot(outlier.size = 1, outlier.colour = "black",
                 width = 0.7, colour = "black", alpha = 0.8) +
    labs(title = paste0(hallmark, "- TCGA"),
         x = "Signature",
         y = "Value Expression",
         fill = "") +
    scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#d62728")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
          axis.text.y = element_text(size = 12, face = "plain"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          panel.grid.major = element_line(size = 0.2, color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    stat_compare_means(method = "kruskal.test", label = "p.signif", label.x = 1.5, size = 5)
  ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Hallmark/TCGA/",
                hallmark, "_boxplot.png"), plot = p, width = 8, height = 6)}


## Heatmap
## OWN
ssGSEA_scores_scale_final_signature <- ssGSEA_scores_scale_final$Signature
ssGSEA_scores_scale_final$Signature <- NULL
## Annotation
annotation_col <- data.frame(Group = ssGSEA_scores_scale_final_signature)
rownames(annotation_col) <- colnames(t(ssGSEA_scores_scale_final))
annotation_colors <- list(Group = c("Low" = "#1f77b4", "High" = "#d62728"))
## Plot
pheatmap(t(ssGSEA_scores_scale_final),
         annotation_col = annotation_col,  
         annotation_colors = annotation_colors,
         scale = "row",  
         fontsize = 12,  
         main = "Hallmark - IRE",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE)

## TCGA
ssGSEA_scores_scale_final_tcga_signature <- ssGSEA_scores_scale_final_tcga$Signature
ssGSEA_scores_scale_final_tcga$Signature <- NULL
## Annotation
annotation_col <- data.frame(Group = ssGSEA_scores_scale_final_tcga_signature)
rownames(annotation_col) <- colnames(t(ssGSEA_scores_scale_final_tcga))
annotation_colors <- list(Group = c("Low" = "#1f77b4", "High" = "#d62728"))
## Plot
pheatmap(t(ssGSEA_scores_scale_final_tcga),
         annotation_col = annotation_col,  
         annotation_colors = annotation_colors,
         scale = "row",  
         fontsize = 12,  
         main = "Hallmark - TCGA",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE)

## ***************************************
## PROGENY
## ***************************************
net <- decoupleR::get_progeny(organism = "human", 
                              top = 500)
























