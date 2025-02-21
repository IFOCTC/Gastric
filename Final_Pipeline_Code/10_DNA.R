## Load script
source("00_Support_ssGSEA.R")
library(TCGAbiolinks)
data_path      <- "/Users/ifo/Desktop/IFO/Gastric/Data"
data_path_tcga <- "/Users/ifo/Desktop/IFO/Gastric/Data/TCGA"
data_path_mutations <- "/Users/ifo/Desktop/IFO/Gastric/Data/msk_met_2021"

## ***************************************
## LOAD DATA
## ***************************************
## MET
df_mutated_genes_met  <- read.delim(paste0(data_path, "/Mutated_genes.txt"))
df_cna_genes_met       <- read.delim(paste0(data_path, "/CNA_Genes.txt"))
df_mutations_met            <- read.delim(paste0(data_path_mutations, "/data_mutations.txt"),
                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
## TCGA
df_signature_tcga            <- read.csv(paste0(data_path_tcga, "/TCGA_Signature.csv"))
# query <- GDCquery(project = c("TCGA-STAD"), data.category = "Simple Nucleotide Variation", 
#                   access = "open", data.type = "Masked Somatic Mutation", 
#                   workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
# GDCdownload(query)
# df_maf_tcga <- GDCprepare(query)
# ## Write csv
# write.csv(df_maf_tcga, file = paste0(data_path_tcga, "/df_maf_tcga.csv"), row.names = TRUE)
## Load data
df_maf_tcga <- read.csv(paste0(data_path_tcga, "/df_maf_tcga.csv"))
df_maf_tcga$X <- NULL

## Wrangling
df_mutated_genes_met$Freq <- sort(df_mutated_genes_met$Freq, decreasing = TRUE)
df_cna_genes_met$Freq     <- sort(df_cna_genes_met$Freq, decreasing = TRUE)
names(df_signature_tcga)[names(df_signature_tcga) == "X"] <- "Patient_ID"

## Extract genes name
df_mutated_genes_met_sorted <- df_mutated_genes_met %>% 
  filter(Freq >= 2)
df_cna_genes_met_sorted     <- df_cna_genes_met %>% 
  filter(Freq >= 2)
genes_mutated <- df_mutated_genes_met_sorted$Gene
genes_cna     <- df_cna_genes_met_sorted$Gene
## Select only useful genes
genes <- c(genes_mutated, genes_cna)
genes <- setdiff(genes, "H3C7")

## Wrangling for tumor sample code
df_maf_tcga_filtered <- df_maf_tcga %>%
  mutate(Patient_ID = sapply(strsplit(as.character(df_maf_tcga$Tumor_Sample_Barcode), "-"),
                             function(x) paste(x[1:3], collapse = "-")))
df_maf_tcga <- df_maf_tcga %>%
  mutate(Patient_ID = sapply(strsplit(as.character(df_maf_tcga$Tumor_Sample_Barcode), "-"),
                             function(x) paste(x[1:3], collapse = "-")))
df_maf_tcga_filtered <- df_maf_tcga_filtered %>% 
  dplyr::filter(Hugo_Symbol %in% genes)

## Work with MAF TCGA file
df_maf_tcga_final <- df_maf_tcga %>% 
  ## Select relevant columns
  dplyr::select(Patient_ID, Hugo_Symbol, Variant_Classification) %>% 
  ## Add column for type of mutations
  group_by(Patient_ID, Hugo_Symbol) %>% 
  summarise(Mutations = paste(Variant_Classification, collapse = ","),
            .groups = "drop") %>% 
  ## Create wide matrix
  tidyr::pivot_wider(names_from = Hugo_Symbol, values_from = Mutations, values_fill = list(Mutations = ""))
## Set rownames
patient_id <- df_maf_tcga_final$Patient_ID
df_maf_tcga_final$Patient_ID <- NULL
rownames(df_maf_tcga_final)  <- patient_id

## Convert values
df_maf_tcga_final_binary <- df_maf_tcga_final %>% 
  mutate_all(~ ifelse(!is.na(.) & . != "", 1, 0))
rownames(df_maf_tcga_final_binary)  <- patient_id
df_maf_tcga_final_binary$Patient_ID <- patient_id

## Merge
df_final      <- merge(df_signature_tcga, df_maf_tcga_final_binary, by = "Patient_ID")
df_final_long <- df_final %>% 
  tidyr::pivot_longer(cols = -c(Patient_ID, Signature, Genes), names_to = "Gene", values_to = "Value")
df_final_long$Value <- factor(df_final_long$Value, levels = c(0, 1))

## Plot only mutated
df_final_long_1 <- df_final_long %>%
  filter(Value == 1)

## Plot only the values for 1 (mutations)
ggplot(df_final_long_1, aes(x = Gene, y = after_stat(count), fill = factor(Signature))) +
  geom_bar(position = "dodge") +  
  scale_fill_manual(values = c("Low" = "#3498DB", "High" = "#E74C3C")) + 
  labs(title = "Gene Mutations by Signature", x = "Gene",  y = "Count") +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "top") + 
  guides(fill = guide_legend(title = "Signatue Status")) + 
  theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12))

## Fisher Test
df_final_test            <- df_final
df_final_test$Patient_ID <- NULL
df_final_test$Genes      <- NULL

## Store the p-values for each gene
fisher_res <- data.frame(Gene = character(), P_value = numeric(), stringsAsFactors = FALSE)
## Loop over each gene (starting from the 2nd column)
for (gene in colnames(df_final_test[-1])) {
  ## Set reference level for Signature
  df_final_test$Signature <- relevel(factor(df_final_test$Signature), ref = "High")
  ## Create a contingency table for the current gene
  contingency_table <- table(df_final_test[[gene]], df_final_test$Signature)
  ## Check if the contingency table has at least 2 rows and 2 columns
  if (nrow(contingency_table) > 1 && ncol(contingency_table) > 1) {
    ## Perform Fisher's exact test
    fisher_test <- fisher.test(contingency_table)
    ## Store the gene name, p-value, and odds ratio (from the test result)
    fisher_res <- rbind(fisher_res, data.frame(Gene = gene, 
                                               P_value = fisher_test$p.value,
                                               NegLogPValue = -log10(fisher_test$p.value),
                                               Odds_Ratio = fisher_test$estimate))
    rownames(fisher_res) <- NULL
  } else {
    ## If the contingency table is not valid, print a message or skip
    message(paste("Skipping gene", gene, "due to insufficient variation"))
  }
}

## Optionally, you can filter the results to see only significant genes
significant_genes <- fisher_res[fisher_res$P_value <= 0.05, ]

## Add color specification
fisher_res$in_list <- ifelse(fisher_res$Gene %in% genes, "In List", "Not in List")
## Add Significance specification
fisher_res$Significance_P_value <- ifelse(fisher_res$P_value <= 0.05, "Significant", "Not Significance")
fisher_res$Significance_Odds_Ratio <- ifelse(fisher_res$Odds_Ratio >= 1, "OR >= 1", "OR < 1")
## Add Association specification
fisher_res <- fisher_res %>%
  mutate(Association = case_when(
    Significance_P_value == "Significant" & Significance_Odds_Ratio == "OR >= 1" ~ "Positive",
    Significance_P_value == "Significant" & Significance_Odds_Ratio == "OR < 1" ~ "Negative",
    TRUE ~ "N.S."))

## Cap Inf values
fisher_res$Odds_Ratio <- ifelse(is.infinite(fisher_res$Odds_Ratio), 20, fisher_res$Odds_Ratio)

## Plot
ggplot(fisher_res, aes(x = Odds_Ratio, y = NegLogPValue)) +
  geom_point(aes(color = factor(
    ifelse(P_value > 0.05, "N.S.", ifelse(Odds_Ratio > 1, "Positive", "Negative")), 
    levels = c("Positive", "Negative", "N.S."))), alpha = 0.7) +
  scale_color_manual(name = "Association",  values = c("Positive" = "#E74C3C",
                                                       "Negative" = "#3498DB",  
                                                       "N.S." = "gray")) +
  labs(x = "Odds Ratio", y = "-log10(P-value)", title = "Fisher's Test",
       subtitle = "Associations of genes with P-value threshold and Odds Ratios") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        axis.ticks.x = element_line()) +
  # geom_text(aes(label = ifelse(P_value <= 0.05, as.character(Gene), "")),
  #           vjust = -0.5, size = 3, check_overlap = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6),
                     labels = scales::comma_format())

## Extract Positive and Negative Association genes
positive_genes <- fisher_res %>%
  filter(P_value <= 0.05, Odds_Ratio >= 1) %>%
  pull(Gene)
negative_genes <- fisher_res %>%
  filter(P_value <= 0.05, Odds_Ratio < 1) %>%
  pull(Gene)
## Check
print(positive_genes);length(positive_genes)
print(negative_genes);length(negative_genes)



