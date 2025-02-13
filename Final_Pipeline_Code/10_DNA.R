## Load script
source("00_Support_ssGSEA.R")

## ***************************************
## LOAD DATA
## ***************************************
df_cna_genes_tcga       <- read.delim(paste0(data_path, "/CNA_Genes.txt"))
df_mutated_genes_tcga   <- read.delim(paste0(data_path, "/Mutated_genes.txt"))
df_fraction_genome_tcga <- read.delim(paste0(data_path, "/Fraction_Genome_Altered.txt"))

## Load alterations across samples for 100 top genes
df_alterations <- read_tsv(paste0(data_path, "/alterations_across_samples.tsv"),
                           show_col_types = FALSE)
## Load Signature
df_signature_tcga <- read_csv(paste0(data_path_own, "/TCGA_Signature.csv"))
names(df_signature_tcga)[names(df_signature_tcga) == "...1"] <- "Patient_ID"
names(df_alterations)[names(df_alterations) == "Patient ID"] <- "Patient_ID"

## Top 100 genes for percentages
df_cna_genes_tcga_sorted <- df_cna_genes_tcga %>% 
  arrange(desc(Freq))
df_cna_genes_tcga_sorted <- head(df_cna_genes_tcga_sorted, 100)

## Merge
df_final <- merge(df_signature_tcga, df_alterations, by = "Patient_ID")

## Extract genes to plot
cols_to_change <- grepl("[: ]", names(df_final))
names(df_final)[cols_to_change] <- gsub("[: ]+", "_", names(df_final)[cols_to_change])

##Select columns that end with _AMP
cols_to_extract <- grep("_AMP$", names(df_final), value = TRUE)
cols_to_extract <- c(cols_to_extract, "Signature")
## Extract the corresponding columns from the dataframe
df_extracted <- df_final[, cols_to_extract]

## Plot
ggplot(df_extracted, aes(x = CCM2L_AMP, y = after_stat(count), fill = Signature)) +
  geom_bar(position = "dodge") +
  labs(title = "", x = "RNU6-442P_AMP", y = "Count", fill = "") +
  scale_fill_manual(values = c("Low" = "#1f77b4", "High" = "#d62728")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(size = 12, face = "plain"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank())

