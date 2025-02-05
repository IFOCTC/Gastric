## Load script
source("00_Support.R")

## ***************************************
## LOAD DATA
## ***************************************
df_rna_counts_tripl <- read_tsv(paste0(data_path,
                            "/salmon.merged.gene_counts_scaled_all.tsv"),
                            show_col_types = FALSE)
df_line_gastric     <- read_excel(paste0(data_path,
                            "/conte_TPM_linnegastrici_notriplkicatoR1.xlsx"))
missed              <- read_tsv(paste0(data_path,
                            "/salmon.merged.gene_counts_length_scaled_02.tsv"),
                            show_col_types = FALSE)

## ***************************************
## WRANGLING
## ***************************************
## Merge
df_rna_counts <- merge(df_rna_counts_tripl, missed, by = c("gene_id", "gene_name"))
## Set replications to delete
replication_to_delete <- c("UNITO.GTR.0125.R", "UNITO.GTR.0210.R", "UNITO.GTR.0221.R",
                           "UNITO.GTR.02459.R", "UNITO.GTR.042.R", "UNITO.GTR.0498.R",
                           "UNITO.GTR.0508.R", "UNITO.GTR.0607.R")
## Removing
df_rna_counts <- df_rna_counts %>% 
  dplyr::select(-all_of(replication_to_delete))
df_rna_counts <- data.frame(df_rna_counts)
## Check
rownames(df_rna_counts) <- make.names(df_rna_counts$gene_name,
                                      unique = TRUE)
## Remove some columns
gene_ids   <- df_rna_counts$gene_id
gene_names <- df_rna_counts$gene_name
df_rna_counts$gene_id   <- NULL
df_rna_counts$gene_name <- NULL
## Check 
dim(df_rna_counts)

## Filtering
row_num_t   <- 5
min_samples <- 10
df_rna_counts_filtered <- df_rna_counts[rowSums(df_rna_counts
                                                >= row_num_t) >= min_samples, ]
dim(df_rna_counts_filtered)

## Load Patients Condition
## Save patient ids and categories
##  RESISTANT: GTR 498, GTR 508, GTR 459; SENSITIVE: GTR210, GTR0042, GTR221, GTR 607, GTR 125
id_cellular_lines <- colnames(df_rna_counts_filtered)
condition_id      <- c("SENSITIVE", "SENSITIVE", "SENSITIVE", "RESISTANT",
                       "RESISTANT", "SENSITIVE", "SENSITIVE", "RESISTANT",
                       "RESISTANT", "RESISTANT", "RESISTANT", "SENSITIVE",
                       "SENSITIVE", "SENSITIVE", "SENSITIVE", "SENSITIVE")

## Create the dataframe
df_condition <- data.frame(Id = id_cellular_lines,
                           Condition = factor(condition_id))
df_condition <- process_df(df_condition, "Id")
## Add a column to the condition according to it
df_condition$Color[df_condition$Condition == "RESISTANT"] <- c("red3")
df_condition$Color[df_condition$Condition == "SENSITIVE"] <- c("blue3")
## Wrangling for TPM
df_line_gastric <- process_df(df_line_gastric, "Gene")
df_line_gastric_plot <- t(df_line_gastric)
df_line_gastric_plot <- data.frame(df_line_gastric_plot)
df_line_gastric_plot <- df_line_gastric_plot %>% 
  dplyr::mutate(ID = rownames(df_line_gastric_plot))
df_condition_plot <- df_condition
df_condition_plot$ID <- rownames(df_condition_plot)
df_final_plot <- merge(df_condition_plot, df_line_gastric_plot,
                       by = "ID")
## GENE
ggplot(df_final_plot,
       aes(x = Condition,
           y = log2(as.numeric(SCD)+1), fill = Condition)) +
  geom_boxplot() +  
  labs(title = "SCD", subtitle = "IRE Cohort Analysis",
       x = "Condition",
       y = "Gene Expression",
       fill = "") +
  scale_fill_manual(values = c("SENSITIVE" = "#1f77b4",  
                               "RESISTANT" = "#d62728")) +
  ylim(0, max(log2(as.numeric(df_final_plot$SCD)+1))) +
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
  stat_compare_means()

## Filtering
row_num_t   <- 5
min_samples <- 9
df_line_gastric_filtered <- df_line_gastric[rowSums(df_line_gastric >= row_num_t) >= min_samples, ]
dim(df_line_gastric_filtered)

## ***************************************
## ANALYSIS - DEA
## ***************************************
dds <- DESeqDataSetFromMatrix(countData = round(df_rna_counts_filtered),
                              colData = df_condition,
                              design = ~ Condition)
## Set reference condition
dds$Condition <- relevel(dds$Condition, ref = "RESISTANT")
## DEA
dds <- DESeq(dds)
## Alpha 0.1
res_1 <- results(dds, alpha = 0.1) ## default
## How many adjusted p-values were less than 0.1
sum(res_1$padj <= 0.1, na.rm = T)
## Check
summary(res_1) 
## Alpha 0.05
res_05 <- results(dds, alpha = 0.05)
## How many adjusted p-values were less than 0.05
sum(res_05$padj <= 0.05, na.rm = T)
## Check
summary(res_05)
## Order
res_ordered_1  <- res_1[order(res_1$pvalue),]
res_ordered_05 <- res_05[order(res_05$pvalue),]

## Different thresholds
t_values <- c(0.8, 1, 1.5, 2)
result_08 <- generate_volcano(t = t_values[1], res = res_1, padj_t = 0.05, top_n = 100)
result_1  <- generate_volcano(t = t_values[2], res = res_1, padj_t = 0.05, top_n = 50)
result_15 <- generate_volcano(t = t_values[3], res = res_1, padj_t = 0.05, top_n = 50)
result_2  <- generate_volcano(t = t_values[4], res = res_1, padj_t = 0.05, top_n = 50)
## Access the plot, volcano data, DEGs data, and expression summary for a specific threshold
## 0.8
volcano_plot_t1 <- result_08$plot
volcano_data_t1 <- result_08$volcano_data
volcano_degs_t1 <- result_08$degs_data
expression_summary_t1 <- result_08$expression_table
print(expression_summary_t1)
## 1
volcano_plot_t2 <- result_1$plot
volcano_data_t2 <- result_1$volcano_data
volcano_degs_t2 <- result_1$degs_data
expression_summary_t2 <- result_1$expression_table
print(expression_summary_t2)
## 1.5
volcano_plot_t3 <- result_15$plot
volcano_data_t3 <- result_15$volcano_data
volcano_degs_t3 <- result_15$degs_data
expression_summary_t3 <- result_15$expression_table
print(expression_summary_t3)
## 2
volcano_plot_t4 <- result_2$plot
volcano_data_t4 <- result_2$volcano_data
volcano_degs_t4 <- result_2$degs_data
expression_summary_t4 <- result_2$expression_table
print(expression_summary_t4)

## ***************************************
## SAVE OUTPUT
## ***************************************

# ## Extract and save only up-regulated
# up_resistant_08 <- result_08$degs_data %>% filter(Expression == "Down-regulated")
# up_resistant_08_genes <- rownames(up_resistant_08)
# up_resistant_1  <- result_1$degs_data  %>% filter(Expression == "Down-regulated")
# up_resistant_1_genes <- rownames(up_resistant_1)
# up_resistant_15 <- result_15$degs_data %>% filter(Expression == "Down-regulated")
# up_resistant_15_genes <- rownames(up_resistant_15)
# up_resistant_2  <- result_2$degs_data  %>% filter(Expression == "Down-regulated")
# up_resistant_2_genes <- rownames(up_resistant_2)
# ## Save in 4 txt files
# writeLines(up_resistant_08_genes, paste0(output_path, "/up_resistant_08_genes.txt"))
# writeLines(up_resistant_1_genes, paste0(output_path, "/up_resistant_1_genes.txt"))
# writeLines(up_resistant_15_genes, paste0(output_path, "/up_resistant_15_genes.txt"))
# writeLines(up_resistant_2_genes, paste0(output_path, "/up_resistant_2_genes.txt"))