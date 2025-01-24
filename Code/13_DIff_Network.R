## Code useful to perform differential co-expression analysis
## considering the set of genes up-regulated within the cellular lines

## Load Useful
library(readr)
library(readxl)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(affy)
library(scales)
library(pheatmap)
library(lsa)
library(vsn)
library(EnhancedVolcano)
library(glue)
library(biomaRt)
library(psych)
library(reshape2)
library(stringr)
library(network)
library(igraph)
library(sna)
library(GGally)
library(gplots)
library(gridExtra)
library(shiny)
library(survival)
library(survminer)
options(scipen = 999)

## Set data path directory
data_path   <- "C:/Users/david/Documents/IFO/Gastric/Data"
output_path <- "C:/Users/david/Documents/IFO/Gastric/Output_Res"

## Load Data
df_clinical <- read_excel(paste0(data_path, "/stomaci_wes_all.xlsx"))
df_tpm  <- read_tsv(paste0(data_path, "/salmon.merged.gene_tpm.tsv"), show_col_types = FALSE)
df_degs <- read.csv(paste0(output_path, "/df_volcano_degs_08.csv"))
## For the analysis, we consider the down-regulated genes, that are down-regulated in sensitive


## **********************************************************************
dim(df_tpm)
## TPM Filtering
row_num_t   <- 5
min_samples <- 50
df_tpm <- df_tpm[ rowSums(df_tpm >= row_num_t) >= min_samples, ]

## Clinical Filtering
columns_survival <- c("ID", "OS_time", "OS_Event", "PFS_time",
                      "PFS_event", "Coorte", "CT1L_1SI", "CT1L_FLOT_1SI")
## Consider only some columns and patients with treatment
df_clinical_survival <- df_clinical %>% 
  dplyr::select(all_of(columns_survival)) %>% 
  filter(CT1L_FLOT_1SI == 1) %>% 
  as.data.frame()

## Patients classification:
## Fast Progressor (PFS <= 6 months & PFS Event == 1), Long Responder (PFS >= 12 months),
## Conventional Responder (remaining part). 
df_clinical_survival <- df_clinical_survival %>% 
  mutate(Response = case_when(PFS_time <= 6 & PFS_event == 1 ~ "Fast Progressor",
                              PFS_time >= 12  ~ "Long Responder",
                              TRUE ~ "Conventional Responder"))
## Filtering
df_clinical_survival <- df_clinical_survival %>%
  filter(Response != "Conventional Responder")

## Check
View(df_clinical_survival)

## Plot
ggplot(df_clinical_survival, aes(x = Response)) +
  geom_bar(aes(fill = Response)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Response Distribution", x = "", y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("#E74C3C", "#F39C12", "#3498DB"))

## Wrangling TPM
df_tpm$gene_id <- NULL
genes_down_regulated <- df_degs %>% 
  filter(Expression == "Down-regulated")

## Filtering on genes considering DEGs Up downregulated and tumor patients
df_tpm_filtered <- df_tpm %>% 
  filter(df_tpm$gene_name %in% genes_down_regulated$X) %>% 
  dplyr::select(gene_name, matches("T$"))

## Transpose ---> Patient x Gene
gene_names <- df_tpm_filtered$gene_name
df_tpm_filtered$gene_name <- NULL
df_tpm_filtered <- data.frame(df_tpm_filtered)
rownames(df_tpm_filtered) <- make.names(gene_names, unique = TRUE)
## Transpose
df_tpm_filtered <- t(df_tpm_filtered)
df_tpm_filtered <- data.frame(df_tpm_filtered)
## Operations for ID
rownames(df_tpm_filtered) <- gsub("\\.1T$", "", rownames(df_tpm_filtered))
rownames(df_tpm_filtered) <- gsub("\\.", "-", rownames(df_tpm_filtered))

## Now we divide the tpm dataframe considering the two reference classes
## Long Responder - Fast Progressor

## Filter id patients w.r.t. the conditions
id_patients_fast_progressor <- df_clinical_survival %>% 
  filter(Response == "Fast Progressor")
id_patients_long_responder  <-df_clinical_survival %>% 
  filter(Response == "Long Responder")

## Filter TPM
df_tpm_filtered_fast_progressor <- df_tpm_filtered %>% 
  filter(rownames(df_tpm_filtered) %in% id_patients_fast_progressor$ID)
df_tpm_filtered_long_responder  <- df_tpm_filtered %>% 
  filter(rownames(df_tpm_filtered) %in% id_patients_long_responder$ID)


## **********************************************************************
## Differential Co-Expressed Network
## Fast Progressor
cor_fast_progressor <- cor(df_tpm_filtered_fast_progressor, method = "pearson")
diag(cor_fast_progressor) <- 0
## Long Responder
cor_long_responder <- cor(df_tpm_filtered_long_responder, method = "pearson")
diag(cor_long_responder) <- 0
## Calculation of differential correlations
## Apply Fisher z-transformation
z_fp <- 0.5 * log((1 + cor_fast_progressor) / (1 - cor_fast_progressor))
z_lr <- 0.5 * log((1 + cor_long_responder) / (1 - cor_long_responder))

## Sample size for each of the condition
n_fp <- ncol(cor_fast_progressor)
n_lr <- ncol(cor_long_responder)

## Z-score to evaluate the correlation
Z <- (z_fp - z_lr) / sqrt(1/(n_fp - 3) + (1/(n_lr - 3)))

## Threshold 
t <- 10

## Adjacency Matrix a_ij = 0, if |Z| < t
adj_differential <- data.frame(ifelse(abs(Z) < t, 0, 1))

## Generate network
net_diff_coex <- graph_from_adjacency_matrix(as.matrix(adj_differential),
                                             mode = "undirected",
                                             diag = FALSE)

## Analysis
## Compute the degree index
degree_diff_coex <- sort(rowSums(adj_differential, na.rm = TRUE),
                         decreasing = TRUE)

## Who is the most connected hub?
(hub_most_connected <- degree_diff_coex[1])

## Plot degree distribution to check if the graph is a scale-free
df5           <- data.frame(cbind(degree_diff_coex))
colnames(df5) <- "Degree"
ggplot(df5, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Differential Co-Expression Network)") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal()

## Compute hubs
## how big is the degree of the most connected nodes?
(q_diff_coex <- quantile(degree_diff_coex[degree_diff_coex > 0],
                         0.80))

## Find the hubs (5% of the nodes with highest degree values)
hubs_diff_coex <- degree_diff_coex[degree_diff_coex >= q_diff_coex]
## Genes
names(hubs_diff_coex)
## How many?
length(hubs_diff_coex)

## **********************************************************************
## Significance genes

## Add ID for tpm data
df_tpm_filtered$ID <- rownames(df_tpm_filtered)
df_tpm_filtered_hubs <- df_tpm_filtered[, colnames(df_tpm_filtered) %in%
                                          names(hubs_diff_coex)]
df_tpm_filtered_hubs$ID <- rownames(df_tpm_filtered_hubs)
## Joing with tpm and clinical features
df_final_survival <- merge(df_clinical_survival, df_tpm_filtered, by = "ID")
df_final_survival_filtered <- merge(df_clinical_survival, df_tpm_filtered_hubs, by = "ID")
## Select columns to exclude for testing
columns_to_exlude <- c("ID", "OS_time", "OS_Event", "PFS_time", "PFS_event",
                       "Coorte", "CT1L_1SI", "CT1L_FLOT_1SI")

## Filter
## Variable: Response
res <- df_final_survival_filtered %>% 
  dplyr::select(-all_of(columns_to_exlude)) %>% 
  kruskal_test_on_genes(condition_col = "Response") %>% 
  filter(p_value <= 0.1) 

## Boxplot
plots_boxplot <- plot_kruskal_boxplot(df_final_survival_filtered, "Response")

## Plot only significantly genes
for (gene in rownames(res)){
  print(plots_boxplot[[gene]])
}











