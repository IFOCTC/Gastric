## Load script
source("00_Support.R")

## ***************************************
## LOAD DATA
## ***************************************
df_degs <- read.csv(paste0(output_path, "/df_volcano_degs_08.csv"))
df_tpm  <- read_tsv(paste0(data_path, "/salmon.merged.gene_tpm.tsv"), show_col_types = FALSE)
df_clinical <- read_excel(paste0(data_path, "/stomaci_wes_all.xlsx"))
ids_marcello <- c("IRE-012", "IRE-016", "IRE-017", "IRE-023", "IRE-034", "IRE-039",
                  "IRE-044", "IRE-048", "IRE-053", "IRE-056", "IRE-059", "IRE-060",
                  "IRE-061", "IRE-063", "IRE-067", "IRE-068", "IRE-071", "IRE-075",
                  "IRE-081", "IRE-083", "IRE-086", "IRE-090", "IRE-098", "IRE-100",
                  "IRE-104", "IRE-105", "IRE-107", "IRE-110", "IRE-111", "IRE-115",
                  "IRE-120", "IRE-121", "IRE-125", "IRE-128", "IRE-130", "IRE-134",
                  "IRE-135", "IRE-136", "IRE-137", "IRE-138", "IRE-139", "IRE-141",
                  "IRE-142", "IRE-143", "IRE-144", "IRE-145", "IRE-148", "IRE-149",
                  "IRE-151", "IRE-153", "IRE-154", "IRE-155", "IRE-156", "IRE-157",
                  "IRE-158", "IRE-159", "IRE-160", "IRE-161", "IRE-162", "IRE-163",
                  "IRE-164", "IRE-165", "IRE-166", "IRE-167", "IRE-168", "IRE-169",
                  "IRE-170", "IRE-171", "IRE-173", "IRE-174", "IRE-175", "IRE-176",
                  "IRE-177", "IRE-178", "IRE-179", "IRE-180", "IRE-181", "IRE-182",
                  "IRE-183", "IRE-184", "IRE-185", "IRE-186", "IRE-187", "IRE-188",
                  "IRE-189", "IRE-190", "IRE-191", "IRE-192", "IRE-193", "IRE-194",
                  "IRE-195", "IRE-196", "IRE-197", "IRE-198", "IRE-199", "IRE-201",
                  "IRE-202", "IRE-203", "IRE-204", "IRE-205", "IRE-206", "IRE-207",
                  "IRE-208", "IRE-209", "IRE-210", "IRE-211", "IRE-213", "IRE-214",
                  "IRE-215", "IRE-216", "IRE-217", "IRE-219", "IRE-220", "IRE-221",
                  "IRE-222", "IRE-224", "IRE-225", "IRE-227", "IRE-228", "IRE-229",
                  "IRE-230", "IRE-231", "IRE-234", "IRE-235", "IRE-236", "IRE-237",
                  "IRE-238", "IRE-239", "IRE-242", "IRE-244", "IRE-246", "IRE-247",
                  "IRE-248", "IRE-249", "IRE-251", "IRE-253", "IRE-254", "IRE-255",
                  "IRE-256", "IRE-258", "IRE-260", "IRE-261", "IRE-262", "IRE-263",
                  "IRE-264", "IRE-265", "IRE-267", "IRE-268", "IRE-269", "IRE-270",
                  "IRE-271", "IRE-272", "IRE-276", "IRE-277", "IRE-280", "IRE-281",
                  "IRE-282", "IRE-284", "IRE-285", "IRE-287", "IRE-288", "IRE-289",
                  "IRE-290", "IRE-291", "IRE-292", "IRE-293", "IRE-294", "IRE-295",
                  "IRE-296", "IRE-297", "IRE-298", "IRE-299", "IRE-300", "IRE-301",
                  "IRE-302")

## ***************************************
## WRANGLING
## ***************************************
## Select only tumor patients in TPM and Clinical
ids_tpm <- df_tpm %>% 
  dplyr::select(matches("T$"))
ids_tpm <- colnames(ids_tpm)
ids_tpm <- cbind(ids_tpm)
colnames(ids_tpm) <- "ID"
ids_tpm[,1] <- gsub("\\.1T$", "", ids_tpm[,1])
ids_tpm[,1] <- gsub("\\.", "-", ids_tpm[,1])
## Check ID patients
cat("IDs Marcello: ", length(ids_marcello))
cat("IDs Clinical: ", length(df_clinical$ID))
cat("IDs Tpm: ", length(ids_tpm))
## IDs missing comparing marcello's patients and df final survival id
ids_missing <- c("IRE-056", "IRE-104", "IRE-105",
                 "IRE-107", "IRE-136", "IRE-153",
                 "IRE-157", "IRE-202", "IRE-291")

## Differences Gens
## EHR2

## CLDN18


## ***************************************
## EXPLORATORY ANALYSIS - CLINICAL FEATURES
## ***************************************
## Clinical Filtering
columns_survival <- c("ID", "OS_time", "OS_Event", "PFS_time",
                      "PFS_event", "Coorte", "CT1L_1SI", "CT1L_FLOT_1SI")
## Consider only some columns and patients with treatment
df_clinical_survival <- df_clinical %>% 
  dplyr::select(all_of(columns_survival))
## Patients classification:
## Fast Progressor (PFS <= 6 months & PFS Event == 1), Long Responder (PFS >= 12 months),
## Conventional Responder (remaining part). 
df_clinical_survival <- df_clinical_survival %>% 
  mutate(Response = case_when(PFS_time <= 6 & PFS_event == 1 ~ "Fast Progressor",
                              PFS_time >= 12  ~ "Long Responder",
                              TRUE ~ "Conventional Responder"),
         OS_Status = case_when(OS_time <= 6 & OS_Event == 1 ~ "Short Survival",
                               OS_time >= 24  ~ "Long Survival",
                               TRUE ~ "Conventional Survival"))
## This is create to make survival analysis above whole population
df_clinical_filtered <- df_clinical %>% 
  dplyr::select(all_of(columns_survival))
## For survival analysis
df_clinical_filtered <- df_clinical_survival %>% 
  mutate(Response = case_when(PFS_time <= 6 & PFS_event == 1 ~ "Fast Progressor",
                              PFS_time >= 12  ~ "Long Responder",
                              TRUE ~ "Conventional Responder"),
         OS_Status = case_when(OS_time <= 6 & OS_Event == 1 ~ "Short Survival",
                               OS_time >= 24  ~ "Long Survival",
                               TRUE ~ "Conventional Survival"))
df_clinical_filtered_2_classes <- df_clinical_survival %>% 
  mutate(Response = case_when(PFS_time <= 6 & PFS_event == 1 ~ "Fast Progressor",
                              PFS_time >= 12  ~ "Long Responder"),
         OS_Status = case_when(OS_time <= 6 & OS_Event == 1 ~ "Short Survival",
                               OS_time >= 24  ~ "Long Survival"))
## Plot whole population
ggplot(df_clinical_filtered, aes(x = OS_Status)) +
  geom_bar(aes(fill = OS_Status)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Overall Status Distribution", x = "", y = "Count", fill = "OS Status") +
  theme_minimal() +
  scale_fill_manual(values = c("#F39C12", "#3498DB", "#E74C3C"))

## ***************************************
## WRANGLING TPM
## ***************************************
df_tpm$gene_id <- NULL
genes_down_regulated <- df_degs %>% 
  filter(Expression == "Down-regulated")
## Filtering on genes considering DEGs Up downregulated and tumor patients
df_tpm_filtered <- df_tpm %>% 
  filter(df_tpm$gene_name %in% genes_down_regulated$X) %>% 
  dplyr::select(gene_name, matches("T$"))
## TPM Filtering
row_num_t   <- 30
min_samples <- 41
df_tpm_filtered <- df_tpm_filtered[rowSums(df_tpm_filtered >= row_num_t) >= min_samples, ]
## Transpose ---> Patient x Gene
gene_names <- df_tpm_filtered$gene_name
df_tpm_filtered$gene_name <- NULL
df_tpm_filtered <- data.frame(df_tpm_filtered)
rownames(df_tpm_filtered) <- make.names(gene_names, unique = TRUE)
## Transpose
df_tpm_filtered <- t(df_tpm_filtered)
df_tpm_filtered <- data.frame(df_tpm_filtered)
dim(df_tpm_filtered)

## ***************************************
## NETWORK APPROACH ANALYSIS
## ***************************************
## 05 - Define a function that computes the hubs of the network
compute_adjacency <- function(data, cor_type = NULL){
  ## Input:  Expression Data
  ## Output: Adjacency matrix
  
  ## Correlation matrix 
  cor_mat       <- cor(t(data), method = cor_type)
  diag(cor_mat) <- 0
  ## Correlation matrix (p-value)
  cor_padj <- corr.p(cor_mat, nrow(cor_mat),
                     adjust = "fdr", ci = FALSE)$p
  cor_padj[lower.tri(cor_padj)] <- t(cor_padj)[lower.tri(cor_padj)]
  ## Build adjacency matrix
  adj_mat_1 <- ifelse(cor_mat >= 0.3, 1,
                      ifelse(cor_mat <= -0.3, -1, 0))
  adj_mat_2 <- ifelse(abs(cor_padj) > 0.05, 0, 1) 
  adj_mat   <- adj_mat_1 * adj_mat_2
  ## Return
  return(adj_mat)
}

## 06 - Define a function that returns the list of the hubs
compute_hubs <- function(adj_matrix, mart = mart_names, quant){
  ## Input: Adjacency Matrix, mart
  ## Output: List of hubs, Degree, Level of quantile
  
  ## Compute degree
  degree <- sort(rowSums(abs(adj_matrix)), decreasing = TRUE)
  ## Compute quantile
  q <- quantile(degree[degree > 0], quant)
  ## Find the hubs (1-q% of the nodes with highest degree values)
  hubs <- degree[degree >= q]
  ## Let's order them by degree
  hubs_ord <- sort(hubs, decreasing = TRUE)
  ## Check names
  final_hubs <- hubs_ord
  ## Return
  return(list("code_hubs" = hubs,
              "final_hubs" = final_hubs,
              "degree" = degree,
              "q" = q))}

## Build Adjacency Matrix
adj_matrix <- compute_adjacency(t(df_tpm_filtered),
                                cor_type = "pearson")
## Compute hubs
hubs <- compute_hubs(adj_matrix, quant = 0.5)
## Degree distribution in order to check if the network is scale free
df1 <- data.frame(cbind(hubs$degree))
colnames(df1) <- "Degree"
(hist_plot <- ggplot(df1, aes(x = Degree)) +
  geom_histogram(fill = "steelblue", color = "black", alpha = 0.7, bins = 20) +  # Added border for bars
  ggtitle("Degree Distribution") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal(base_size = 14) +  # Base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  ## Center title and increase size
    axis.title.x = element_text(face = "bold", size = 14),             ## Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 14),             ## Bold y-axis title
    axis.text = element_text(size = 12),                               ## Increase size of axis text
    panel.grid.major = element_line(color = "gray", size = 0.5),       ## Adjust grid color/size
    panel.grid.minor = element_blank()                                 ## Remove minor grid lines for a cleaner look
  ))
## Compute Network
net <- network(adj_matrix, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "weights",
               directed = FALSE)
## Compute density
network.density(net)
## Giant Component
nrow(component.largest(net, result = "graph"))
## How many positive/negative correlations?
sum(adj_matrix == 1)
sum(adj_matrix == -1)
sum(adj_matrix == 0)
## Plot Hubs
plot_graph(net, hubs$code_hubs,
           title = "Differential Co-Expression Network Hubs")
## Identify the hubs
length(hubs$code_hubs)
hubs$code_hubs
names(hubs$code_hubs)











