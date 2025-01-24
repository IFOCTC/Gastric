## Code useful to perform analysis at patients level

## Load scripts
source("00_Support.R")

## Load Data
df_tpm      <- read_tsv(paste0(output_path, "/tpm_cleaned.tsv"), show_col_types = FALSE)
df_clinical <- read_tsv(paste0(output_path, "/clinical_cleaned.tsv"), show_col_types = FALSE)
df_survival <- read_tsv(paste0(output_path, "/survival_cleaned.tsv"), show_col_types = FALSE)
## Convert to df
df_tpm      <- data.frame(df_tpm)
df_clinical <- data.frame(df_clinical)
df_survival <- data.frame(df_survival)

## RUN THIS CHUNK ONLY IF WE WANT TO CONSIDER PATIENTS WITH STAGE III-IV
## Consider only patients with TNMN or stage III-IV
df_clinical_filtered <- df_clinical %>%
  filter(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA",
                                                "Stage IIIB", "Stage IIIC",
                                                "Stage IV") )
## Filter survival considering this information
df_survival_filtered <- df_survival %>%
  filter(sample %in% df_clinical_filtered$sample)
## Filter tpm
df_tpm_filtered <- df_tpm %>%
  filter(ID %in% df_clinical_filtered$submitter_id)

## Compute mean for each patient (column)
df_tpm <- df_tpm %>%
  mutate(Mean = rowMeans(across(-ID), na.rm = TRUE))

## Compute 75% (Q3) and median value
q3_value     <- quantile(df_tpm$Mean, 0.75, na.rm = TRUE)
median_value <- median(df_tpm$Mean, na.rm = TRUE)

## Plot
ggplot(df_tpm, aes(x = Mean, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "black", bins = 30) +
  geom_vline(aes(xintercept = median_value), color = "steelblue",
             linetype = "dashed", linewidth = 1) +
  ggtitle("Mean Patients Distribution x Gene") +
  xlab("Mean") +
  ylab("Count") +
  theme_minimal()

## Thresholding
df_tpm <- df_tpm %>%
  mutate(Condition_Q3 = ifelse(Mean >= q3_value, "High", "Low"),
         Condition_Mean = ifelse(Mean >= median_value, "High", "Low"))
## Count
table(df_tpm$Condition_Mean)

## CLINICAL
## Merge
df_final_survival <- merge(df_survival, df_tpm, by = "ID")
## Create month
df_final_survival$OS.time.month <- df_final_survival$OS.time/30
df_final_survival <- df_final_survival %>% 
  mutate(OS_Status = case_when(OS.time.month <= 6 & OS == 1 ~ "Short Survival",
                               OS.time.month >= 24  ~ "Long Survival",
                               TRUE ~ "Conventional Survival"))

## Plot
ggplot(df_final_survival, aes(x = OS_Status)) +
  geom_bar(aes(fill = OS_Status)) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Overall Status Distribution", x = "", y = "Count", fill = "OS Status") +
  theme_minimal() +
  scale_fill_manual(values = c("#F39C12", "#3498DB", "#E74C3C"))

## Survival Analysis
## OS Survival
fit_os <- survfit(Surv(OS.time.month, OS) ~ OS_Status,
                  data = df_final_survival)
## Plot
ggsurv_os <- ggsurvplot(fit_os, data = df_final_survival,
                        pval = TRUE, conf.int = FALSE, palette = c("#F39C12", "#3498DB", "#E74C3C"),
                        xlim = c(0, 48), title = "Overall Survival Analysis",
                        xlab = "Months",
                        ylab = "OS Probability", break.time.by = 4,
                        ggtheme = theme_light(), risk.table = TRUE,
                        risk.table.y.text.col = TRUE,
                        risk.table.height = 0.25,
                        risk.table.y.text = TRUE,
                        conf.int.style = "step",
                        surv.median.line = "hv")
ggsurv_os

## TESTING GENES
## Select columns to exclude for testing
columns_to_exlude <- c("Mean", "Condition_Mean", "Condition_Q3", "OS.time.month")

## Filter - Exclude also the conventional responder
## Variable: Response
df_final_survival_filtered <- df_final_survival %>%
  dplyr::select(-all_of(columns_to_exlude))
#   filter(OS_Status != "Conventional Survival")

## Plot
plots_boxplot <- plot_kruskal_boxplot(df_final_survival_filtered, "OS_Status")
plots_boxplot
