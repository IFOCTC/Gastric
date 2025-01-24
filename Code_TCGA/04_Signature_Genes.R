## Code to create a signature starting from significant and correct oriented genes

## Consider the survival df
dim(df_final_survival);View(df_final_survival)

## Top Genes:
## Filtering
df_final_survival_signature <- df_final_survival
col_names <- colnames(df_final_survival_signature)
## Check
dim(df_final_survival_signature);View(df_final_survival_signature)
## Setting rownames
rownames(df_final_survival_signature) <- df_final_survival_signature$ID
## Compute metrics for signature
df_final_survival_signature <- df_final_survival_signature[,top_genes_02] %>%
  mutate(Mean = rowMeans(across(everything()), na.rm = TRUE),  
         Median = apply(across(everything()), 1, median, na.rm = TRUE), 
         Q1 = apply(across(everything()), 1, quantile, probs = 0.25, na.rm = TRUE),  
         Q3 = apply(across(everything()), 1, quantile, probs = 0.75, na.rm = TRUE))
## Compute 75% (Q3) on mean dataset
q3_value     <- quantile(df_final_survival_signature$Mean, 0.75, na.rm = TRUE)
q1_value     <- quantile(df_final_survival_signature$Mean, 0.25, na.rm = TRUE)
median_value <- median(df_final_survival_signature$Mean)
## Thresholding
df_final_survival_signature <- df_final_survival_signature %>%
  mutate(Condition_Q3 = ifelse(Mean >= q3_value, "High", "Low"),
         Condition_Q1 = ifelse(Mean >= q1_value, "High", "Low"),
         Condition_Mean = ifelse(Mean >= median_value, "High", "Low"))
## Attach ID column
df_final_survival_signature$ID <- rownames(df_final_survival_signature)
## Merge
df_final_survival_signature <- merge(df_final_survival_signature, df_survival,
                                     by = "ID")
## Filtering operation for survival analysis
columns_survival <- c("ID", "Condition_Mean", "OS.time", "OS")
columns_survival <- c(columns_survival, top_genes)
df_final_survival_signature <- df_final_survival_signature %>% 
  dplyr::select(all_of(columns_survival))
## Create OS Time Month
df_final_survival_signature$OS.time.month <- df_final_survival_signature$OS.time/30

## Plot
## OS
fit_os <- survfit(Surv(OS.time.month, as.integer(OS)) ~ Condition_Mean,
                  data = df_final_survival_signature)
os_plot <- ggsurvplot(fit_os, data = df_final_survival_signature,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      xlim = c(0, 48),
                      title = "Overall Survival Analysis",
                      subtitle = "Top Genes (n: 12)",
                      xlab = "Time in Months",  ylab = "OS Probability",
                      break.time.by = 4, ggtheme = theme_light(), 
                      risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                      risk.table.y.text = TRUE, conf.int.style = "step",
                      surv.median.line = "hv")
os_plot
