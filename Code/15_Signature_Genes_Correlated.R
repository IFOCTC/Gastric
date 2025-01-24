## Code to create a signature starting from significant and correct oriented genes

## Top Genes: 30
## Filtering
df_final_survival_signature_corr <- df_final_survival %>% 
  dplyr::select(all_of(columns_to_mantain_corr))
## Check
dim(df_final_survival_signature_corr);View(df_final_survival_signature_corr)
## Setting rownames
rownames(df_final_survival_signature_corr) <- df_final_survival_signature_corr$ID
## Compute metrics for signature
df_final_survival_signature_corr <- df_final_survival_signature_corr[,top_genes_corr] %>%
  mutate(Mean = rowMeans(across(everything()), na.rm = TRUE),  
         Median = apply(across(everything()), 1, median, na.rm = TRUE), 
         Q1 = apply(across(everything()), 1, quantile, probs = 0.25, na.rm = TRUE),  
         Q3 = apply(across(everything()), 1, quantile, probs = 0.75, na.rm = TRUE))
## Compute 75% (Q3) on mean dataset
q3_value     <- quantile(df_final_survival_signature_corr$Mean, 0.75, na.rm = TRUE)
q1_value     <- quantile(df_final_survival_signature_corr$Mean, 0.25, na.rm = TRUE)
median_value <- median(df_final_survival_signature_corr$Mean)
## Thresholding
df_final_survival_signature_corr <- df_final_survival_signature_corr %>%
  mutate(Condition_Q3 = ifelse(Mean >= q3_value, "High", "Low"),
         Condition_Q1 = ifelse(Mean >= q1_value, "High", "Low"),
         Condition_Mean = ifelse(Mean >= median_value, "High", "Low"))
## Attach ID column
df_final_survival_signature_corr$ID <- rownames(df_final_survival_signature_corr)
## Merge
df_final_survival_signature_corr <- merge(df_final_survival_signature_corr, df_clinical,
                                     by = "ID")
## Filtering operation for survival analysis
columns_survival <- c("ID", "CD44", "LGALS1", "CTSZ", "PDLIM1",
                      "GNAI2", "CD68", "IFITM3", "RSU1", "RAB8B",
                      "TAF10", "SMAD3", "Condition_Mean", "OS_time", "OS_Event",
                      "PFS_time", "PFS_event", "CT1L_FLOT_1SI")
df_final_survival_signature_corr <- df_final_survival_signature_corr %>%
  dplyr::select(all_of(columns_survival))
## Check
dim(df_final_survival_signature_corr);View(df_final_survival_signature_corr)

## Filtering
df_final_survival_signature_filtered_CT1L_FLOT_corr <- df_final_survival_signature_corr %>% 
  dplyr::select(all_of(colnames(df_final_survival_signature_corr))) %>% 
  filter(df_final_survival_signature_corr$CT1L_FLOT_1SI == 1)
## Check
dim(df_final_survival_signature_filtered_CT1L_FLOT_corr);View(df_final_survival_signature_filtered_CT1L_FLOT_corr)

## Plot
## OS
fit_os <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
                  data = df_final_survival_signature_corr)
os_plot <- ggsurvplot(fit_os, data = df_final_survival_signature_corr,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      xlim = c(0, 48),
                      title = "Overall Survival Analysis",
                      subtitle = "Top Genes correlation >= 0.6 (n: 11)",
                      xlab = "Time in Months",  ylab = "OS Probability",
                      break.time.by = 4, ggtheme = theme_light(),
                      risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                      risk.table.y.text = TRUE, conf.int.style = "step",
                      surv.median.line = "hv")
os_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_os_filtered <- survfit(Surv(Survival_Time, as.integer(OS_Event)) ~ Condition_Mean,
                           data = df_final_survival_signature_filtered_CT1L_FLOT_corr)
os_plot_filtered <- ggsurvplot(fit_os_filtered, data = df_final_survival_signature_filtered_CT1L_FLOT_corr,
                               risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                               palette = c("#E7B800", "#2E9FDF"),
                               xlim = c(0, 48),
                               title = "Overall Survival Analysis",
                               subtitle = "Top Genes correlation >= 0.6 (n: 11), with CT1L_FLOT_1SI == 1",
                               xlab = "Time in Months",  ylab = "OS Probability",
                               break.time.by = 4, ggtheme = theme_light(),
                               risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                               risk.table.y.text = TRUE, conf.int.style = "step",
                               surv.median.line = "hv")
os_plot_filtered

## PFS
fit_pfs <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                   data = df_final_survival_signature_corr)
pfs_plot <- ggsurvplot(fit_pfs, data = df_final_survival_signature_corr,
                       risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                       palette = c("#E7B800", "#2E9FDF"),
                       xlim = c(0, 48),
                       title = "Progression Free Survival Analysis",
                       subtitle = "Top Genes correlation >= 0.6 (n: 11)",
                       xlab = "Time in Months",  ylab = "PFS Probability",
                       break.time.by = 4, ggtheme = theme_light(),
                       risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                       risk.table.y.text = TRUE, conf.int.style = "step",
                       surv.median.line = "hv")
pfs_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_pfs_filtered <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                            data = df_final_survival_signature_filtered_CT1L_FLOT_corr)
pfs_plot_filtered <- ggsurvplot(fit_pfs_filtered, data = df_final_survival_signature_filtered_CT1L_FLOT_corr,
                                risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                                palette = c("#E7B800", "#2E9FDF"),
                                xlim = c(0, 48),
                                title = "Progression Free Survival Analysis",
                                subtitle = "Top Genes correlation >= 0.6 (n: 11), with CT1L_FLOT_1SI == 1",
                                xlab = "Time in Months",  ylab = "PFS Probability",
                                break.time.by = 4, ggtheme = theme_light(),
                                risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                                risk.table.y.text = TRUE, conf.int.style = "step",
                                surv.median.line = "hv")
pfs_plot_filtered


## Top Genes: 16
## Filtering
df_final_survival_signature_corr_02 <- df_final_survival %>% 
  dplyr::select(all_of(columns_to_mantain_corr_02))
## Check
dim(df_final_survival_signature_corr_02);View(df_final_survival_signature_corr_02)
## Setting rownames
rownames(df_final_survival_signature_corr_02) <- df_final_survival_signature_corr_02$ID
## Compute metrics for signature
df_final_survival_signature_corr_02 <- df_final_survival_signature_corr_02[,top_genes_corr_02] %>%
  mutate(Mean = rowMeans(across(everything()), na.rm = TRUE),  
         Median = apply(across(everything()), 1, median, na.rm = TRUE), 
         Q1 = apply(across(everything()), 1, quantile, probs = 0.25, na.rm = TRUE),  
         Q3 = apply(across(everything()), 1, quantile, probs = 0.75, na.rm = TRUE))
## Compute 75% (Q3) on mean dataset
q3_value     <- quantile(df_final_survival_signature_corr_02$Mean, 0.75, na.rm = TRUE)
q1_value     <- quantile(df_final_survival_signature_corr_02$Mean, 0.25, na.rm = TRUE)
median_value <- median(df_final_survival_signature_corr_02$Mean)
## Thresholding
df_final_survival_signature_corr_02 <- df_final_survival_signature_corr_02 %>%
  mutate(Condition_Q3 = ifelse(Mean >= q3_value, "High", "Low"),
         Condition_Q1 = ifelse(Mean >= q1_value, "High", "Low"),
         Condition_Mean = ifelse(Mean >= median_value, "High", "Low"))
## Attach ID column
df_final_survival_signature_corr_02$ID <- rownames(df_final_survival_signature_corr_02)
## Merge
df_final_survival_signature_corr_02 <- merge(df_final_survival_signature_corr_02, df_clinical,
                                          by = "ID")
## Filtering operation for survival analysis
columns_survival <- c("ID", "CD44", "CTSZ", "GNAI2", "RSU1", "RAB8B",
                      "TAF10", "SMAD3", "Condition_Mean", "OS_time", "OS_Event",
                      "PFS_time", "PFS_event", "CT1L_FLOT_1SI")
df_final_survival_signature_corr_02 <- df_final_survival_signature_corr_02 %>%
  dplyr::select(all_of(columns_survival))
## Check
dim(df_final_survival_signature_corr_02);View(df_final_survival_signature_corr_02)

## Filtering
df_final_survival_signature_filtered_02_CT1L_FLOT_corr <- df_final_survival_signature_corr_02 %>% 
  dplyr::select(all_of(colnames(df_final_survival_signature_corr_02))) %>% 
  filter(df_final_survival_signature_corr_02$CT1L_FLOT_1SI == 1)
## Check
dim(df_final_survival_signature_filtered_02_CT1L_FLOT_corr);View(df_final_survival_signature_filtered_02_CT1L_FLOT_corr)

## Plot
## OS
fit_os <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
                  data = df_final_survival_signature_corr_02)
os_plot <- ggsurvplot(fit_os, data = df_final_survival_signature_corr_02,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      xlim = c(0, 48),
                      title = "Overall Survival Analysis",
                      subtitle = "Top Genes correlation >= 0.6 (n: 7)",
                      xlab = "Time in Months",  ylab = "OS Probability",
                      break.time.by = 4, ggtheme = theme_light(),
                      risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                      risk.table.y.text = TRUE, conf.int.style = "step",
                      surv.median.line = "hv")
os_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_os_filtered <- survfit(Surv(Survival_Time, as.integer(OS_Event)) ~ Condition_Mean,
                           data = df_final_survival_signature_filtered_02_CT1L_FLOT_corr)
os_plot_filtered <- ggsurvplot(fit_os_filtered, data = df_final_survival_signature_filtered_02_CT1L_FLOT_corr,
                               risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                               palette = c("#E7B800", "#2E9FDF"),
                               xlim = c(0, 48),
                               title = "Overall Survival Analysis",
                               subtitle = "Top Genes correlation >= 0.6 (n: 11), with CT1L_FLOT_1SI == 1",
                               xlab = "Time in Months",  ylab = "OS Probability",
                               break.time.by = 4, ggtheme = theme_light(),
                               risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                               risk.table.y.text = TRUE, conf.int.style = "step",
                               surv.median.line = "hv")
os_plot_filtered

## PFS
fit_pfs <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                   data = df_final_survival_signature_corr_02)
pfs_plot <- ggsurvplot(fit_pfs, data = df_final_survival_signature_corr_02,
                       risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                       palette = c("#E7B800", "#2E9FDF"),
                       xlim = c(0, 48),
                       title = "Progression Free Survival Analysis",
                       subtitle = "Top Genes correlation >= 0.6 (n: 7)",
                       xlab = "Time in Months",  ylab = "PFS Probability",
                       break.time.by = 4, ggtheme = theme_light(),
                       risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                       risk.table.y.text = TRUE, conf.int.style = "step",
                       surv.median.line = "hv")
pfs_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_pfs_filtered <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                            data = df_final_survival_signature_filtered_02_CT1L_FLOT_corr)
pfs_plot_filtered <- ggsurvplot(fit_pfs_filtered, data = df_final_survival_signature_filtered_02_CT1L_FLOT_corr,
                                risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                                palette = c("#E7B800", "#2E9FDF"),
                                xlim = c(0, 48),
                                title = "Progression Free Survival Analysis",
                                subtitle = "Top Genes correlation >= 0.6 (n: 7), with CT1L_FLOT_1SI == 1",
                                xlab = "Time in Months",  ylab = "PFS Probability",
                                break.time.by = 4, ggtheme = theme_light(),
                                risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                                risk.table.y.text = TRUE, conf.int.style = "step",
                                surv.median.line = "hv")
pfs_plot_filtered
