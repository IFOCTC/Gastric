## Code to create a signature starting from significant and correct oriented genes

## Condier the survival df
dim(df_final_survival);View(df_final_survival)

## Top Genes: 30
## Filtering
df_final_survival_signature <- df_final_survival %>% 
  dplyr::select(all_of(columns_to_mantain))
col_names <- colnames(df_final_survival_signature)
## Check
dim(df_final_survival_signature);View(df_final_survival_signature)
## Setting rownames
rownames(df_final_survival_signature) <- df_final_survival_signature$ID
## Compute metrics for signature
df_final_survival_signature <- df_final_survival_signature[,top_genes] %>%
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
df_final_survival_signature <- merge(df_final_survival_signature, df_clinical,
                                     by = "ID")
## Filtering operation for survival analysis
columns_survival <- c("ID", "CD44", "PFKP", "FAM50A", "CD82",
                      "LGALS1", "HIF1A", "CTSZ", "CTSH", "HSPB1",
                      "PDLIM1", "OAS2", "DBN1", "GNAI2", "ECE1",
                      "CTSD", "MORF4L2",  "CD68", "PHF11", "IFITM3",
                      "NONO", "RSU1", "CDC123", "MCU", "INPPL1",
                      "RAB8B", "TAF10", "SMAD3", "SART1", "FAM3C",
                      "PSMB9", "Condition_Mean", "OS_time", "OS_Event",
                      "PFS_time", "PFS_event", "CT1L_FLOT_1SI")
df_final_survival_signature <- df_final_survival_signature %>% 
  dplyr::select(all_of(columns_survival))
## Check
dim(df_final_survival_signature);View(df_final_survival_signature)

## Filtering
df_final_survival_signature_filtered_CT1L_FLOT <- df_final_survival_signature %>% 
  dplyr::select(all_of(colnames(df_final_survival_signature))) %>% 
  filter(df_final_survival_signature$CT1L_FLOT_1SI == 1)

## Plot
## OS
fit_os <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
               data = df_final_survival_signature)
os_plot <- ggsurvplot(fit_os, data = df_final_survival_signature,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      xlim = c(0, 48),
                      title = "Overall Survival Analysis",
                      subtitle = "Top Genes (n: 30)",
                      xlab = "Time in Months",  ylab = "OS Probability",
                      break.time.by = 4, ggtheme = theme_light(), 
                      risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                      risk.table.y.text = TRUE, conf.int.style = "step",
                      surv.median.line = "hv")
os_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_os_filtered <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
                  data = df_final_survival_signature_filtered_CT1L_FLOT)
os_plot_filtered <- ggsurvplot(fit_os_filtered, data = df_final_survival_signature_filtered_CT1L_FLOT,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      palette = c("#E7B800", "#2E9FDF"),
                      xlim = c(0, 48),
                      title = "Overall Survival Analysis",
                      subtitle = "Top Genes (n: 30), with CT1L_FLOT_1SI == 1",
                      xlab = "Time in Months",  ylab = "OS Probability",
                      break.time.by = 4, ggtheme = theme_light(),
                      risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                      risk.table.y.text = TRUE, conf.int.style = "step",
                      surv.median.line = "hv")
os_plot_filtered

## PFS
fit_pfs <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                   data = df_final_survival_signature)
pfs_plot <- ggsurvplot(fit_pfs, data = df_final_survival_signature,
                       risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                       palette = c("#E7B800", "#2E9FDF"),
                       xlim = c(0, 48),
                       title = "Progression Free Survival Analysis",
                       subtitle = "Top Genes (n: 30)",
                       xlab = "Time in Months",  ylab = "PFS Probability",
                       break.time.by = 4, ggtheme = theme_light(),
                       risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                       risk.table.y.text = TRUE, conf.int.style = "step",
                       surv.median.line = "hv")
pfs_plot

## Filtered for CT1L_FLOT_1SI == 1
fit_pfs_filtered <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                   data = df_final_survival_signature_filtered_CT1L_FLOT)
pfs_plot_filtered <- ggsurvplot(fit_pfs_filtered, data = df_final_survival_signature_filtered_CT1L_FLOT,
                       risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                       palette = c("#E7B800", "#2E9FDF"),
                       xlim = c(0, 48),
                       title = "Progression Free Survival Analysis",
                       subtitle = "Top Genes (n: 30) , with CT1L_FLOT_1SI == 1",
                       xlab = "Time in Months",  ylab = "PFS Probability",
                       break.time.by = 4, ggtheme = theme_light(),
                       risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                       risk.table.y.text = TRUE, conf.int.style = "step",
                       surv.median.line = "hv")
pfs_plot_filtered


## Top Genes 02: 16
## Filtering
df_final_survival_signature_02 <- df_final_survival %>% 
  dplyr::select(all_of(columns_to_mantain_02))
## Check
dim(df_final_survival_signature_02);View(df_final_survival_signature_02)
## Setting rownames
rownames(df_final_survival_signature_02) <- df_final_survival_signature_02$ID

## Compute metrics for signature
df_final_survival_signature_02 <- df_final_survival_signature_02[,top_genes_02] %>%
  mutate(Mean = rowMeans(across(everything()), na.rm = TRUE),  
         Median = apply(across(everything()), 1, median, na.rm = TRUE), 
         Q1 = apply(across(everything()), 1, quantile, probs = 0.25, na.rm = TRUE),  
         Q3 = apply(across(everything()), 1, quantile, probs = 0.75, na.rm = TRUE))
## Compute 75% (Q3) on mean dataset
q3_value     <- quantile(df_final_survival_signature_02$Mean, 0.75, na.rm = TRUE)
q1_value     <- quantile(df_final_survival_signature_02$Mean, 0.25, na.rm = TRUE)
median_value <- median(df_final_survival_signature_02$Mean)

## Thresholding
df_final_survival_signature_02 <- df_final_survival_signature_02 %>%
  mutate(Condition_Q3 = ifelse(Mean >= q3_value, "High", "Low"),
         Condition_Q1 = ifelse(Mean >= q1_value, "High", "Low"),
         Condition_Mean = ifelse(Mean >= median_value, "High", "Low"))
## Attach ID column
df_final_survival_signature_02$ID <- rownames(df_final_survival_signature_02)
## Merge
df_final_survival_signature_02 <- merge(df_final_survival_signature_02, df_clinical,
                                    by = "ID")
## Filtering operation for survival analysis
columns_survival_02 <- c("ID", "CD44", "PFKP", "FAM50A", "CD82",
                         "CTSZ", "CTSH", "HSPB1", "PDLIM1", "DBN1",
                         "GNAI2", "IFITM3", "RSU1", "CDC123", 
                         "RAB8B", "TAF10", "SMAD3", "Condition_Mean", "OS_time",
                         "OS_Event", "PFS_time", "PFS_event", "CT1L_FLOT_1SI")
df_final_survival_signature_02 <- df_final_survival_signature_02 %>% 
  dplyr::select(all_of(columns_survival_02))
## Check
dim(df_final_survival_signature_02);View(df_final_survival_signature_02)

## Filtering
df_final_survival_signature_02_filtered_CT1L_FLOT <- df_final_survival_signature_02 %>% 
  dplyr::select(all_of(colnames(df_final_survival_signature_02))) %>% 
  filter(df_final_survival_signature_02$CT1L_FLOT_1SI == 1)

## Plot

## OS
fit_os_02 <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
                     data = df_final_survival_signature_02)
os_plot_02 <- ggsurvplot(fit_os_02, data = df_final_survival_signature_02,
                         risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                         palette = c("#E7B800", "#2E9FDF"),
                         xlim = c(0, 48),
                         title = "Overall Survival Analysis",
                         subtitle = "Top Genes (n:16)",
                         xlab = "Time in Months",  ylab = "OS Probability",
                         break.time.by = 4, ggtheme = theme_light(),
                         risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                         risk.table.y.text = TRUE, conf.int.style = "step",
                         surv.median.line = "hv")
os_plot_02

## Filtered for CT1L_FLOT_1SI == 1
fit_os_filtered <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Condition_Mean,
                            data = df_final_survival_signature_02_filtered_CT1L_FLOT)
os_plot_filtered <- ggsurvplot(fit_os_filtered, data = df_final_survival_signature_02_filtered_CT1L_FLOT,
                                risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                                palette = c("#E7B800", "#2E9FDF"),
                                xlim = c(0, 48),
                                title = "Overall Survival Analysis",
                                subtitle = "Top Genes (n:16) , with CT1L_FLOT_1SI == 1",
                                xlab = "Time in Months",  ylab = "OS Probability",
                                break.time.by = 4, ggtheme = theme_light(),
                                risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                                risk.table.y.text = TRUE, conf.int.style = "step",
                                surv.median.line = "hv")
os_plot_filtered

## PFS
fit_pfs_02 <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                               data = df_final_survival_signature_02)
pfs_plot_02 <- ggsurvplot(fit_pfs_02, data = df_final_survival_signature_02,
                                   risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                                   palette = c("#E7B800", "#2E9FDF"),
                                   xlim = c(0, 48),
                                   title = "Progression Free Survival Analysis",
                                   subtitle = "Top Genes (n:16)",
                                   xlab = "Time in Months",  ylab = "PFS Probability",
                                   break.time.by = 4, ggtheme = theme_light(),
                                   risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                                   risk.table.y.text = TRUE, conf.int.style = "step",
                                   surv.median.line = "hv")
pfs_plot_02

## Filtered for CT1L_FLOT_1SI == 1
fit_pfs_02_filtered <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Condition_Mean,
                      data = df_final_survival_signature_02_filtered_CT1L_FLOT)
pfs_plot_02_filtered <- ggsurvplot(fit_pfs_02_filtered, data = df_final_survival_signature_02_filtered_CT1L_FLOT,
                          risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                          palette = c("#E7B800", "#2E9FDF"),
                          xlim = c(0, 48),
                          title = "Progression Free Survival Analysis",
                          subtitle = "Top Genes (n:16), with CT1L_FLOT_1SI == 1",
                          xlab = "Time in Months",  ylab = "PFS Probability",
                          break.time.by = 4, ggtheme = theme_light(),
                          risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                          risk.table.y.text = TRUE, conf.int.style = "step",
                          surv.median.line = "hv")
pfs_plot_02_filtered
