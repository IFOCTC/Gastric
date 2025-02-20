## Load useful
source("00_Support_ssGSEA.R")
library(GEOquery)
library(limma)
library(umap)
library(maptools)
library(readxl)

## List of genes/probs
## Genes
genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14", "FYN", "TUBA1A", "TNFRSF1A", "CTSZ",
           "CLIC4", "TNIP1", "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC", "SERPINH1", "ECE1",
           "FBN1", "TMSB10", "DCTN2", "PDLIM1", "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO",
           "ITGB1", "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R", "SPON2", "ABCA1",
           "TAF10", "EMP3", "CSGALNACT2", "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2",
           "SMAD3", "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B", "CD82", "ARFGAP1",
           "NRP2", "VPS26A", "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
           "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
           "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9")
## Probs
probs <- c("201040_at", "230490_x_at", "212298_at", "212334_at", "212203_x_at", "202828_s_at", "212486_s_at", "209118_s_at", "207643_s_at", "210042_s_at",
           "201560_at", "207196_s_at", "208747_s_at", "200989_at", "227068_at", "212829_at", "201105_at", "212667_at", "207714_s_at", "201749_at",
           "235318_at", "217733_s_at", "200932_s_at", "208690_s_at", "211858_x_at", "202806_at", "201508_at", "200766_at", "217767_at", "200057_s_at",
           "211945_s_at", "227308_x_at", "225673_at", "209082_s_at", "212097_at", "203507_at", "202214_s_at", "212067_s_at", "218638_s_at", "203505_at",
           "200055_at", "203729_at", "222235_s_at", "212063_at", "211284_s_at", "201136_at", "222026_at", "221816_s_at", "202796_at", "201864_at",
           "202295_s_at", "218854_at", "208869_s_at", "239952_at", "211662_s_at", "218284_at", "201889_at", "201426_s_at", "201994_at", "217765_at",
           "201300_s_at", "209191_at", "200833_s_at", "228910_at", "202211_at", "203413_at", "201807_at", "211012_s_at", "205547_s_at", "200051_at",
           "1554149_at", "202628_s_at", "226633_at", "200743_s_at", "205480_s_at", "202472_at", "201850_at", "203789_s_at", "216620_s_at", "225522_at",
           "225320_at", "209360_s_at", "205227_at", "201082_s_at", "202607_at", "208970_s_at", "223192_at", "228415_at", "216438_s_at", "202800_at",
           "202572_s_at", "212974_at", "202957_at", "209270_at", "202041_s_at", "217738_at", "213671_s_at", "200661_at", "208697_s_at", "210731_s_at",
           "201289_at", "203236_s_at")
## Combine dataframe
genes_set <- data.frame(Genes = probs)
genes_set <- as.list(genes_set)

## Load data
df_01 <- read_excel(paste0(data_path_own,
                           "/41591_2015_BFnm3850_MOESM34_ESM.xls"))
# df_02 <- read_excel(paste0(data_path_own,
#                            "/41591_2015_BFnm3850_MOESM37_ESM.xlsx"))
# df_03 <- read_excel(paste0(data_path_own,
#                            "/41591_2015_BFnm3850_MOESM38_ESM.xlsx"))
# ## Check
# dim(df_01);View(df_01) ## Clinical
# dim(df_02);View(df_02)
# dim(df_03);View(df_03)

## Load microarray
## Set data path
# dest_dir = "C:/Users/david/Documents/IFO/Final_Pipeline_Code/"
# gset_GSE62254 <- getGEO("GSE62254", destdir = dest_dir)
# gset_GSE62254 <- getGEO(filename = paste0(dest_dir, "GSE62254_series_matrix.txt.gz"))

## Extract expression
# df_expression_GSE62254 <- exprs(gset_GSE62254)
# df_expression_GSE62254 <- data.frame(df_expression_GSE62254)
## Save matrix expression
# write.csv(df_expression_GSE62254, paste(data_path_own, "df_expression_GSE62254.csv"),
#           row.names = FALSE)
## Load Data
df_expression_GSE62254 <- read.csv(paste(data_path_own, "df_expression_GSE62254.csv"),
                                   header = TRUE)
rownames(df_expression_GSE62254) <- df_expression_GSE62254$X
df_expression_GSE62254$X <- NULL
df_expression_GSE62254 <- data.frame(df_expression_GSE62254)

## Extract metadata
# df_metadata_GSE62254   <- pData(gset_GSE62254)
# write.csv(df_metadata_GSE62254, paste(data_path_own, "df_metadata_GSE62254.csv"),
#           row.names = FALSE)
df_metadata_GSE62254 <- read.csv(paste(data_path_own, "df_metadata_GSE62254.csv"),
                                   header = TRUE)
df_metadata_GSE62254 <- data.frame(df_metadata_GSE62254)
## Select only useful columns
df_metadata_GSE62254 <- df_metadata_GSE62254 %>% 
  dplyr::select(c("title", "geo_accession"))
df_metadata_GSE62254$title <- gsub("^T", "", df_metadata_GSE62254$title)
## Join
names(df_01)[names(df_01) == "Tumor ID"] <- "title"
df_clinical <- merge(df_metadata_GSE62254, df_01, by = "title")

## ***************************************
## SSGSEA - GSE14210
## ***************************************
## Compute own tpm
system.time(assign("res_GSE62254", ssgsea(as.matrix(df_expression_GSE62254), genes_set,
                                          scale = TRUE, norm = TRUE)))
## Transpose
res_transposed_GSE62254  <- t(res_GSE62254)
## Z-Score the ssgsea output for comparative analysis
res_own_GSE62254 <- (res_GSE62254 - rowMeans(res_GSE62254))/
  (rowSds(as.matrix(res_GSE62254)))[row(res_GSE62254)]
## Transpose
res_own_transposed_GSE62254  <- t(res_own_GSE62254)
## Df
res_own_transposed_GSE62254  <- data.frame(res_own_transposed_GSE62254)

## ***************************************
## SURVIVAL ANALYSIS
## ***************************************
## Set t
max_months <- 44
t <- 0
res_own_transposed_GSE62254 <- res_own_transposed_GSE62254 %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))
## Check
table(res_own_transposed_GSE62254$Signature)
## Merge
res_own_transposed_GSE62254$geo_accession <- rownames(res_own_transposed_GSE62254)
res_own_transposed_merged_GSE62254 <- merge(res_own_transposed_GSE62254, df_clinical,
                                            by = "geo_accession")
names(res_own_transposed_merged_GSE62254)[names(res_own_transposed_merged_GSE62254) == "OS\n(months)"] <- "OS_Months"
res_own_transposed_merged_GSE62254$OS_Status <- res_own_transposed_merged_GSE62254$`FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss`
res_own_transposed_merged_GSE62254 <- res_own_transposed_merged_GSE62254 %>%
  mutate(OS_Status = case_when(OS_Status %in% c(0, 1) ~ 0,
                               OS_Status %in% c(2, 3, 4) ~ 1,
                               TRUE ~ OS_Status))

## Survival Analysis
## OS
## Fit model
fit_os_GSE62254 <- survfit(Surv(OS_Months, as.integer(OS_Status)) ~ Signature,
                           data = res_own_transposed_merged_GSE62254)
## Plot
os_plot_GSE62254 <- ggsurvplot(fit_os_GSE62254, data = res_own_transposed_merged_GSE62254,
                               risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                               palette = c("#E74C3C", "#3498DB"),   
                               title = "OS",  subtitle = "GSE62254",
                               xlab = "Time (Months)",  ylab = "Survival Probability", 
                               break.time.by = 4,  
                               ggtheme = theme_classic(base_size = 16), 
                               risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                               risk.table.y.text = TRUE, conf.int.style = "step",  
                               surv.median.line = "hv",  
                               pval.size = 6,  
                               font.title = c(18, "bold"), 
                               font.subtitle = c(14, "italic"), 
                               font.x = c(14, "plain"),  
                               font.y = c(14, "plain"),  
                               font.tickslab = c(12, "plain"),  
                               legend.title = "Group",  
                               legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                               legend = c(0.85, 0.85),  
                               risk.table.fontsize = 3.5)
os_plot_GSE62254

## OS Stage III-IV
## Fit model
res_own_transposed_merged_GSE62254_stage_III_IV <- res_own_transposed_merged_GSE62254 %>% 
  dplyr::filter(pStage %in% c("III", "IV"))
fit_os_GSE62254_Stage_III_IV <- survfit(Surv(OS_Months, as.integer(OS_Status)) ~ Signature,
                                        data = res_own_transposed_merged_GSE62254_stage_III_IV)
## Plot
os_plot_GSE62254_Stage_III_IV <- ggsurvplot(fit_os_GSE62254_Stage_III_IV, data = res_own_transposed_merged_GSE62254_stage_III_IV,
                               risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                               palette = c("#E74C3C", "#3498DB"),   
                               title = "OS",  subtitle = "GSE62254 - Stage III/IV",
                               xlab = "Time (Months)",  ylab = "Survival Probability", 
                               xlim = c(0, 48),
                               break.time.by = 4,  
                               ggtheme = theme_classic(base_size = 16), 
                               risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                               risk.table.y.text = TRUE, conf.int.style = "step",  
                               surv.median.line = "hv",  
                               pval.size = 6,  
                               font.title = c(18, "bold"), 
                               font.subtitle = c(14, "italic"), 
                               font.x = c(14, "plain"),  
                               font.y = c(14, "plain"),  
                               font.tickslab = c(12, "plain"),  
                               legend.title = "Group",  
                               legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                               legend = c(0.85, 0.85),  
                               risk.table.fontsize = 3.5)
os_plot_GSE62254_Stage_III_IV

## OS Stage IV
## Fit model
res_own_transposed_merged_GSE62254_stage_IV <- res_own_transposed_merged_GSE62254 %>% 
  dplyr::filter(pStage %in% c("IV"))
fit_os_GSE62254_Stage_IV <- survfit(Surv(OS_Months, as.integer(OS_Status)) ~ Signature,
                                        data = res_own_transposed_merged_GSE62254_stage_IV)
## Plot
os_plot_GSE62254_Stage_IV <- ggsurvplot(fit_os_GSE62254_Stage_IV, data = res_own_transposed_merged_GSE62254_stage_IV,
                                            risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                                            palette = c("#E74C3C", "#3498DB"),   
                                            title = "OS",  subtitle = "GSE62254 - Stage IV",
                                            xlab = "Time (Months)",  ylab = "Survival Probability", 
                                            xlim = c(0, 48),
                                            break.time.by = 4,  
                                            ggtheme = theme_classic(base_size = 16), 
                                            risk.table.y.text.col = TRUE, risk.table.height = 0.25,  
                                            risk.table.y.text = TRUE, conf.int.style = "step",  
                                            surv.median.line = "hv",  
                                            pval.size = 6,  
                                            font.title = c(18, "bold"), 
                                            font.subtitle = c(14, "italic"), 
                                            font.x = c(14, "plain"),  
                                            font.y = c(14, "plain"),  
                                            font.tickslab = c(12, "plain"),  
                                            legend.title = "Group",  
                                            legend.labs = c("Group 1 - High", "Group 2 - Low"), 
                                            legend = c(0.85, 0.85),  
                                            risk.table.fontsize = 3.5)
os_plot_GSE62254_Stage_IV

## Boxplot for stage
res_own_transposed_merged_GSE62254_filtered <- res_own_transposed_merged_GSE62254 %>% 
  dplyr::filter(pStage %in% c("I", "II", "III", "IV"))
custom_colors <- c("I"   = "#3498DB",
                   "II"  = "#3498DB",
                   "III" = "#E74C3C",
                   "IV"  = "#E74C3C")
ggplot(res_own_transposed_merged_GSE62254_filtered, aes(x = pStage, y = Genes, fill = pStage)) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.7, outlier.shape = NA) +
  labs(title = "Gene Expression by AJCC Pathologic Stage",
       x = "Pathologic Stage", y = "Gene Expression",
       fill = "Stage", subtitle = "GSE62254 Microarray Cohort") +
  scale_fill_manual(values = custom_colors) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank()) +
  stat_compare_means(label = "p.signif", ref.group = "IV") +
  stat_compare_means(method = "kruskal.test", label.y = 3.1) +
  geom_jitter(shape = 16, color = "black", size = 1, width = 0.15, alpha = 0.5)

## Boxplot for Stage I-II Vs- III-IV
res_own_transposed_merged_GSE62254_filtered_stage <- res_own_transposed_merged_GSE62254 %>% 
  dplyr::filter(pStage %in% c("I", "II", "III", "IV"))
res_own_transposed_merged_GSE62254_filtered_stage <- res_own_transposed_merged_GSE62254_filtered_stage %>%
  mutate(pStage = case_when(pStage %in% c("I", "II") ~ "I-II",
                            pStage %in% c("III", "IV") ~ "III-IV",
                            TRUE ~ pStage))
custom_colors <- c("I-II"  = "#3498DB",
                   "III-IV" = "#E74C3C")
ggplot(res_own_transposed_merged_GSE62254_filtered_stage, aes(x = pStage, y = Genes, fill = pStage)) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.7, outlier.shape = NA) +
  labs(title = "Gene Expression by AJCC Pathologic Stage",
       x = "Pathologic Stage", y = "Gene Expression",
       fill = "Stage", subtitle = "GSE62254 Microarray Cohort") +
  scale_fill_manual(values = custom_colors) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank()) +
  stat_compare_means(method = "kruskal.test", label.y = 3.1) +
  geom_jitter(shape = 16, color = "black", size = 1, width = 0.15, alpha = 0.5)


