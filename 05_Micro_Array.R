## Load useful
source("00_Support_ssGSEA.R")
library(GEOquery)
library(limma)
library(umap)
library(maptools)
library(affy)

## Download
## GSE14210 (N = 145)
## GSE15459 (N = 200)
## GSE22377 (N = 43)
## GSE29272 (N = 268)
## GSE51105 (N = 94)
## GSE62254 (N = 300)

## *********************
## LOAD DATA
## *********************
## Set data path
dest_dir = "C:/Users/david/Documents/IFO/Final_Pipeline_Code/"

## Load data expression
options(timeout = max(300, getOption("timeout")))
gset_GSE14210 <- getGEO("GSE14210", destdir = dest_dir)
gset_GSE15459 <- getGEO("GSE15459", destdir = dest_dir)
gset_GSE22377 <- getGEO("GSE22377", destdir = dest_dir)
gset_GSE29272 <- getGEO("GSE29272", destdir = dest_dir)
gset_GSE51105 <- getGEO("GSE62254", destdir = dest_dir)
gset_GSE62254 <- getGEO("GSE62254", destdir = dest_dir)

## Clinical
df_clinical_GSE14210 <- read.table("GSE14208_SurvivalData.txt", header = TRUE)
df_clinical_GSE15459 <- read_excel("GSE15459_outcome.xls")
df_clinical_GSE22377
df_clinical_GSE29272
df_clinical_GSE51105

## Format ID
df_clinical_GSE14210$ID <- gsub("\\.CEL$", "", df_clinical_GSE14210$ID)

## *********************
## WRANGLING
## *********************

## Extract
gset_GSE14210 <- gset_GSE14210[[1]]
gset_GSE15459 <- gset_GSE15459[[1]]
gset_GSE22377 <- gset_GSE22377[[1]]
gset_GSE29272 <- gset_GSE29272[[1]]
gset_GSE51105 <- gset_GSE51105[[1]]

## Expression
df_expression_GSE14210 <- exprs(gset_GSE14210)
df_expression_GSE15459 <- exprs(gset_GSE15459)
df_expression_GSE22377 <- exprs(gset_GSE22377)
df_expression_GSE29272 <- exprs(gset_GSE29272)
df_expression_GSE51105 <- exprs(gset_GSE51105)

## Check
dim(df_expression_GSE14210)
dim(df_expression_GSE15459)
dim(df_expression_GSE22377)
dim(df_expression_GSE29272)
dim(df_expression_GSE51105)

## Format
df_expression_GSE14210 <- data.frame(df_expression_GSE14210)
df_expression_GSE15459 <- data.frame(df_expression_GSE15459)
df_expression_GSE22377 <- data.frame(df_expression_GSE22377)
df_expression_GSE29272 <- data.frame(df_expression_GSE29272)
df_expression_GSE51105 <- data.frame(df_expression_GSE51105)

## Add ID Column
df_expression_GSE14210$ID <- rownames(df_expression_GSE14210)
df_expression_GSE15459$ID <- rownames(df_expression_GSE15459)
df_expression_GSE22377$ID <- rownames(df_expression_GSE22377)
df_expression_GSE29272$ID <- rownames(df_expression_GSE29272)
df_expression_GSE51105$ID <- rownames(df_expression_GSE51105)

## Add label dataset
df_expression_GSE14210$Label_GSE14210 <- rep("GSE14210", rownames(df_expression_GSE14210))
df_expression_GSE15459$Label_GSE15459 <- rep("GSE15459", rownames(df_expression_GSE15459))
df_expression_GSE22377$Label_GSE22377 <- rep("GSE22377", rownames(df_expression_GSE22377))
df_expression_GSE29272$Label_GSE29272 <- rep("GSE29272", rownames(df_expression_GSE29272))
df_expression_GSE51105$Label_GSE51105 <- rep("GSE51105", rownames(df_expression_GSE51105))

## Merge
dfs <- list(df_expression_GSE14210, df_expression_GSE15459, df_expression_GSE22377,
            df_expression_GSE29272, df_expression_GSE51105)
df_final_merged <- Reduce(function(x, y) merge(x, y, by = "ID"), dfs)

## Check
dim(df_final_merged)

## Genes
genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3",
           "MMP14", "FYN", "TUBA1A", "TNFRSF1A", "CTSZ",
           "CLIC4", "TNIP1", "C1S", "HIF1A", "PGK1",
           "PIP4K2A", "LGALS1", "SPARC", "SERPINH1", "ECE1",
           "FBN1", "TMSB10", "DCTN2", "PDLIM1", "GNAS",
           "DBN1", "IGFBP4", "CTSD", "C3", "NONO",
           "ITGB1", "LTBP3", "MYADM", "COL18A1", "CAV1",
           "CD68", "CUL4B", "C1R", "SPON2", "ABCA1",
           "TAF10", "EMP3", "CSGALNACT2", "CD44", "GRN",
           "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
           "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2",
           "SMAD3", "FAM3C", "VIM", "MORF4L2", "NRBP1",
           "PRNP", "TUBB6", "RAP1B", "CD82", "ARFGAP1",
           "NRP2", "VPS26A", "PML", "TAGLN", "SART1",
           "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
           "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
           "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
           "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
           "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP",
           "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8",
           "CYR61", "LGALS9")
## Probs
probs <- c("201040_at", "230490_x_at", "212298_at", "212334_at", "212203_x_at",
           "202828_s_at", "212486_s_at", "209118_s_at", "207643_s_at", "210042_s_at",
           "201560_at", "207196_s_at", "208747_s_at", "200989_at", "227068_at",
           "212829_at", "201105_at", "212667_at", "207714_s_at", "201749_at",
           "235318_at", "217733_s_at", "200932_s_at", "208690_s_at", "211858_x_at",
           "202806_at", "201508_at", "200766_at", "217767_at", "200057_s_at",
           "211945_s_at", "227308_x_at", "225673_at", "209082_s_at", "212097_at",
           "203507_at", "202214_s_at", "212067_s_at", "218638_s_at", "203505_at",
           "200055_at", "203729_at", "222235_s_at", "212063_at", "211284_s_at",
           "201136_at", "222026_at", "221816_s_at", "202796_at", "201864_at",
           "202295_s_at", "218854_at", "208869_s_at", "239952_at", "211662_s_at",
           "218284_at", "201889_at", "201426_s_at", "201994_at", "217765_at",
           "201300_s_at", "209191_at", "200833_s_at", "228910_at", "202211_at",
           "203413_at", "201807_at", "211012_s_at", "205547_s_at", "200051_at",
           "1554149_at", "202628_s_at", "226633_at", "200743_s_at", "205480_s_at",
           "202472_at", "201850_at", "203789_s_at", "216620_s_at", "225522_at",
           "225320_at", "209360_s_at", "205227_at", "201082_s_at", "202607_at",
           "208970_s_at", "223192_at", "228415_at", "216438_s_at", "202800_at",
           "202572_s_at", "212974_at", "202957_at", "209270_at", "202041_s_at",
           "217738_at", "213671_s_at", "200661_at", "208697_s_at", "210731_s_at",
           "201289_at", "203236_s_at")

## Combine dataframe
genes_set <- data.frame(Genes = probs)
genes_set <- as.list(genes_set)

## Gene Expression
df_gene_expression_GSE14210 <- exprs(gset_GSE14210)
df_gene_expression_GSE14210 <- data.frame(df_gene_expression_GSE14210)

## Clinical Information
sample_info_GSE14210 <- pData(gset_GSE14210)
names(sample_info_GSE14210)[names(sample_info_GSE14210) == "title"] <- "ID"
## Merge
df_clinical_merged_GSE14210 <- merge(df_clinical_GSE14210, sample_info_GSE14210, by = "ID")

## ***************************************
## SSGSEA
## ***************************************
## Compute own tpm
system.time(assign("res_GSE14210", ssgsea(as.matrix(df_gene_expression_GSE14210), genes_set,
                                 scale = TRUE, norm = TRUE)))
## Transpose
res_transposed_GSE14210  <- t(res_GSE14210)
## Z-Score the ssgsea output for comparative analysis
res_own_GSE14210 <- (res_GSE14210 - rowMeans(res_GSE14210))/
  (rowSds(as.matrix(res_GSE14210)))[row(res_GSE14210)]
## Transpose
res_own_transposed_GSE14210  <- t(res_own_GSE14210)
## Df
res_own_transposed_GSE14210  <- data.frame(res_own_transposed_GSE14210)

## ***************************************
## SURVIVAL ANALYSIS
## ***************************************
## Set t
max_months <- 44
t <- 0
res_own_transposed_GSE14210 <- res_own_transposed_GSE14210 %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))
res_own_transposed_GSE14210$geo_accession <- rownames(res_own_transposed_GSE14210)
## Merge
res_own_transposed_merged_GSE14210 <- merge(res_own_transposed_GSE14210, df_clinical_merged_GSE14210, by = "geo_accession")
## Check
table(res_own_transposed_GSE14210$Signature)

## OS
p_os_status_distribution_GSE14210 <- ggplot(res_own_transposed_merged_GSE14210, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),
           fill = c("#E74C3C", "#3498DB"), color = c("#E74C3C", "#3498DB"), alpha = 0.8) + 
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), 
            stat = "count", vjust = -0.5, size = 5, fontface = "bold") + 
  scale_y_continuous(labels = percent_format()) + 
  labs(title = "OS Status Distribution - Micro-Array",
       subtitle = "GSE14210",
       x = "Signature",
       y = "Percentage") +
  theme_classic(base_size = 14, base_family = "Arial") + 
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),  
        axis.title.x = element_text(vjust = -0.5),  
        axis.title.y = element_text(vjust = 1.5),   
        axis.text.x = element_text(face = "bold"),  
        axis.text.y = element_text(face = "bold"), 
        panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),  
        legend.position = "none")
p_os_status_distribution_GSE14210
## Fit model
fit_os_GSE14210 <- survfit(Surv(overallSurvival_mo, as.integer(death1_censor0)) ~ Signature,
                           data = res_own_transposed_merged_GSE14210)
## Plot
os_plot_GSE14210 <- ggsurvplot(fit_os_GSE14210, data = res_own_transposed_merged_GSE14210,
                          risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                          palette = c("#E74C3C", "#3498DB"),  
                          xlim = c(0, max_months),  
                          title = "OS - IRE",  subtitle = "GSE14210",
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
os_plot_GSE14210

## PFS
## Fit model
fit_pfs_GSE14210 <- survfit(Surv(timeToProgression_mo, as.integer(progression1_censor0)) ~ Signature,
                           data = res_own_transposed_merged_GSE14210)
## Plot
pfs_plot_GSE14210 <- ggsurvplot(fit_pfs_GSE14210, data = res_own_transposed_merged_GSE14210,
                               risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                               palette = c("#E74C3C", "#3498DB"),  
                               xlim = c(0, max_months),  
                               title = "PFS - IRE",  subtitle = "GSE14210",
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
pfs_plot_GSE14210
















