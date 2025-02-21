## Load script and useful
source("00_Support_ssGSEA.R")
library(maftools)
library(mclust)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)

## Data path
data_path      <- "/Users/ifo/Desktop/Docs/IFO/Gastric/Data"
data_path_tcga <- "/Users/ifo/Desktop/Docs/IFO/Gastric/Data/TCGA"
data_path_mutations <- "/Users/ifo/Desktop/Docs/IFO/Gastric/Data/msk_met_2021"

## ***************************************
## LOAD DATA
## ***************************************
## MAF File
df_maf_stad       <- read.csv(paste0(data_path_tcga, "/df_maf_tcga.csv"))
## Signature
df_signature_stad <- read.csv(paste0(data_path_tcga, "/TCGA_Signature.csv"))
## Clinical 
df_clinical_stad  <- read_tsv(paste0(data_path_tcga, "/TCGA-STAD.clinical.tsv"))
## Mutations
df_mutated_genes  <- read.delim(paste0(data_path_tcga, "/Mutated_genes.txt"))
## Survival
df_survival_stad  <- read_tsv(paste0(data_path_tcga, "/TCGA-STAD.survival.tsv"))

## Wrangling
df_maf_stad$X <- NULL
names(df_signature_stad)[names(df_signature_stad) == "X"] <- "Tumor_Sample_Barcode"
df_signature_stad$Genes <- NULL

## Subset genes
genes <- df_mutated_genes %>% 
  dplyr::select(Gene, Is.Cancer.Gene..source..OncoKB.) %>% 
  dplyr::filter(Is.Cancer.Gene..source..OncoKB. == "Yes")

## Filter maf
df_maf_stad_filtered <- df_maf_stad %>% 
  dplyr::filter(Hugo_Symbol %in% genes$Gene)

## Subset
df_signature_stad_high <- df_signature_stad %>% 
  dplyr::filter(Signature == "High")
df_signature_stad_low  <- df_signature_stad %>% 
  dplyr::filter(Signature == "Low")
## Subset
df_maf_stad_high <- df_maf_stad_filtered
df_maf_stad_high$Tumor_Sample_Barcode <- gsub("^([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1",
                                              df_maf_stad_high$Tumor_Sample_Barcode)
df_maf_stad_low  <- df_maf_stad_filtered
df_maf_stad_low$Tumor_Sample_Barcode  <- gsub("^([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1",
                                              df_maf_stad_low$Tumor_Sample_Barcode)
df_maf_stad_filtered$Tumor_Sample_Barcode <- gsub("^([A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+).*", "\\1",
                                                  df_maf_stad_filtered$Tumor_Sample_Barcode)
## Filtering
df_maf_stad_high <- df_maf_stad_high %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% df_signature_stad_high$Tumor_Sample_Barcode)
df_maf_stad_low  <- df_maf_stad_low %>% 
  dplyr::filter(Tumor_Sample_Barcode %in% df_signature_stad_low$Tumor_Sample_Barcode)

## Merge
df_final      <- maftools::read.maf(maf = df_maf_stad_filtered, clinicalData = df_signature_stad)
df_final_high <- maftools::read.maf(maf = df_maf_stad_high, clinicalData = df_signature_stad_high)
df_final_low  <- maftools::read.maf(maf = df_maf_stad_low, clinicalData = df_signature_stad_low)

## Check
df_final

## Sample summary
getSampleSummary(df_final)
## Gene summary
getGeneSummary(df_final)
## Clinical data associated with samples
getClinicalData(df_final)

## Visualization
## Summary
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Summary_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
plotmafSummary(maf = df_final, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
dev.off()

## Barplot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Barplot_Mutated_Genes_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
mafbarplot(maf = df_final)
dev.off()

## Oncoplot
## Top ten mutated genes
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Oncoplot_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
oncoplot(maf = df_final, top = 20, draw_titv = TRUE)
dev.off()
## add signature

## Transition and Transversions
df_final_titv <- titv(maf = df_final, plot = FALSE, useSyn = TRUE)
## Plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Titv_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
plotTiTv(res = df_final_titv)
dev.off()

## Lollipop amino acid changes
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_TP53_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
lollipopPlot(maf = df_final, gene = "TP53")
dev.off()

## Rainfall plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Rainfall_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
rainfallPlot(maf = df_final, detectChangePoints = TRUE, pointSize = 0.4)
dev.off()

## Mutation load against TCGA cohorts
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Mutation_Load_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
tcgaCompare(maf = df_final, cohortName = "Example-STAD", logscale = TRUE, capture_size = 50)
dev.off()

## Somatic Interactions
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Somatic_Interactions_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
somaticInteractions(maf = df_final, top = 30, pvalue = c(0.05, 0.1))
dev.off()

## Comparing two cohorts
high_low <- mafCompare(m1 = df_final_high, m2 = df_final_low,
                       m1Name = "High", m2Name = "Low")
## Forest plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Forest_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
forestPlot(mafCompareRes = high_low, pVal = 0.05)
dev.off()

## Co onco plots
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Co_Onco_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
coOncoplot(m1 = df_final_high, m2 = df_final_low, m1Name = "High", m2Name = "Low", removeNonMutated = TRUE)
dev.off()

## Co bar plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Bar_Co_Onco_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
coBarplot(m1 = df_final_high, m2 = df_final_low, m1Name = "High", m2Name = "Low")
dev.off()

## Oncogenic signaling pathways
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Signaling_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
pws <- pathways(maf = df_final, plotType = "treemap")
dev.off()

## Plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Signaling_Res_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
plotPathways(maf = df_final, pathlist = pws)
dev.off()

## Enrichment analysis
fab_ce <- clinicalEnrichment(maf = df_final, clinicalFeature = "Signature")
## Check
fab_ce$groupwise_comparision[p_value < 0.05]
## Plot
png(filename = "/Users/ifo/Desktop/Docs/IFO/Gastric/Output/Maf_Analysis/Is_Cancer_Genes/Maf_Enrichment_Is_Cancer.png",
    width = 4000, height = 3500, units = "px", res = 300)
plotEnrichmentResults(enrich_res = fab_ce, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
dev.off()













