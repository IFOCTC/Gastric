## Load useful
library(GEOquery)
library(limma)
library(umap)
library(maptools)

## Download
## GSE14210 (N = 145)
## GSE15459 (N = 200)
## GSE22377 (N = 43)
## GSE29272 (N = 268)
## GSE51105 (N = 94)
## GSE62254 (N = 300)

## Load data
gset <- getGEO("GSE14210")[[1]]

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

## Gene Expression
df_gene_expression <- exprs(gset)
View(df_gene_expression)

## Gene Set Annotation
df_gene_set <- fData(gset)
View(df_gene_set)

## Clinical Information
df_clinical <- pData(gset)
View(df_clinical)














