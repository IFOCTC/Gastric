## Code useful to perform GSVA

## ***********************************************************
## DATA LOADING
## ***********************************************************

## Load useful
source("00_Support.R")

## Load tcga tpm
df_tpm_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.star_tpm.tsv"),
                        show_col_types = FALSE)
## Load own tpm
df_tpm_own  <- read_tsv(paste0(data_path_own, "/salmon.merged.gene_tpm.tsv"),
                        show_col_types = FALSE)

## ***********************************************************
## WRANGLING
## ***********************************************************

## Log transformation
exclude_cols <- c("gene_id", "gene_name")
df_tpm_own_filtered <- df_tpm_own %>% 
  dplyr::select(-all_of(exclude_cols)) %>%
  mutate(across(everything(), ~ log2(. + 1)))
df_tpm_own <- df_tpm_own %>%
  dplyr::select(all_of(exclude_cols)) %>%
  bind_cols(df_tpm_own_filtered)

## Load survival data tcga
df_survival_tcga <- read_tsv(paste0(output_path, "/survival_cleaned.tsv"), show_col_types = FALSE)
df_clinical_tcga <- read_tsv(paste0(data_path, "/TCGA-STAD.clinical.tsv"), show_col_types = FALSE)
df_clinical_tcga_filtered <- df_clinical_tcga %>%
  filter(ajcc_pathologic_stage.diagnoses %in% c("Stage III", "Stage IIIA",
                                                "Stage IIIB", "Stage IIIC", "Stage IV"))
ids_df_clinical_tcga <- df_clinical_tcga_filtered$sample
df_survival_tcga_filtered     <- df_survival_tcga %>%
  filter(sample %in% ids_df_clinical_tcga)

## Load survival data own
df_survival_own <- read_excel(paste0(data_path_own, "/stomaci_wes_all.xlsx"))
df_survival_own_filtered_os  <- df_survival_own %>%
  filter(CT1L_FLOT_1SI == 1)
df_survival_own_filtered_pfs <- df_survival_own %>%
  filter(CT1L_1SI == 1)

## Ensembl id gene name
Ensembl_ID                     <- df_tpm_tcga$Ensembl_ID
Ensembl_ID_cleaned             <- sub("\\..*", "", Ensembl_ID) 
df_tpm_tcga_cleaned            <- df_tpm_tcga
df_tpm_tcga_cleaned$Ensembl_ID <- Ensembl_ID_cleaned
## Attach the gene name of own genes to tpm tcga
df_support_gene_name <- df_tpm_own %>% 
  dplyr::select(gene_id, gene_name)
## Change name column for merge
df_support_gene_name <- df_support_gene_name %>% 
  rename(Ensembl_ID = gene_id)
## Merge
df_tpm_tcga_final <- df_tpm_tcga_cleaned %>% 
  left_join(df_support_gene_name, by = "Ensembl_ID")
## Move the column
df_tpm_tcga_final <- df_tpm_tcga_final %>% 
  dplyr::select(1, gene_name, everything())
## Remove gene id
df_tpm_tcga_final$Ensembl_ID <- NULL
df_tpm_own$gene_id <- NULL
## Export row names
genes_names_tcga <- df_tpm_tcga_final$gene_name
genes_names_tpm  <- df_tpm_own$gene_name
## Remove column
df_tpm_tcga_final$gene_name <- NULL
df_tpm_own$gene_name        <- NULL

## TCGA TPM
# row_num_t   <- 5
# min_samples <- 10
# df_tpm_tcga_filtered <- df_tpm_tcga_final[rowSums(df_tpm_tcga_final >= row_num_t) >= min_samples, ]
# df_tpm_tcga_final <- df_tpm_tcga_final[rowSums(df_tpm_tcga_final >= row_num_t), ]
## Consider only tumor samples
selected_colnames_tcga <- colnames(df_tpm_tcga_final)[grep("01A$", colnames(df_tpm_tcga_final))]
## Filtering
df_tpm_tcga_cleaned_tumor <- df_tpm_tcga_final %>% 
  dplyr::select(any_of(selected_colnames_tcga))
## Remove the last part of barcode
selected_colnames_cleaned_tcga <- sub("-[^-]+$", "", colnames(df_tpm_tcga_cleaned_tumor))
colnames(df_tpm_tcga_cleaned_tumor) <- selected_colnames_cleaned_tcga

## OWN TPM
# row_num_t   <- 5
# min_samples <- 10
# df_tpm_own_filtered <- df_tpm_own[rowSums(df_tpm_own >= row_num_t) >= min_samples, ]
# df_tpm_own <- df_tpm_own[rowSums(df_tpm_own >= row_num_t), ]
## Consider only tumor patients
df_tpm_own_filtered_tumor <- df_tpm_own %>% 
  dplyr::select(matches("T$"))

## Setting row names
rownames(df_tpm_tcga_cleaned_tumor) <- make.names(genes_names_tcga, unique = TRUE) 
rownames(df_tpm_own_filtered_tumor) <- make.names(genes_names_tpm, unique = TRUE) 

## Define the list of genes to use

## q = 10%, 178 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R",
#            "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2", "CD44",
#            "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH",
#            "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3", "FAM3C",
#            "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B",
#            "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML", "TAGLN",
#            "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
#            "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1", "MCU", "RUNX1",
#            "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2",
#            "TMSB4X", "SLC1A3", "DLGAP4", "DENND3", "HCLS1", "LAMB3",
#            "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61",
#            "LGALS9", "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4",
#            "CTGLF12P", "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48",
#            "MAT2A", "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1", "TNFAIP3",
#            "LOXL2", "SLC39A13", "TNFAIP2", "DPYD", "ZDHHC6", "FAM50A",
#            "R3HDM2", "PFKFB3", "SVIL", "PFKP", "WFS1", "INPPL1", "OAS2",
#            "EHBP1", "IL1RN", "SLC25A37", "OR7E38P", "RBM34", "LAMA3",
#            "F3", "ELL2", "NDUFV1", "AC138783.12", "SULF2", "DCBLD2",
#            "PDE8A", "PITRM1", "ST3GAL3", "TOM1L2", "RRAS2", "WARS",
#            "ST5", "IL8", "NET1", "F8A1", "AGAP9", "FER1L4", "TAZ",
#            "SGK1", "SPATS2", "RARRES3", "TMEM185A", "HIBADH", "INO80B",
#            "IFI35", "CNOT2", "RALGPS2", "WDR45", "NAMPTL", "TMSB4XP6",
#            "FAM107B", "SHMT2")
## only > 
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R",
#            "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2", "CD44",
#            "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH",
#            "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3", "FAM3C",
#            "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B",
#            "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML", "TAGLN",
#            "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
#            "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1", "MCU", "RUNX1",
#            "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2",
#            "TMSB4X", "SLC1A3", "DLGAP4", "DENND3", "HCLS1", "LAMB3",
#            "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61",
#            "LGALS9", "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4",
#            "CTGLF12P", "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48",
#            "MAT2A", "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1", "TNFAIP3",
#            "LOXL2", "SLC39A13", "TNFAIP2", "DPYD", "ZDHHC6", "FAM50A",
#            "R3HDM2", "PFKFB3", "SVIL", "PFKP", "WFS1", "INPPL1", "OAS2",
#            "EHBP1", "IL1RN", "SLC25A37", "OR7E38P", "RBM34", "LAMA3",
#            "F3", "ELL2", "NDUFV1", "AC138783.12", "SULF2", "DCBLD2",
#            "PDE8A", "PITRM1", "ST3GAL3", "TOM1L2", "RRAS2", "WARS",
#            "ST5", "IL8", "NET1", "F8A1", "AGAP9", "FER1L4", "TAZ",
#            "SGK1", "SPATS2", "RARRES3", "TMEM185A", "HIBADH", "INO80B",
#            "IFI35", "CNOT2", "RALGPS2", "WDR45", "NAMPTL", "TMSB4XP6")

## q = 15%, 168 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R",
#            "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2", "CD44",
#            "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH",
#            "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3", "FAM3C",
#            "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B",
#            "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML", "TAGLN",
#            "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
#            "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1", "MCU", "RUNX1",
#            "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2",
#            "TMSB4X", "SLC1A3", "DLGAP4", "DENND3", "HCLS1", "LAMB3",
#            "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61",
#            "LGALS9", "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4",
#            "CTGLF12P", "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48",
#            "MAT2A", "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1", "TNFAIP3",
#            "LOXL2", "SLC39A13", "TNFAIP2", "DPYD", "ZDHHC6", "FAM50A",
#            "R3HDM2", "PFKFB3", "SVIL", "PFKP", "WFS1", "INPPL1", "OAS2",
#            "EHBP1", "IL1RN", "SLC25A37", "OR7E38P", "RBM34", "LAMA3",
#            "F3", "ELL2", "NDUFV1", "AC138783.12", "SULF2", "DCBLD2",
#            "PDE8A", "PITRM1", "ST3GAL3", "TOM1L2", "RRAS2", "WARS",
#            "ST5", "IL8", "NET1", "F8A1", "AGAP9", "FER1L4", "TAZ",
#            "SGK1", "SPATS2", "RARRES3", "TMEM185A")

## q = 20%, 162 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R",
#            "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2", "CD44",
#            "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH",
#            "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3", "FAM3C",
#            "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B",
#            "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML", "TAGLN",
#            "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
#            "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1", "MCU", "RUNX1",
#            "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2",
#            "TMSB4X", "SLC1A3", "DLGAP4", "DENND3", "HCLS1", "LAMB3",
#            "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61",
#            "LGALS9", "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4",
#            "CTGLF12P", "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48",
#            "MAT2A", "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1", "TNFAIP3",
#            "LOXL2", "SLC39A13", "TNFAIP2", "DPYD", "ZDHHC6", "FAM50A",
#            "R3HDM2", "PFKFB3", "SVIL", "PFKP", "WFS1", "INPPL1", "OAS2",
#            "EHBP1", "IL1RN", "SLC25A37", "OR7E38P", "RBM34", "LAMA3",
#            "F3", "ELL2", "NDUFV1", "AC138783.12", "SULF2", "DCBLD2",
#            "PDE8A", "PITRM1", "ST3GAL3", "TOM1L2", "RRAS2", "WARS",
#            "ST5", "IL8", "NET1", "F8A1", "AGAP9")
## only > 
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B", "C1R",
#            "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2", "CD44",
#            "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1", "CTSH",
#            "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3", "FAM3C",
#            "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6", "RAP1B",
#            "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML", "TAGLN",
#            "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1", "UGP2",
#            "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1", "MCU", "RUNX1",
#            "IL1RAP", "DCTN1", "NDST1", "UROD", "SLC25A28", "AP1S2",
#            "TMSB4X", "SLC1A3", "DLGAP4", "DENND3", "HCLS1", "LAMB3",
#            "FIBP", "NAMPT", "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61",
#            "LGALS9", "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4",
#            "CTGLF12P", "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48",
#            "MAT2A", "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1", "TNFAIP3",
#            "LOXL2", "SLC39A13", "TNFAIP2", "DPYD", "ZDHHC6", "FAM50A",
#            "R3HDM2", "PFKFB3", "SVIL", "PFKP", "WFS1", "INPPL1", "OAS2",
#            "EHBP1", "IL1RN", "SLC25A37", "OR7E38P", "RBM34", "LAMA3",
#            "F3", "ELL2", "NDUFV1", "AC138783.12", "SULF2", "DCBLD2",
#            "PDE8A", "PITRM1", "ST3GAL3", "TOM1L2")

## q = 25%, 150 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6",
#            "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML",
#            "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1",
#            "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
#            "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD",
#            "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3", "DLGAP4",
#            "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS",
#            "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9", "SEMA4A",
#            "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48", "MAT2A",
#            "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1",
#            "TNFAIP3", "LOXL2", "SLC39A13", "TNFAIP2", "DPYD",
#            "ZDHHC6", "FAM50A", "R3HDM2", "PFKFB3", "SVIL",
#            "PFKP", "WFS1", "INPPL1", "OAS2", "EHBP1", "IL1RN",
#            "SLC25A37", "OR7E38P", "RBM34", "LAMA3", "F3", "ELL2",
#            "NDUFV1", "AC138783.12", "SULF2")

## q = 30%, 138 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6",
#            "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML",
#            "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1",
#            "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
#            "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD",
#            "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3", "DLGAP4",
#            "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS",
#            "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9", "SEMA4A",
#            "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48", "MAT2A",
#            "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1",
#            "TNFAIP3", "LOXL2", "SLC39A13", "TNFAIP2", "DPYD",
#            "ZDHHC6", "FAM50A", "R3HDM2", "PFKFB3", "SVIL",
#            "PFKP", "WFS1", "INPPL1")
## only > 


## q = 35%, 130 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6",
#            "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML",
#            "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1",
#            "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
#            "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD",
#            "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3", "DLGAP4",
#            "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS",
#            "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9", "SEMA4A",
#            "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48", "MAT2A",
#            "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9", "ANKIB1", "C16orf62", "HSPB1",
#            "TNFAIP3", "LOXL2", "SLC39A13", "TNFAIP2", "DPYD")

# ## q = 40%, 122 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6",
#            "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML",
#            "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1",
#            "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
#            "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD",
#            "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3", "DLGAP4",
#            "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS",
#            "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9", "SEMA4A",
#            "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48", "MAT2A",
#            "FMNL1", "CDK17", "CHFR", "IRF1", "TRAFD1", "CDC123",
#            "IFITM1", "PSMB9")
## only > 115 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP", "TUBB6",
#            "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A", "PML",
#            "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B", "TPP1",
#            "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10", "AAK1",
#            "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1", "UROD",
#            "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3", "DLGAP4",
#            "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT", "MARS",
#            "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9", "SEMA4A",
#            "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2", "ATP2B4", "OPTN", "BAG3", "C15orf48", "MAT2A",
#            "FMNL1")

## q = 45%, 109 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
#            "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
#            "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT",
#            "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9",
#            "SEMA4A", "MAP7D1", "TNFSF10", "PLCG1", "HERC4", "CTGLF12P",
#            "TPM2")
## only > 115 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
#            "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
#            "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT",
#            "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9",
#            "SEMA4A")

## q = 50%, 103 genes
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
#            "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
#            "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT",
#            "MARS", "CTSA", "EIF3E", "LGALS8", "CYR61", "LGALS9",
#            "SEMA4A")

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


## only > 
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
#            "UROD", "SLC25A28", "AP1S2", "TMSB4X", "SLC1A3",
#            "DLGAP4", "DENND3", "HCLS1", "LAMB3", "FIBP", "NAMPT",
#            "MARS")

## q = 55%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1", "NDST1",
#            "UROD", "SLC25A28", "AP1S2", "TMSB4X")
## only > 

## q = 60%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI", "CAPG", "SEMA3C", "ARHGEF10",
#            "AAK1", "MCU", "RUNX1", "IL1RAP", "DCTN1")
## only >
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2", "NRBP1", "PRNP",
#            "TUBB6", "RAP1B", "CD82", "ARFGAP1", "NRP2", "VPS26A",
#            "PML", "TAGLN", "SART1", "CLDND1", "SERPINE1", "RAB8B",
#            "TPP1", "UGP2", "MPI")

## q = 70%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3","MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10", "EMP3", "CSGALNACT2",
#            "CD44", "GRN", "PLP2", "RBM3", "PHF11", "SYNPO", "GDI1",
#            "CTSH", "DSE", "GABARAPL1", "ZEB1", "VDAC2", "SMAD3",
#            "FAM3C", "VIM", "MORF4L2")

## q = 75%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14", "FYN",
#            "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1", "C1S",
#            "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC", "SERPINH1",
#            "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1", "GNAS", "DBN1",
#            "IGFBP4", "CTSD", "C3", "NONO", "ITGB1", "LTBP3", "MYADM",
#            "COL18A1", "CAV1", "CD68", "CUL4B", "C1R", "SPON2", "ABCA1",
#            "TAF10", "EMP3", "CSGALNACT2", "CD442", "GRN", "PLP2", "RBM3",
#            "PHF11", "SYNPO", "GDI1")

## q = 80%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1", "CAV1", "CD68", "CUL4B",
#            "C1R", "SPON2", "ABCA1", "TAF10")
## only >
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1", "ECE1", "FBN1", "TMSB10", "DCTN2", "PDLIM1",
#            "GNAS", "DBN1", "IGFBP4", "CTSD", "C3", "NONO", "ITGB1",
#            "LTBP3", "MYADM", "COL18A1")

## q = 90%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3", "MMP14",
#            "FYN", "TUBA1A", "TNFRSF1A", "CTSZ", "CLIC4", "TNIP1",
#            "C1S", "HIF1A", "PGK1", "PIP4K2A", "LGALS1", "SPARC",
#            "SERPINH1")

## q = 95%
# genes <- c("GNAI2", "RSU1", "NRP1", "GNS", "IFITM3",
#            "MMP14", "FYN", "TUBA1A", "TNFRSF1A", "CTSZ")



genes_set <- data.frame(Genes = genes)
genes_set <- as.list(genes_set)


## *********************************************************** 
## SSGSEA
## ***********************************************************

## Compute own tpm
system.time(assign("res_own_tpm", ssgsea(as.matrix(df_tpm_own_filtered_tumor), genes_set,
                                 scale = TRUE, norm = TRUE)))
## Compute tcga tpm
system.time(assign("res_tcga_tpm", ssgsea(as.matrix(df_tpm_tcga_cleaned_tumor), genes_set,
                                 scale = TRUE, norm = TRUE)))
## Transpose
res_own_tpm_transposed  <- t(res_own_tpm)
res_tcga_tpm_transposed <- t(res_tcga_tpm)

## Z-Score the ssgsea output for comparative analysis
mat_own_tpm  <- (res_own_tpm - rowMeans(res_own_tpm))/
                (rowSds(as.matrix(res_own_tpm)))[row(res_own_tpm)]
mat_tcga_tpm <- (res_tcga_tpm - rowMeans(res_tcga_tpm))/
                (rowSds(as.matrix(res_tcga_tpm)))[row(res_tcga_tpm)]
## Transpose
mat_own_tpm_transposed  <- t(mat_own_tpm)
mat_tcga_tpm_transposed <- t(mat_tcga_tpm)

## df
mat_own_tpm_transposed  <- data.frame(mat_own_tpm_transposed)
mat_tcga_tpm_transposed <- data.frame(mat_tcga_tpm_transposed)

## Number of patients before signature
## OWN  -> 182
## TCGA -> 410

## Different values for t
## Survival Analysis - threshold for t 
# perform_survival_analysis(t = 0.5,
#                           df_1 = mat_own_tpm_transposed, df_2 = mat_tcga_tpm_transposed,
#                           max_months = 48)

## ***********************************************************
## SURVIVAL ANALYSIS
## ***********************************************************

## Set t
max_months <- 48
t <- 0
mat_own_tpm_transposed <- mat_own_tpm_transposed %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))
mat_tcga_tpm_transposed <- mat_tcga_tpm_transposed %>%
  mutate(Signature = case_when(Genes >= t ~ "High",
                               Genes <  t ~ "Low"))
## Check
table(mat_own_tpm_transposed$Signature)
table(mat_tcga_tpm_transposed$Signature)

## Wrangling for survival analysis
## TCGA Tpm
mat_tcga_tpm_transposed$ID <- rownames(mat_tcga_tpm_transposed)
mat_tcga_tpm_transposed_merged <- merge(mat_tcga_tpm_transposed, df_survival_tcga, by = "ID")
mat_tcga_tpm_transposed_merged_filtered <- merge(mat_tcga_tpm_transposed,
                                                 df_survival_tcga_filtered, by = "ID")
mat_tcga_tpm_transposed_merged$OS.time.month <- mat_tcga_tpm_transposed_merged$OS.time / 30
mat_tcga_tpm_transposed_merged_filtered$OS.time.month <- mat_tcga_tpm_transposed_merged_filtered$OS.time / 30
  
## Own TPM
mat_own_tpm_transposed$ID <- rownames(mat_own_tpm_transposed)
mat_own_tpm_transposed$ID <- gsub("\\.1T$", "", mat_own_tpm_transposed$ID)
mat_own_tpm_transposed$ID <- gsub("\\.", "-", mat_own_tpm_transposed$ID)
mat_own_tpm_transposed_merged_os  <- merge(mat_own_tpm_transposed, df_survival_own_filtered_os, by = "ID")
mat_own_tpm_transposed_merged_pfs <- merge(mat_own_tpm_transposed, df_survival_own_filtered_pfs, by = "ID")

## Now, export our table with signature created, to use
## in script 12_Multivariate_Analysis.R in Code Gastric
# writexl::write_xlsx(mat_own_tpm_transposed_merged_os,
#                     "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Data/df_signature_os_filtered.xlsx")
# writexl::write_xlsx(mat_own_tpm_transposed_merged_pfs,
#                     "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Data/df_signature_pfs_filtered.xlsx")
# writexl::write_xlsx(mat_own_tpm_transposed,
#                     "C:/Users/david/Documents/IFO/Gastric_TCGA/Output_Data/df_signature.xlsx")

## After the merging the number of patients is:
## OWN OS           -> 155
## OWN PFS          -> 148
## TCGA ALL SAMPLES -> 383
## TCGA FILTERED    -> 197
  
## Signature Density Distribution
## OS IRE 
p_signature_os <- ggplot(mat_own_tpm_transposed_merged_os, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature OS IRE") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## PFS IRE 
p_signature_pfs <- ggplot(mat_own_tpm_transposed_merged_pfs, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature PFS IRE") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## TCGA
p_signature_tcga <- ggplot(mat_tcga_tpm_transposed_merged, aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature TCGA") +
  xlab("") +
  ylab("Count") +
  theme_minimal()
## TCGA Filtered
p_signature_filtered_tcga <- ggplot(mat_tcga_tpm_transposed_merged_filtered,
                                    aes(x = Genes, y = after_stat(count))) +
  geom_histogram(fill = "orange", color = "white", bins = 30) +
  ggtitle("Signature Filetered TCGA") +
  xlab("") +
  ylab("Count") +
  theme_minimal()

## SIgnature Barplot Distribution
## OS Status distribution
p_os_status_distribution <- ggplot(mat_own_tpm_transposed_merged_os, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  ## Calculate percentage relative to total
           fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), stat = "count", vjust = -0.5) +
  scale_y_continuous(labels = scales::percent_format()) +  ## Scale y-axis to percentage
  labs(title = "Category OS IRE", x = "", y = "Percentage") +
  theme_minimal() 
p_pfs_status_distribution <- ggplot(mat_own_tpm_transposed_merged_pfs, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  ## Calculate percentage relative to total
           fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), stat = "count", vjust = -0.5) +
  scale_y_continuous(labels = scales::percent_format()) +  ## Scale y-axis to percentage
  labs(title = "Category PFS IRE", x = "", y = "Percentage") +
  theme_minimal() 
## TCGA TPM plot with percentage on the y-axis
p_os_status_distribution_tcga <- ggplot(mat_tcga_tpm_transposed_merged, aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  ## Calculate percentage relative to total
           fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), stat = "count", vjust = -0.5) +
  scale_y_continuous(labels = scales::percent_format()) +  ## Scale y-axis to percentage
  labs(title = "Category OS TCGA", x = "", y = "Percentage") +
  theme_minimal()
## OS Status Distribution Filtered
p_os_status_distribution_filtered_tcga <- ggplot(mat_tcga_tpm_transposed_merged_filtered,
                                                 aes(x = Signature)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  ## Calculate percentage relative to total
           fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(aes(label = paste0(..count.., " (", round((..count..)/sum(..count..) * 100, 1), "%)"), 
                y = (..count..)/sum(..count..)), stat = "count", vjust = -0.5) +
  scale_y_continuous(labels = scales::percent_format()) +  # Scale y-axis to percentage
  labs(title = "Category OS Filtered TCGA", x = "", y = "Percentage") +
  theme_minimal()
  
  
## ****************************
## Survival Analysis - IRE
fit_os_own <- survfit(Surv(OS_time, as.integer(OS_Event)) ~ Signature,
                      data = mat_own_tpm_transposed_merged_os)
os_plot_own <- ggsurvplot(fit_os_own, data = mat_own_tpm_transposed_merged_os,
                           risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                           palette = c("#E74C3C", "#3498DB", "#F39C12"),
                           xlim = c(0, max_months),
                           title = "OS IRE",
                           subtitle = paste0("Top Genes:", length(genes_set$Gene),
                                             ", t:", t, " - CT1L_FLOT_1SI Samples"),
                           xlab = "Time (Months)",  ylab = "OS Probability",
                           break.time.by = 4, ggtheme = theme_light(),
                           risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                           risk.table.y.text = TRUE, conf.int.style = "step",
                           surv.median.line = "hv")

## ****************************
## Progression Free Survival - IRE
fit_pfs_own <- survfit(Surv(PFS_time, as.integer(PFS_event)) ~ Signature,
                       data = mat_own_tpm_transposed_merged_pfs)
pfs_plot_own <- ggsurvplot(fit_pfs_own, data = mat_own_tpm_transposed_merged_pfs,
                            risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                            palette = c("#E74C3C", "#3498DB", "#F39C12"),
                            xlim = c(0, max_months),
                            title = "PFS IRE",
                            subtitle = paste0("Top Genes:", length(genes_set$Gene),
                                              ", t:", t, " - CT1L_1SI Samples"),
                            xlab = "Time (Months)",  ylab = "PFS Probability",
                            break.time.by = 4, ggtheme = theme_light(),
                            risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                            risk.table.y.text = TRUE, conf.int.style = "step",
                            surv.median.line = "hv")

## ****************************
## Survival Analysis - TCGA 
fit_os_tcga <- survfit(Surv(OS.time.month, as.integer(OS)) ~ Signature,
                       data = mat_tcga_tpm_transposed_merged)
os_plot_tcga <- ggsurvplot(fit_os_tcga, data = mat_tcga_tpm_transposed_merged,
                            risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                            palette = c("#E74C3C", "#3498DB", "#F39C12"),
                            xlim = c(0, max_months),
                            title = "OS TCGA",
                            subtitle = paste0("Top Genes:", length(genes_set$Gene),
                                              ", t:", t, " - All Samples"),
                            xlab = "Time (Months)",  ylab = "OS Probability",
                            break.time.by = 4, ggtheme = theme_light(),
                            risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                            risk.table.y.text = TRUE, conf.int.style = "step",
                            surv.median.line = "hv")

## ****************************
## Survival Analysis - TCGA Filtered
fit_os_tcga_filtered <- survfit(Surv(OS.time.month, as.integer(OS)) ~ Signature,
                                data = mat_tcga_tpm_transposed_merged_filtered)
os_plot_tcga_filtered <- ggsurvplot(fit_os_tcga_filtered, data = mat_tcga_tpm_transposed_merged_filtered,
                                     risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                                     palette = c("#E74C3C", "#3498DB", "#F39C12"),
                                     xlim = c(0, max_months),
                                     title = "OS TCGA",
                                     subtitle = paste0("Top Genes:", length(genes_set$Gene),
                                                       ", t:", t, " - Stage III/IV"),
                                     xlab = "Time (Months)",  ylab = "OS Probability",
                                     break.time.by = 4, ggtheme = theme_light(),
                                     risk.table.y.text.col = TRUE, risk.table.height = 0.25,
                                     risk.table.y.text = TRUE, conf.int.style = "step",
                                     surv.median.line = "hv")

## Plot everything together
blank <- grid::nullGrob()
grid.arrange(p_os_status_distribution, os_plot_own$plot, p_pfs_status_distribution, pfs_plot_own$plot,
             p_os_status_distribution_tcga, os_plot_tcga$plot, p_os_status_distribution_filtered_tcga, os_plot_tcga_filtered$plot, 
             ncol = 4)



















## ************************************************************
## EDA

## Boxplot of expression conditioning on all stages
names(df_clinical_tcga)[names(df_clinical_tcga) == "sample"] <- "ID"
df_clinical_tcga_plot <- df_clinical_tcga %>%  
  filter(ID %in% selected_colnames_tcga)
df_clinical_tcga_plot$ID <- sub("-[^-]+$", "", df_clinical_tcga_plot$ID)

## Merge
df_tcga_boxplot_stages <- merge(mat_tcga_tpm_transposed_merged, df_clinical_tcga_plot, by = "ID")

## Group by
df_plot_stage <- df_tcga_boxplot_stages %>%
  filter(!is.na(ajcc_pathologic_stage.diagnoses)) %>%  
  mutate(ajcc_pathologic_stage.diagnoses = recode(ajcc_pathologic_stage.diagnoses,
                                                  "Stage I"    = "Stage I-II",
                                                  "Stage IA"   = "Stage I-II",
                                                  "Stage IB"   = "Stage I-II",
                                                  "Stage II"   = "Stage I-II",
                                                  "Stage IIA"  = "Stage I-II",
                                                  "Stage IIB"  = "Stage I-II",
                                                  "Stage III"  = "Stage III-IV",
                                                  "Stage IIIA" = "Stage III-IV",
                                                  "Stage IIIB" = "Stage III-IV",
                                                  "Stage IIIC" = "Stage III-IV",
                                                  "Stage IV"   = "Stage III-IV"))

## Boxplot Stage I-II Vs. Stage III-Iv
custom_colors <- c("Stage I-II"    = "steelblue", 
                   "Stage III-IV"  = "red")
ggplot(df_plot_stage, aes(x = ajcc_pathologic_stage.diagnoses, y = Genes, fill = ajcc_pathologic_stage.diagnoses)) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.7) +
  labs(title = "Genes Expression By Stage", x = "Stage", y = "Genes", fill = "Stage",
       subtitle = "TCGA") +
  scale_fill_manual(values = custom_colors) + 
  theme_minimal() +
  theme(legend.position = "right") +
  stat_compare_means(method = "kruskal.test", label.y = 2.5)

## Molecular Subtype TCGA
## Load clinical from bioportal where we consider the subtype variable
df_clinical_bioportal <- read_tsv(paste0(data_path, "/stad_tcga_pan_can_atlas_2018_clinical_data.tsv"),
                                  show_col_types = FALSE)

df_clinical_bioportal_filtered <- df_clinical_bioportal %>%  
  filter(`Patient ID` %in% mat_tcga_tpm_transposed_merged$ID)
names(df_clinical_bioportal_filtered)[names(df_clinical_bioportal_filtered) == "Patient ID"] <- "ID"
df_clinical_bioportal_filtered <- merge(mat_tcga_tpm_transposed_merged, df_clinical_bioportal_filtered, by = "ID")

## Check
table(df_clinical_bioportal_filtered$Subtype)

## Filter NA
df_clinical_bioportal_filtered_plot <- df_clinical_bioportal_filtered %>%
  filter(!is.na(Subtype))

## Tests and Plots
## Calculate percentage within each Signature category
df_percent <- df_clinical_bioportal_filtered_plot %>%
  group_by(Signature) %>%
  count(Subtype) %>%
  mutate(Percentage = n / sum(n))

## Barplot Molecular Subtypes with p valu and position fill
ggplot(df_percent, aes(x = Signature, y = Percentage, fill = Subtype)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Molecular Subtype By Signature", 
       x = "Signature", y = "Percentage", fill = "Subtype") + 
  scale_y_continuous(labels = scales::percent_format()) +  
  scale_fill_manual(values = c("STAD_CIN" = "steelblue", "STAD_EBV" = "red", "STAD_GS" = "green", 
                               "STAD_MSI" = "purple", "STAD_POLE" = "orange"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +  
  theme_minimal() +  
  theme(legend.position = "right")

## Conditioning
## Remove NA values for relevant columns
df_clinical_bioportal_filtered_plot_no_na <- df_clinical_bioportal_filtered_plot %>%
  filter(!is.na(Signature), !is.na(Subtype), !is.na(Genes))

## Define pairwise comparisons
my_comparison <- list(c("STAD_CIN", "STAD_EBV"), c("STAD_CIN", "STAD_GS"),
                      c("STAD_CIN", "STAD_MSI"), c("STAD_CIN", "STAD_POLE"),
                      c("STAD_EBV", "STAD_GS"), c("STAD_EBV", "STAD_MSI"),
                      c("STAD_EBV", "STAD_POLE"), c("STAD_GS", "STAD_MSI"),
                      c("STAD_GS", "STAD_POLE"),  c("STAD_MSI", "STAD_POLE"))
## Plot
plot_genes_signature <- ggplot(df_clinical_bioportal_filtered_plot_no_na, aes(x = Signature, y = Genes, fill = Subtype)) +
  geom_boxplot() + 
  labs(title = "Genes Expression by Signature", 
       x = "Signature", y = "Genes", fill = "Subtype") + 
  scale_fill_manual(values = c("STAD_CIN" = "steelblue", "STAD_EBV" = "red", "STAD_GS" = "green", 
                               "STAD_MSI" = "purple", "STAD_POLE" = "orange"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +  
  theme_minimal() + 
  theme(legend.position = "right") 

## Plot without conditioning
ggplot(df_clinical_bioportal_filtered_plot_no_na, aes(x = Subtype, y = Genes, fill = Subtype)) +
  geom_boxplot() + 
  labs(title = "Genes Expression by Subtype", 
       x = "Subtype", y = "Genes", fill = "Subtype") + 
  scale_fill_manual(values = c("STAD_CIN" = "steelblue", "STAD_EBV" = "red", "STAD_GS" = "green", 
                               "STAD_MSI" = "purple", "STAD_POLE" = "orange"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +  
  theme_minimal() + 
  theme(legend.position = "right") +
  stat_compare_means(method = "kruskal.test", label = "p.format")


## Plot without signature condition
df_percent <- df_clean %>%
  mutate(Subtype_Signature = paste(Subtype, Signature, sep = " "))
df_percent$Subtype_Signature <- factor(df_percent$Subtype_Signature, 
                                       levels = c("STAD_CIN High", "STAD_CIN Low", 
                                                  "STAD_EBV High", "STAD_EBV Low", 
                                                  "STAD_GS High", "STAD_GS Low", 
                                                  "STAD_MSI High", "STAD_MSI Low", 
                                                  "STAD_POLE High", "STAD_POLE Low"))
## Plot
my_comparisons<- list(c("STAD_CIN High", "STAD_CIN Low"),
                      c("STAD_EBV High", "STAD_EBV Low"),
                      c("STAD_GS High", "STAD_GS Low"),
                      c("STAD_MSI High", "STAD_MSI Low"),
                      c("STAD_POLE High", "STAD_POLE Low"))
ggplot(df_percent, aes(x = Subtype_Signature, y = Genes, fill = Subtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Molecular Subtype By Signature",
       x = "Signature", y = "Expression", fill = "Subtype") +
  scale_fill_manual(values = c("STAD_CIN" = "steelblue", "STAD_EBV" = "red", "STAD_GS" = "green",
                               "STAD_MSI" = "purple", "STAD_POLE" = "orange"),
                    labels = c("CIN", "EBV", "GS", "MSI", "POLE")) +
  theme_minimal() +
  theme(legend.position = "right") +
  stat_compare_means(comparisons = my_comparisons, bracket.size = .6, size = 4) +
  stat_compare_means(label.y = 4.5, method = "kruskal.test") + theme(legend.position = "none")
  
## Barplot
my_comparisons<- list(c("STAD_CIN High", "STAD_CIN Low"),
                      c("STAD_EBV High", "STAD_EBV Low"),
                      c("STAD_GS High", "STAD_GS Low"),
                      c("STAD_MSI High", "STAD_MSI Low"),
                      c("STAD_POLE High", "STAD_POLE Low"))
ggplot(df_percent, aes(x = Subtype_Signature, y = after_stat(count), fill = Subtype_Signature)) +
  geom_bar() + 
  labs(title = "Molecular Subtype By Signature", 
       x = "Signature", y = "Count", fill = "Subtype") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  
  scale_fill_manual(values = c("STAD_CIN High" = "steelblue", "STAD_CIN Low" = "steelblue",
                               "STAD_EBV High" = "red", "STAD_EBV Low" = "red",
                               "STAD_GS High" = "green", "STAD_GS Low" = "green", 
                               "STAD_MSI High" = "purple", "STAD_MSI Low" = "purple",
                               "STAD_POLE High" = "orange", "STAD_POLE Low" = "orange")) +  
  theme_minimal() +  
  theme(legend.position = "right") +
  stat_compare_means(comparisons = my_comparisons, bracket.size = .6, size = 4) +
  stat_compare_means(label.y = 4.5, method = "kruskal.test") + theme(legend.position = "none")


stat_compare_means(df_percent, comparisons = my_comparisons, bracket.size = .6, size = 4)


