## Load libraries
library(tidyr)
library(dplyr)
library(readxl)
library(survival)
library(broom)
library(forestplot)
library(survminer)

## Set data path directory
data_path   <- "C:/Users/david/Documents/IFO/Gastric/Data"
output_path <- "C:/Users/david/Documents/IFO/Gastric/Output_Res"

## ***************************************
## LOAD DATA
## ***************************************
df_prospettici   <- read_excel(paste0(data_path, "/Prospettici_DB_aggiornato_cate.xlsx"))
df_retrospettivi <- read_excel(paste0(data_path, "/Retrospettivi_DB_aggiornato_cate.xlsx"))
## Consider the final dataset with signature
df_signature      <- read_excel(paste0(data_path, "/df_signature.xlsx"))
df_signature_os   <- read_excel(paste0(data_path, "/df_signature_os_filtered.xlsx"))
df_signature_pfs  <- read_excel(paste0(data_path, "/df_signature_pfs_filtered.xlsx"))

## ***************************************
## WRANGLING
## ***************************************
## Change column names format
colnames(df_prospettici)   <- gsub(" ", "_", colnames(df_prospettici))
colnames(df_retrospettivi) <- gsub(" ", "_", colnames(df_retrospettivi))
## PROSPETTICI
## Change name columns
names(df_prospettici)[names(df_prospettici) == "Sample"] <- "ID"
names(df_prospettici)[names(df_prospettici) == "PAZIENTE_codice"] <- "NOME_COGNOME"
names(df_prospettici)[names(df_prospettici) == "Performance_Status_ECOG"] <- "PS"
names(df_prospettici)[names(df_prospettici) == "Metastatico_LocAvanzato_operabile"] <- "MLA"
names(df_prospettici)[names(df_prospettici) == "Intervento_Chirurgico_si_no"] <- "SURGERY"
names(df_prospettici)[names(df_prospettici) == "N_Metastasi"] <- "METASTASIS_NUMBER"
names(df_prospettici)[names(df_prospettici) == "Sesso"] <- "SEX"
## Convert datetime variables and create age
df_prospettici$Data_diagnosi <- as.Date(df_prospettici$Data_diagnosi, format = "%d/%m/%Y")
df_prospettici$Data_nascita  <- as.Date(df_prospettici$Data_nascita, format = "%d/%m/%Y")
## Age variable
df_prospettici$AGE <- as.numeric(df_prospettici$Data_diagnosi - df_prospettici$Data_nascita)
df_prospettici$AGE <- df_prospettici$AGE / 365
## Populate Age variable
ids_age <- c("IRE-190", "IRE-191", "IRE-192", "IRE-193", "IRE-194",
             "IRE-195", "IRE-196", "IRE-197", "IRE-198", "IRE-199",
             "IRE-200", "IRE-256", "IRE-257", "IRE-258", "IRE-259",
             "IRE-260", "IRE-261", "IRE-262", "IRE-263", "IRE-289",
             "IRE-290")
ids_age_replace <- c(63, 71, 45, 55, 62,
                     69, 56, 73, 49, 58,
                     66, 65, 69, 53, 59,
                     68, 71, 73, 68, 67,
                     71)
df_support      <- data.frame(ID = ids_age, AGE = ids_age_replace)
matched_indices <- match(df_prospettici$ID, ids_age)
df_prospettici$AGE[!is.na(matched_indices)] <- ids_age_replace[matched_indices[!is.na(matched_indices)]]
df_prospettici$AGE <- ifelse(df_prospettici$AGE >
                               quantile(df_prospettici$AGE, na.rm = T, probs = c(0.5)),
                             "High", "Low")
## Metastatico Localmente Avanzato
df_prospettici$MLA <- ifelse(df_prospettici$MLA == 2, "LOCALLY ADVANCED", "METASTATIC")
## Numero Metastasi
df_prospettici$METASTASIS_NUMBER <- ifelse(df_prospettici$METASTASIS_NUMBER == 0, "1 SITE", "+ 1 SITE")

## Intervento Chirurgico
df_prospettici$SURGERY <- ifelse(df_prospettici$SURGERY == 1, "SI", "NO")

## RETROSPETTIVI
## Change name columns
names(df_retrospettivi)[names(df_retrospettivi) == "ID_DNA"] <- "ID"
names(df_retrospettivi)[names(df_retrospettivi) == "Nome_e_Cognome"] <- "NOME_COGNOME"
names(df_retrospettivi)[names(df_retrospettivi) == "Performance_Status_ECOG"] <- "PS"
names(df_retrospettivi)[names(df_retrospettivi) == "Int_Chir"] <- "SURGERY"
names(df_retrospettivi)[names(df_retrospettivi) == "I_linea_CT"] <- "MLA"
names(df_retrospettivi)[names(df_retrospettivi) == "Sesso"] <- "SEX"
## Convert datetime variables and create age
df_retrospettivi$`Data diagn`       <- as.Date(df_retrospettivi$Data_diagn, format = "%Y-%m-%d")
df_retrospettivi$`Data di nascita`  <- as.Date(df_retrospettivi$Data_di_nascita, format = "%Y-%m-%d")
## Age variable
df_retrospettivi$AGE <- as.numeric(df_retrospettivi$Data_diagn - df_retrospettivi$Data_di_nascita)
df_retrospettivi$AGE <- df_retrospettivi$AGE / 365
df_retrospettivi$AGE <- ifelse(df_retrospettivi$AGE >
                                 quantile(df_retrospettivi$AGE, na.rm = T, probs = c(0.5)),
                               "High", "Low")
## Metastasis Number
df_retrospettivi$METASTASIS_NUMBER <- ifelse(grepl("\\+|,", df_retrospettivi$Sede_M), "+ 1 SITE", "1 SITE")
## Intervento Chirurgico
df_retrospettivi$SURGERY <- ifelse(df_retrospettivi$SURGERY == 1, "SI", "NO")
## MLA
df_retrospettivi$MLA <- ifelse(df_retrospettivi$MLA == 2, "LOCALLY ADVANCED ", "METASTATIC")
## Select columns
columns_clinical  <- c("ID", "NOME_COGNOME", "AGE", "SEX",
                       "PS", "MLA", "METASTASIS_NUMBER", "SURGERY")

columns_signature <- c("ID", "Genes", "Signature", "OS_time", "OS_Event",
                       "PFS_time", "PFS_event", "Coorte")

## Filter dataset
## PROSPETTICI
df_prospettici_filtered <- df_prospettici %>% 
  dplyr::select(all_of(columns_clinical))
## RETROSPETTIVI
df_retrospettivi_filtered <- df_retrospettivi %>% 
  dplyr::select(all_of(columns_clinical))
## Filter signature
## OS
df_signature_os_filtered  <- df_signature_os %>% 
  dplyr::select(all_of(columns_signature))
## PFS
df_signature_pfs_filtered <- df_signature_pfs %>% 
  dplyr::select(all_of(columns_signature))
## Merge
df_clinical_dataset <- rbind(df_prospettici_filtered, df_retrospettivi_filtered)
## OS Signature
df_final_os <- df_signature_os_filtered %>%
  merge(df_clinical_dataset, by = "ID")
## PFS Signature
df_final_pfs <- df_signature_pfs_filtered %>%
  merge(df_clinical_dataset, by = "ID")


## ## ***************************************
## MULTIVARIATE ANALYSIS
## ***************************************
covariate_names <- c(AGE = "AGE", PS = "PS", MLA = "MLA", SEX = "SEX",
                     METASTASIS_NUMBER = "METASTASIS_NUMBER", SURGERY = "SURGERY")
## Reveal
## OS
df_final_os$PS  <- relevel(factor(df_final_os$PS), ref = "0")
df_final_os$MLA <- relevel(factor(df_final_os$MLA), ref = "LOCALLY ADVANCED")
df_final_os$METASTASIS_NUMBER <- relevel(factor(df_final_os$METASTASIS_NUMBER), ref = "1 SITE")
df_final_os$SURGERY <- relevel(factor(df_final_os$SURGERY), ref = "NO")
df_final_os$SIGNATURE <- relevel(factor(df_final_os$Signature), ref = "Low")
df_final_os$AGE <- relevel(factor(df_final_os$AGE), ref = "Low")
df_final_os$SEX <- relevel(factor(df_final_os$SEX), ref = "M")
## PFS
df_final_pfs$PS  <- relevel(factor(df_final_pfs$PS), ref = "0")
df_final_pfs$MLA <- relevel(factor(df_final_pfs$MLA), ref = "LOCALLY ADVANCED")
df_final_pfs$METASTASIS_NUMBER <- relevel(factor(df_final_pfs$METASTASIS_NUMBER), ref = "1 SITE")
df_final_pfs$SURGERY <- relevel(factor(df_final_pfs$SURGERY), ref = "NO")
df_final_pfs$SIGNATURE <- relevel(factor(df_final_pfs$Signature), ref = "Low")
df_final_pfs$AGE <- relevel(factor(df_final_pfs$AGE), ref = "Low")
df_final_pfs$SEX <- relevel(factor(df_final_pfs$SEX), ref = "M")

## Regression
## OS
reg_os  <- coxph(Surv(OS_time) ~ AGE + PS + SEX + METASTASIS_NUMBER + MLA + SIGNATURE,
                 data = df_final_os)
## PFS
reg_pfs  <- coxph(Surv(PFS_time) ~ AGE + PS + SEX + METASTASIS_NUMBER + MLA + SIGNATURE,
                  data = df_final_pfs)

## Forest Plot
## OS
ggforest(reg_os, df_final_os, main = "Overall Survival Multivariate Analysis")
## PFS
ggforest(reg_pfs, df_final_pfs, main = "Progression Free Survival Multivariate Analysis")
