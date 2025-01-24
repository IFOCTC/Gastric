## Code useful to perform clinical EDA

## Load scripts
source("00_Support.R")

## Load Data
df_clinical_tcga <- read_tsv(paste0(output_path, "/clinical_cleaned.tsv"), show_col_types = FALSE)

## Plot some distributions
## Disease Type
ggplot(df_clinical_tcga, aes(x = disease_type)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Disease Type Distribution", x = "", y = "Count") +
  theme_minimal() 

## Race
ggplot(df_clinical_tcga, aes(x = race.demographic)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Race Distribution", x = "", y = "Count") +
  theme_minimal()

## Gender
ggplot(df_clinical_tcga, aes(x = gender.demographic)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Gender Distribution", x = "", y = "Count") +
  theme_minimal()

## Ethnicity 
ggplot(df_clinical_tcga, aes(x = ethnicity.demographic)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Ethnicity Distribution", x = "", y = "Count") +
  theme_minimal()

## Vital Status
ggplot(df_clinical_tcga, aes(x = vital_status.demographic)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Vital Status Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()

## Age At Index
ggplot(df_clinical_tcga, aes(x = age_at_index.demographic)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Age At Index Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()

## Ajcc Stage
ggplot(df_clinical_tcga, aes(x = ajcc_pathologic_stage.diagnoses)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Ajcc Stage Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()

## Ajcc M
ggplot(df_clinical_tcga, aes(x = ajcc_pathologic_m.diagnoses)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Ajcc M Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()

## Origin
ggplot(df_clinical_tcga, aes(x = tissue_or_organ_of_origin.diagnoses)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Origin Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()

## Primary Diagnosis
ggplot(df_clinical_tcga, aes(x = primary_diagnosis.diagnoses)) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Primary Diagnosis Distribution", x = "", y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()
