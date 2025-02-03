## Load Useful
source("00_Support.R")

## Load data
df_farmaci_20  <- read_excel(paste0(data_path,
                                    "/farmaci_20nM.xlsx"))
df_farmaci_200 <- read_excel(paste0(data_path,
                                   "/farmaci_200nM.xlsx"))
## Wrangling
ids_resistant <- c("GTR459",  "GTR00498", "GTR00508")
ids_sensitive <- c("GTR00210", "GTR0042", "GTR125", "GTR221", "GTR607")
columns <- c("FARMACO", ids_resistant, ids_sensitive)

## Create Condition Column
## 20
df_farmaci_filtered_20  <- df_farmaci_20 %>% 
  dplyr::select(all_of(columns))
## 200
df_farmaci_filtered_200 <- df_farmaci_200 %>% 
  dplyr::select(all_of(columns))

## Reshape data
## 20
df_farmaci_filtered_long_20  <- df_farmaci_filtered_20 %>%
  tidyr::pivot_longer(cols = all_of(c(ids_resistant, ids_sensitive)),  
                      names_to = "Variable", 
                      values_to = "Value") %>%
  mutate(Group = ifelse(Variable %in% ids_resistant, "Resistant", "Sensitive"))
## 200
df_farmaci_filtered_long_200 <- df_farmaci_filtered_200 %>%
  tidyr::pivot_longer(cols = all_of(c(ids_resistant, ids_sensitive)),  
                      names_to = "Variable", 
                      values_to = "Value") %>%
  mutate(Group = ifelse(Variable %in% ids_resistant, "Resistant", "Sensitive"))

## Plot
## 20
drugs_20 <- unique(df_farmaci_filtered_long_20$FARMACO)
for (drug in drugs_20) {
  df_drug <- subset(df_farmaci_filtered_long_20, FARMACO == drug)
  p <- ggplot(df_drug, aes(x = factor(Group, levels = c("Sensitive", "Resistant")),
                           y = Value, fill = Group)) +
    geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
                 colour = "black", alpha = 0.8) +  
    labs(title = paste(drug, "Drug Distribution by Response"),
         x = "Response",
         y = "Alive Cellular Concentration [%]",
         fill = "Response") +
    scale_fill_manual(values = c("Sensitive" = "#1f77b4",  
                                 "Resistant" = "#d62728")) + 
    ylim(0, max(df_farmaci_filtered_long_20$Value)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "",  
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
          axis.text.y = element_text(size = 12, face = "plain"), 
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
          axis.title = element_text(size = 14, face = "bold"), 
          panel.grid.major = element_line(size = 0.2, color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    stat_compare_means(method = "kruskal.test", label.x = 1.5) 
  ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Boxplot_Drug_20/",
                drug, "_boxplot_20.png"), plot = p, width = 8, height = 6)}

## 200
drugs_200 <- unique(df_farmaci_filtered_long_200$FARMACO)
for (drug in drugs_200) {
  df_drug <- subset(df_farmaci_filtered_long_200, FARMACO == drug)
  p <- ggplot(df_drug, aes(x = factor(Group, levels = c("Sensitive", "Resistant")),
                           y = Value, fill = Group)) +
    geom_boxplot(outlier.size = 1, outlier.colour = "black", width = 0.7, 
                 colour = "black", alpha = 0.8) +  
    labs(title = paste(drug, "Drug Distribution by Response"),
         x = "Response",
         y = "Alive Cellular Concentration [%]",
         fill = "Response") +
    scale_fill_manual(values = c("Sensitive" = "#1f77b4",  
                                 "Resistant" = "#d62728")) + 
    ylim(0, max(df_farmaci_filtered_long_200$Value)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "",  
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
          axis.text.y = element_text(size = 12, face = "plain"), 
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
          axis.title = element_text(size = 14, face = "bold"), 
          panel.grid.major = element_line(size = 0.2, color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white")) +
    stat_compare_means(method = "kruskal.test", label.x = 1.5) 
  ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Boxplot_Drug_200/",
                drug, "_boxplot_200.png"), plot = p, width = 8, height = 6)}
