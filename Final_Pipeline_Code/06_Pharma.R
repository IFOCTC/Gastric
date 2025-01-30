## Load Useful
source("00_Support.R")

## Load data
df_farmaci <- read_excel(paste0(data_path,
                                "/farmaci_20nM.xlsx"))

## Wrangling
ids_resistant <- c("GTR459",  "GTR00498", "GTR00508")
ids_sensitive <- c("GTR00210", "GTR0042", "GTR125", "GTR221", "GTR607")
columns <- c("FARMACO", ids_resistant, ids_sensitive)

## Create Condition Column
df_farmaci_filtered <- df_farmaci %>% 
  dplyr::select(all_of(columns))

## Reshape data
df_farmaci_filtered_long <- df_farmaci_filtered %>%
  pivot_longer(cols = all_of(c(ids_resistant, ids_sensitive)),  
               names_to = "Variable", 
               values_to = "Value") %>%
  mutate(Group = ifelse(Variable %in% ids_resistant, "Resistant", "Sensitive"))


## Perform Kruskal-Wallis test for each unique FARMACO
res_test <- lapply(unique(df_farmaci_filtered_long$FARMACO), function(Farmaco) {
  farmaco_data <- subset(df_farmaci_filtered_long, FARMACO == Farmaco)
  test_result <- kruskal.test(Value ~ Group, data = farmaco_data)
  list(Farmaco = Farmaco, p_value = test_result$p.value)})

## Combine res
res_test_df <- do.call(rbind,
                       lapply(res_test,
                       function(res_test) data.frame(FARMACO = res_test$Farmaco,
                                                     p_value = res_test$p_value)))
## P-value <= 0.1
res_test_df  <- res_test_df %>%
  filter(p_value <= 0.1)

## Merge
df_farmaci_filtered_long_plot <- df_farmaci_filtered_long %>%
  filter(FARMACO %in%  res_test_df$FARMACO)

## Plot
drugs <- unique(df_farmaci_filtered_long$FARMACO)
for (drug in drugs) {
  df_drug <- subset(df_farmaci_filtered_long, FARMACO == drug)
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
    theme_minimal(base_size = 14) +
    theme(legend.position = "",  
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  
          axis.text.y = element_text(size = 12, face = "plain"), 
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
          axis.title = element_text(size = 14, face = "bold"), 
          panel.grid.major = element_line(size = 0.2, color = "gray90"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = "white"))
    #stat_compare_means(method = "kruskal.test", label.x = 1.5)}
  ggsave(paste0("C:/Users/david/Documents/IFO/Final_Pipeline_Code/Output_Images/Boxplot_Drug/",
                drug, "_boxplot.png"), plot = p, width = 8, height = 6)}



































