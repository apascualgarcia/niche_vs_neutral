library(tidyverse)


#### HERE I PLOT THE P-P PLOT DIFFERENCE AGAINST THE ASSOCIATIONS ########


to_matrix <- function(edgelist)
{
  names = unique(c(edgelist$taxa_1, edgelist$taxa_2))
  mat <- matrix(0, nrow = length(names), ncol = length(names))
  colnames(mat) = rownames(mat) = names
  
  for (i in 1:nrow(edgelist))
    mat[edgelist[i,"taxa_1"], edgelist[i,"taxa_2"]] = edgelist[i,"value"]
  
  mat = mat + t(mat)
  return(mat)
}




cors <- read.table("../6.Correlation_Distribution/Correlations.tsv", header = T, sep = "\t")
meta <- read.table("../5.Real_data/Single_tree/Metadata_singleTree_wDirich-K10.tsv", header = T, sep = "\t", row.names = 1)
data <- t(read.table(file = "../6.Real_data/Single_tree/otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, row.names = 1, sep = "\t"))

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))



plot_data <- data.frame()
all_analysis <- data.frame()
for (i in c(1,2,3,4,5,6,7,8,10))
{
  cors_red <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i) & abs(value) > 0.25)
  
  cors_red_mat <- to_matrix(cors_red)
  
  # NOTE: the original code used a Top 30, even though the variables are called 'top_20'. The Top 30 is kept.
  top_30_asvs <- names(sort(abs(colSums(cors_red_mat)), decreasing = T))[1:30]
  
  data_meta <- data_rel[rownames(meta[meta$Dirichlet == i & !is.na(meta$Dirichlet), ]),]
  
  if (nrow(data_meta) > 2){
    data_meta <- data_meta[,colSums(data_meta > 0) == nrow(data_meta)]
  } else{
    data_meta <- data_meta[data_meta > 0]}
  
  rescaled_data <- as.data.frame(t((t(data_meta) - colMeans(data_meta))/apply(data_meta,FUN = sd, MARGIN = 2)))
  rescaled_data <- as.data.frame(rescaled_data)
  
  long_data <- rescaled_data %>% 
    rownames_to_column(var = "Samples") %>%  
    pivot_longer(
      cols = -Samples,
      names_to = "ASVs",
      values_to = "Value"
    )
  
  long_data_data <- as.data.frame(data_meta) %>% 
    rownames_to_column(var = "Samples") %>%  
    pivot_longer(
      cols = -Samples,
      names_to = "ASVs",
      values_to = "Value"
    )
  
  long_data$Rel_abun = long_data_data$Value
  
  long_data$Metaloc <- i
  
  all_analysis <- rbind(all_analysis, long_data)
  
  data_metaCom_all <- all_analysis %>% filter(Metaloc == i)
  

  
  z <- data_metaCom_all$Value
  shift <- abs(min(z)) + 0.001 
  z_shifted <- z + shift
  
  x_vals <- seq(min(z), max(z), length.out = 1000)
  # For gamma
  shape_val <- mean(z + abs(min(z)))^2 / var(z) 
  rate_val <- shape_val / mean(z + abs(min(z)))
  
  # For lognormal
  meanlog_val <- mean(log(z_shifted))
  sdlog_val <- sd(log(z_shifted))
  
  # --- MAIN CHANGE 1: data preparation and reordering for the Z-order ---
  data_metaCom <- data_metaCom_all %>%
    mutate(
      Shifted_Value = Value + shift,
      # The original variable name is kept for consistency, although it is a Top 30
      is_top20 = ifelse(ASVs %in% top_30_asvs, "Top", "Others")
    ) %>%
    arrange(Shifted_Value) %>% 
    mutate(
      emp_prob = ppoints(n()), 
      theo_prob_gamma = pgamma(Shifted_Value, shape = shape_val, rate = rate_val),
      theo_prob_LN = plnorm(Shifted_Value, meanlog = meanlog_val, sdlog = sdlog_val)
    ) %>%
    # Crucial! Reorder so that the 'Top' points are drawn LAST (on top)
    arrange(desc(is_top20)) 
  
  top_20_abun <- data_metaCom %>%
    group_by(ASVs) %>%
    summarise(Total_Rel_abun = sum(Rel_abun, na.rm = TRUE)) %>%
    slice_max(order_by = Total_Rel_abun, n = 20) %>%
    pull(ASVs)
  
  no_top_20_asvs <- data_metaCom %>%
    filter(! ASVs %in% top_30_asvs)
  top_20_asvs <- data_metaCom %>%
    filter(ASVs %in% top_30_asvs)
  top_20_asvs_abun <- data_metaCom %>%
    filter(ASVs %in% top_20_abun)
  no_top_20_asvs_abun <- data_metaCom %>%
    filter(!ASVs %in% top_20_abun)
  
  plot_data <- data_metaCom
  plot_data <- no_top_20_asvs
  plot_data <- top_20_asvs
  plot_data <- top_20_asvs_abun
  plot_data <- no_top_20_asvs_abun
  
  # --- MAIN CHANGE 2: P-P plot for the Gamma ---
  pp_plot_gamma <- ggplot(plot_data, aes(x = theo_prob_gamma, y = emp_prob, color = is_top20)) +
    # alpha is added as an aesthetic (aes), not as a fixed value
    geom_point(aes(alpha = is_top20), size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
    
    # Define the colours (changed to 'Top' to match)
    scale_color_manual(values = c("Top" = "red", "Others" = "grey80")) +
    
    # --- Manual transparencies ---
    # 'Others' gets a very low alpha (0.1) and 'Top' a high alpha (1)
    scale_alpha_manual(values = c("Top" = 1, "Others" = 0.1)) +
    
    # Remove the alpha legend so it is not confusing; keep only the colour legend
    guides(alpha = "none") +
    
    labs(
      title = paste("P-P Plot Gamma -", i),
      x = "Theoretical probabilities",
      y = "Empirical probabilities",
      color = "ASV Status"
    ) +
    theme_bw()
  
  
  # --- MAIN CHANGE 3: P-P plot for the Lognormal ---
  pp_plot_lnorm <- ggplot(plot_data, aes(x = theo_prob_LN, y = emp_prob, color = is_top20)) +
    # alpha is added as an aesthetic (aes)
    geom_point(aes(alpha = is_top20), size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
    
    # Define the colours (changed to 'Top')
    scale_color_manual(values = c("Top" = "#377EB8", "Others" = "grey85")) +
    
    # --- Manual transparencies ---
    # Change these values to make the grey points more or less visible
    scale_alpha_manual(values = c("Top" = 1, "Others" = 0.15)) +
    
    # Remove the alpha legend
    guides(alpha = "none") +
    
    labs(
      title = paste("P-P Plot Lognormal -", i),
      x = "Theoretical probabilities",
      y = "Empirical probabilities",
      color = "ASV Status"
    ) +
    theme_bw()
  
  # Show the plots
  print(pp_plot_gamma)
  print(pp_plot_lnorm)
}
