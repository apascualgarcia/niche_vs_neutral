library(tidyverse)
library(nleqslv)
library(fitdistrplus)
library(fitteR)
library(flexsurv)

library(goftest)

neg_log_lik_trunc_norm <- function(params, data, cutoff) {
  mu <- params[1]
  sigma <- params[2]
  
  if (sigma <= 0) return(Inf)
  
  
  densities <- dnorm(data, mean = mu, sd = sigma)
  trunc_factor <- 1 - pnorm(cutoff, mean = mu, sd = sigma)
  
  ll <- sum(log(densities / trunc_factor))
  
  return(-ll) 
}

dgamma_shifted <- function(x, shape = 6.25, rate = 2.5, shift = -2.5) {
  # The Gamma distribution only exists to the right of the 'shift'
  # If x is smaller than the shift, the density is 0
  ifelse(x > shift, dgamma(x - shift, shape = shape, rate = rate), 0)
}


prepare_data <- function(data_meta, Meta)
{
  if(is.matrix(data_meta))
    rescaled_data <- as.data.frame(t((t(data_meta) - colMeans(data_meta))/apply(data_meta,FUN = sd, MARGIN = 2)))
  else{
    rescaled_data <- as.data.frame(t(as.data.frame((log_data_meta - mean(log_data_meta))/sd(log_data_meta))))
    rownames(rescaled_data)[1] = "Sam1"
    }
  
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
  
  long_data$Metaloc <- Meta
  
  return(long_data)
}
# Helper function to compute the Vuong test without repeating code
calc_vuong <- function(ll_1, ll_2, name_1, name_2) {
  diff_ll <- ll_1 - ll_2
  vuong_stat <- mean(diff_ll) / (sd(diff_ll) / sqrt(n))
  p_val <- 2 * pnorm(-abs(vuong_stat))
  
  cat("--- Vuong Test:", name_1, "vs", name_2, "---\n")
  cat("Statistic:", vuong_stat, "\n")
  cat("P-value:", p_val, "\n")
  if(p_val < 0.05) {
    cat("Winner:", ifelse(vuong_stat > 0, name_1, name_2), "\n\n")
  } else {
    cat("Result: statistical tie (neither is significantly better)\n\n")
  }
}


data <- t(read.table(file = "../otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = "../Metadata_singleTree_wDirich-K10.tsv", header = T, row.names = 1, sep = "\t")
data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))
meta <- meta %>% filter(Location %in% names(table(meta$Location)[table(meta$Location) > 1]))

#### Then, I will represent the distribution of the metacommunities with the tree holes distributions


log_all_analysis <- data.frame()

for (i in unique(meta$Dirichlet[!is.na(meta$Dirichlet)]))
{
  meta_Dir = meta[meta$Dirichlet == i & !is.na(meta$Dirichlet), ]

  TH = unique(meta_Dir$Location)

  # TH = TH[  table(meta_Dir$Location) > 1]

  data_meta <- data_rel[rownames(meta_Dir),]

  if (nrow(data_meta) >= 2){
    data_meta <- data_meta[,colSums(data_meta > 0) == nrow(data_meta)]
  } else{
    data_meta <- data_meta[data_meta > 0]}

  for (j in TH)
  {
    Meta_Dir_TH <- meta_Dir[meta_Dir$Location == j & !is.na(meta_Dir$Location), ]

    data_meta_TH <- data_meta[rownames(Meta_Dir_TH),]

    log_data_meta <- log(data_meta_TH)

    log_long_data <- prepare_data(log_data_meta, i)

    log_long_data$TH = j

    log_all_analysis <- rbind(log_all_analysis, log_long_data)

  }
}

### Now with the data ordered I try to do the figures and the analysis:

indval <- read.csv("indval_M7.csv", header = T, row.names = 1)

for (i in unique(log_all_analysis$Metaloc))
{
  data_metaCom <- log_all_analysis %>% filter(Metaloc == i)
  
  z <- data_metaCom$Value
  shift <- abs(min(z)) + 0.001 
  z_shifted <- z + shift
  
  x_vals <- seq(min(z), max(z), length.out = 1000)
  
  data_metaCom$TH <- factor(data_metaCom$TH, 
                            levels = names(sort(table(data_metaCom$TH), decreasing = TRUE)))
  
  
  data_metaCom$indval <- indval[data_metaCom$ASVs, "Location"]
  
  N <- nrow(data_metaCom)
  
  dens_data <- density(data_metaCom$Value)
  df_real_points <- data.frame(x = dens_data$x, y = dens_data$y)
  step <- ceiling(nrow(df_real_points) / 50) # Change the '30' to get more or fewer points
  df_real_points <- df_real_points[seq(1, nrow(df_real_points), by = step), ]
  
  dens_real <- density(data_metaCom$Value)
  df_dens_real <- data.frame(x = dens_real$x, y = dens_real$y)
  
  step <- ceiling(nrow(df_dens_real) / 50) 
  df_trajectory_points <- df_dens_real[seq(1, nrow(df_dens_real), by = step), ]
  
 
  p <- ggplot() +
    # A. Global empirical distribution (grey background)
    geom_area(data = df_dens_real, aes(x = x, y = y), 
              fill = "grey90") +
    geom_line(data = df_dens_real, aes(x = x, y = y), 
              color = "darkgrey", linetype = "dashed") +
    
    # B. Spaced points
    geom_point(data = df_trajectory_points, aes(x = x, y = y), 
               color = "black", size = 2) +
    
    # C. Distributions per combination of TH and indval (NAs omitted)
    geom_density(data = data_metaCom %>% filter(!is.na(indval)),
                 aes(x = Value, 
                     y = after_stat(count) / N, 
                     fill = paste(TH, indval, sep = "-")), 
                 alpha = 0.5, 
                 color = "black") +
    
    labs(
      title = paste("Metalocation Distributions -", i, "  Samples:   ", 
                    paste(paste("Hole", names(table((meta %>% filter(Dirichlet == i))$Location))), 
                          table((meta %>% filter(Dirichlet == i))$Location), sep = ": ", collapse = ", ")),
      x = "Transformed log relative abundance",
      y = "Density",
      color = "Models",           # Legend title of the theoretical lines
      fill = "TH - Indval"        # Updated to reflect the new combination
    ) +
    theme_bw()
  
  print(p)

}


