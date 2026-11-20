library(tidyverse)
library(nleqslv)
library(fitdistrplus)
library(sads)
library(vegan)
library(EnvStats)

dweibull_shifted <- function(x, shape, scale = 1, shift = 0, log = FALSE) {
  dweibull(x - shift, shape = shape, scale = scale, log = log)
}
pweibull_shifted <- function(q, shape, scale = 1, shift = 0, lower.tail = TRUE, log.p = FALSE) {
  pweibull(q - shift, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}
qweibull_shifted <- function(p, shape, scale = 1, shift = 0, lower.tail = TRUE, log.p = FALSE) {
  qweibull(p, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p) + shift
}
rweibull_shifted <- function(n, shape, scale = 1, shift = 0) {
  rweibull(n, shape = shape, scale = scale) + shift
}

neg_log_lik_trunc_norm <- function(params, data, cutoff) {
  mu <- params[1]
  sigma <- params[2]
  if (sigma <= 0) return(Inf)
  densities <- dnorm(data, mean = mu, sd = sigma)
  trunc_factor <- 1 - pnorm(cutoff, mean = mu, sd = sigma)
  
  ll <- sum(log(densities / trunc_factor))
  
  return(-ll) 
}

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

compute_vuong_winner <- function(ll1, ll2, name1, name2, loc) {
  n <- length(ll1)
  diff_ll <- ll1 - ll2
  LR <- sum(diff_ll)
  omega <- sd(diff_ll)
  
  if (omega == 0) {
    z_stat <- 0
    p_val <- 1
  } else {
    z_stat <- LR / (sqrt(n) * omega)
    p_val <- 2 * (1 - pnorm(abs(z_stat)))
  }
  
  # Pick the winner (5% significance)
  if (p_val < 0.05) {
    winner <- ifelse(z_stat > 0, name1, name2)
  } else {
    winner <- "Tie"
  }
  
  data.frame(
    Metaloc = loc,
    Model_1 = name1,
    Model_2 = name2,
    Z_stat = z_stat,
    P_value = p_val,
    Winner = winner
  )
}

dgamma_shifted <- function(x, shape = 6.25, rate = 2.5, shift = -2.5) {
  ifelse(x > shift, dgamma(x - shift, shape = shape, rate = rate), 0)
}

pgamma_shifted <- function(q, shape = 6.25, rate = 2.5, shift = -2.5) {
  ifelse(q > shift, pgamma(q - shift, shape = shape, rate = rate), 0)
}

data <- t(read.table(file = "../otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = "../Metadata_singleTree_wDirich-K10.tsv", header = T, row.names = 1, sep = "\t")

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

data <- data[rownames(meta)[rownames(meta) %in% rownames(data)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

all_analysis <- data.frame()
all_analysis_counts <- data.frame()
for (i in unique(meta$Dirichlet[!is.na(meta$Dirichlet)]))
{
  data_meta <- data_rel[rownames(meta[meta$Dirichlet == i & !is.na(meta$Dirichlet), ]),]
  data_meta <- data_meta[,colSums(data_meta > 0) == nrow(data_meta)]
  data_norel_meta <- data_meta * 100000000
  
  
  #The logarithm of the relative abundance
  log_means <- log(colMeans(data_meta))
  #The rarefied counts
  means <- colMeans(data_norel_meta)

### To estimate the correction parameters
log_c <- min(log_means)
  
  # Simple initial estimates (empirical mean and standard deviation)
  init_mu <- mean(log_means)
  init_sigma <- sd(log_means)
  
  # Optimisation to find the true mu and sigma
  fit <- optim(par = c(init_mu, init_sigma), 
               fn = neg_log_lik_trunc_norm, 
               data = log_means, 
               cutoff = log_c,
               method = "L-BFGS-B", # Method that allows bounds
               lower = c(-Inf, 0.001)) # sigma must be greater than 0
  
  mu_est <- fit$par[1]
  sigma_est <- fit$par[2]
### End of the estimation

  z_scores <- (log_means - mu_est) / sigma_est

  point_densities <- approx(density(z_scores)$x, density(z_scores)$y, xout = z_scores)$y
  
  long_data <- as.data.frame(t(t(z_scores)))            #To test the transformed data
  
  colnames(long_data)[1] = "Value"
  long_data$ASVs = rownames(long_data)
  long_data$Metaloc <- i
  long_data$Density = point_densities
  long_data$log_rel = log_means             #To test the logarithm of the relative abundance
  long_data$rel = colMeans(data_meta)       #To test the relative abundance
  # long_data$Count_Meta = round(means)
  
  all_analysis <- rbind(all_analysis, long_data)
  long_data_counts <- as.data.frame(round(means))       #To test the rarefied counts
  
  colnames(long_data_counts)[1] = "Value"
  long_data_counts$ASVs = rownames(long_data_counts)
  long_data_counts$Metaloc <- i
  
  all_analysis_counts <- rbind(all_analysis_counts, long_data_counts)
  
}



vuong_list <- list()
gof_list <- list()
parameter_results <- data.frame()
complete_fits <- list()

indval <- read.csv("indval_M7.csv", header = T, row.names = 1)

for (i in unique(all_analysis$Metaloc))
{
  all_analysis_metaloc <- all_analysis %>% filter(Metaloc == i)
  
  values <- all_analysis_metaloc$Value
  
  fit_n <- fitdist(values, "norm")
  
  norm_mean <- fit_n$estimate["mean"]
  norm_sd   <- fit_n$estimate["sd"]
  
  # B. Shifted Weibull
  start_values_w <- list(
    shape = 1, 
    scale = mean(values) - min(values), 
    shift = min(values) - 0.1
  )
  fit_w <- fitdist(values, "weibull_shifted", start = start_values_w)
  
  weibull_shape <- fit_w$estimate["shape"]
  weibull_scale <- fit_w$estimate["scale"]
  weibull_shift <- fit_w$estimate["shift"]
  
  # C. Shifted Gamma (NEW: fitted with fitdist so that gofstat can be used)
  # The method-of-moments estimates are used as the starting point of the algorithm
  M <- mean(values, na.rm = TRUE)
  V <- var(values, na.rm = TRUE)
  S <- skewness(values, na.rm = TRUE)
  
  calc_shape <- 4 / (S^2)
  calc_rate  <- sqrt(calc_shape / V)
  calc_shift <- M - (calc_shape / calc_rate)
  
  start_values_g <- list(shape = calc_shape, rate = calc_rate, shift = calc_shift)
  fit_g <- fitdist(values, "gamma_shifted", start = start_values_g)
  
  all_analysis_metaloc$indval <- indval[all_analysis_metaloc$ASVs, "Location"]
  
  pdf(paste0( "MAD_", i, ".pdf"), width = 15, height = 12)
  # 1. Build the text block with the parameters (rounded to 3 decimals)
  parameters_label <- sprintf(
    "Fitted Parameters:\nNormal: Mean = %.3f, SD = %.3f\nShifted Gamma: Shape = %.3f, Rate = %.3f, Shift = %.3f\nShifted Weibull: Shape = %.3f, Scale = %.3f, Shift = %.3f",
    norm_mean, norm_sd,
    calc_shape, calc_rate, calc_shift, # Or use fit_g$estimate if the Gamma was fitted with fitdist
    weibull_shape, weibull_scale, weibull_shift
  )
  
  N <- nrow(all_analysis_metaloc)
  
  # 2. Build the plot
  p <- ggplot(all_analysis_metaloc, aes(x = Value)) +
    
    # Empirical density plot
    geom_density(fill = "grey", color = "gray30", alpha = 0.6) +
    
    # Add Theoretical Shifted Gamma
    stat_function(
      aes(color = "Shifted Gamma"), 
      fun = dgamma_shifted, 
      args = list(shape = calc_shape, rate = calc_rate, shift = calc_shift), 
      linewidth = 1.2, 
      linetype = "solid"
    ) +
    
    # Add Fitted Normal (Updated from Standard Normal)
    stat_function(
      aes(color = "Fitted Normal"), 
      fun = dnorm, 
      args = list(mean = norm_mean, sd = norm_sd), 
      linewidth = 1.2, 
      linetype = "solid"
    ) +
    
    # Add Shifted Weibull
    stat_function(
      aes(color = "Shifted Weibull"), 
      fun = dweibull_shifted, 
      args = list(shape = weibull_shape, scale = weibull_scale, shift = weibull_shift), 
      linewidth = 1.2, 
      linetype = "solid"
    ) +
    
    # Add parameter text box to the plot
    annotate(
      "text", 
      x = -Inf, y = Inf,          # Top-Left corner of the plot
      hjust = -0.05, vjust = 1.2, # Padding so it doesn't touch the axis lines
      label = parameters_label, 
      size = 3.5,                 # Text size
      family = "mono",            # Monospaced font makes numbers align neatly
      color = "black"
    ) +
    
    # Unified legend aesthetic control
    scale_color_manual(
      name = "Theoretical Models",
      values = c(
        "Shifted Gamma" = "blue",
        "Fitted Normal" = "red",
        "Shifted Weibull" = "green" 
      ) 
    ) +
    
    # ---> THIS IS THE CHANGE <---
    # Distributions per combination of TH and indval (NAs omitted)
    geom_density(data = all_analysis_metaloc %>% filter(!is.na(indval)),
                 aes(x = Value, 
                     y = after_stat(count) / N, 
                     fill = factor(indval)), 
                 alpha = 0.5, 
                 color = "black") +
    
    # Plot aesthetics
    theme_minimal() +
    theme(
      legend.position = c(0.8, 0.8), # Legend inside the plot (Top-Right)
      # Added a subtle border to the legend so it stands out from the grid lines
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3) 
    ) +
    labs(
      title = "Third law: Mean Abundance Distribution",
      subtitle = "Comparison of Theoretical Fits",
      x = "Transformed log relative abundance",
      y = "Probability Density",
      fill = "TH - Indval"   
    )
  
  # 3. Save and display
  plot(p)
  
  dev.off()
  
   
}
