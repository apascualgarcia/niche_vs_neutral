library(tidyverse)
library(nleqslv)
library(fitdistrplus)
library(fitteR)

library(goftest)

#### HERE I PLOT THE NUMBER OF POSITIVE OR NEGATIVE ASSOCIATIONS AGAINST THE FITS OBTAINED FOR GRILLI'S FIRST LAW######

get_r2_pp <- function(fit) {
  # 1. Get the empirical probabilities
  n <- fit$n
  p_empirical <- (1:n - 0.5) / n
  
  # 2. Compute the theoretical probabilities from the fit
  # fit$distname gives the name of the distribution (norm, gamma, etc.)
  dist_func <- paste0("p", fit$distname)
  p_theoretical <- do.call(dist_func, c(list(sort(fit$data)), as.list(fit$estimate)))
  
  # 3. Compute R2 against the identity line (y = x)
  rss <- sum((p_theoretical - p_empirical)^2) # Residual sum of squares
  tss <- sum((p_empirical - mean(p_empirical))^2) # Total sum of squares
  r_squared <- 1 - (rss / tss)
  
  return(r_squared)
}


data <- t(read.table(file = "../5.Real_data/Single_tree/otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = "../5.Real_data/Single_tree/Metadata_singleTree_wDirich-K10.tsv", header = T, row.names = 1, sep = "\t")
cor <- read.table("../6.Correlation_Distribution/Correlation_Meta_Com_025.tsv", header = T, sep = "\t")

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

all_analysis <- data.frame()
for (i in unique(meta$Dirichlet[!is.na(meta$Dirichlet)]))
{
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
        cols = -Samples,                       # 2. Transform all the columns EXCEPT the one just created
        names_to = "ASVs",
        values_to = "Value"
      )
 

  long_data_data <- as.data.frame(data_meta) %>% 
    rownames_to_column(var = "Samples") %>%  
    pivot_longer(
      cols = -Samples,                       # 2. Transform all the columns EXCEPT the one just created
      names_to = "ASVs",
      values_to = "Value"
    )
  
  long_data$Rel_abun = long_data_data$Value
  
  
  long_data$Metaloc <- i
  
  all_analysis <- rbind(all_analysis, long_data)
  
  }


for (i in unique(all_analysis$Metaloc))
{
  data_metaCom <- all_analysis %>% filter(Metaloc == i)
  
  z <- data_metaCom$Value
  shift <- abs(min(z)) + 0.001 
  z_shifted <- z + shift
  
  x_vals <- seq(min(z), max(z), length.out = length(z))
  # For gamma
  shape_val <- mean(z + abs(min(z)))^2 / var(z) 
  rate_val <- shape_val / mean(z + abs(min(z)))
  prob_density_gamma <- dgamma(x_vals + abs(min(z)), shape = shape_val, rate = rate_val)

  # For lognormal
  meanlog_val <- mean(log(z_shifted))
  sdlog_val <- sd(log(z_shifted))
  prob_density_LN <- dlnorm(x_vals + shift, meanlog = meanlog_val, sdlog = sdlog_val)
  
  
  real_density_approx <- density(z, from = min(z), to = max(z), n = length(z))$y
  residuals_gamma <- real_density_approx - rate_val
  residuals_ln    <- real_density_approx - prob_density_LN
  
  fit_n <- fitdist(z_shifted, "lnorm")
  fit_g <- fitdist(z_shifted, "gamma")
  
  a<- (ks.test(z_shifted, "plnorm", mean = meanlog_val, sd = sdlog_val))
  
  
  cor[cor$MetaCom ==paste0("Meta_", i),"Var_log"] = sdlog_val
  cor[cor$MetaCom ==paste0("Meta_", i),"Mean_log"] = meanlog_val
  cor[cor$MetaCom ==paste0("Meta_", i),"Shape_g"] = shape_val
  cor[cor$MetaCom ==paste0("Meta_", i),"Rate_g"] = rate_val
  cor[cor$MetaCom ==paste0("Meta_", i),"Dev_g"] = get_r2_pp(fit_g)
  cor[cor$MetaCom ==paste0("Meta_", i),"Dev_LN"] = get_r2_pp(fit_n)
  cor[cor$MetaCom ==paste0("Meta_", i),"KS_log"] = a$statistic
}

z <- all_analysis$Value
shift <- abs(min(z)) + 0.001 
z_shifted <- z + shift

x_vals <- seq(min(z), max(z), length.out = length(z))
# For gamma
shape_val <- mean(z + abs(min(z)))^2 / var(z) 
rate_val <- shape_val / mean(z + abs(min(z)))
prob_density_gamma <- dgamma(x_vals + abs(min(z)), shape = shape_val, rate = rate_val)
# For lognormal
meanlog_val <- mean(log(z_shifted))
sdlog_val <- sd(log(z_shifted))
prob_density_LN <- dlnorm(x_vals + shift, meanlog = meanlog_val, sdlog = sdlog_val)


real_density_approx <- density(z, from = min(z), to = max(z), n = length(z))$y
residuals_gamma <- real_density_approx - rate_val
residuals_ln    <- real_density_approx - prob_density_LN

fit_n <- fitdist(z_shifted, "lnorm")
fit_g <- fitdist(z_shifted, "gamma")

cat("R2 for Normal:", get_r2_pp(fit_n), "\n")
cat("R2 for Gamma:", get_r2_pp(fit_g), "\n")
par(mfrow = c(1, 2))
ppcomp(list(fit_n, fit_g), legendtext = c("Normal", "Gamma"), main = "P-P Plot")

a<- (ks.test(z_shifted, "plnorm", mean = meanlog_val, sd = sdlog_val))


cor[cor$MetaCom =="Total","Var_log"] = sdlog_val
cor[cor$MetaCom =="Total","Mean_log"] = meanlog_val
cor[cor$MetaCom =="Total","Shape_g"] = shape_val
cor[cor$MetaCom =="Total","Rate_g"] = rate_val
cor[cor$MetaCom ==paste0("Meta_", i),"Dev_g"] = get_r2_pp(fit_g)
cor[cor$MetaCom ==paste0("Meta_", i),"Dev_LN"] = get_r2_pp(fit_n)
cor[cor$MetaCom ==paste0("Meta_", i),"KS_log"] = a$statistic

cor_re <- cor %>% filter(MetaCom != "Total")

ggplot(cor_re, aes(x = KS_log, y =  Positive + Negative, color = MetaCom)) + 
  geom_point(size = 5) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = F) +
  scale_x_log10() +
  scale_y_log10()

ggplot(cor_re, aes(x = Var_log, y = Positive + Negative, color = MetaCom)) + 
  geom_point(size = 5) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = F) +
  scale_x_log10() +
  scale_y_log10()

ggplot(cor_re) + 
  geom_point(aes(x = as.numeric(Mean_log), y = Negative + Positive, color = MetaCom), size = 5) + 
  scale_x_log10() +
  scale_y_log10()  

ggplot(cor_re) + 
  geom_point(aes(x = as.numeric(Shape_g), y = Negative + Positive, color = MetaCom), size = 5) + 
  scale_x_log10() +
  scale_y_log10()  

ggplot(cor_re) + 
  geom_point(aes(x = as.numeric(Rate_g), y = Negative + Positive, color = MetaCom), size = 5) + 
  scale_x_log10() +
  scale_y_log10()  



### I will plot a score of the abundances against the variance
cors <- read.table("../6.Correlation_Distribution/Correlations.tsv", header = T, sep = "\t")

rownames(cor) = cor$MetaCom

for (i in c(1,2,3,4,5,6,7,8,10))
{
  cors_red <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i) & abs(value) > 0.25)
  cor[paste0("Meta_",i), "Sum"] = sum(abs(cors_red$value))
  cor[paste0("Meta_",i), "Mean"] = mean(abs(cors_red$value))
}

cor_re <- cor %>% filter(MetaCom != "Total")

ggplot(cor_re, aes(x = Var_log, y = Mean, color = MetaCom)) + 
  geom_point(size = 5) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = F) +
  scale_x_log10() +
  scale_y_log10()


ggplot(cor_re, aes(x = Positive + Negative, y = Mean, color = MetaCom)) + 
  geom_point(size = 5) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = F) +
  scale_x_log10() +
  scale_y_log10()

