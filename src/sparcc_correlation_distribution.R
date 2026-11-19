###############################################################################
# Distribution of the SparCC scores per metacommunity and per tree hole
#
# Gathers the SparCC edge lists (and the Pearson correlations of the relative
# abundances, for comparison), computes their z-score within each metacommunity,
# plots their density and cumulative distribution and counts the links above
# |r| = 0.25, 0.5 and |z| = 2..10. Repeated for the tree-hole networks.
#
# Input:
#   4.3_Niche_Interactions/SparCC/{Meta_*,Hole_*}/{correlations.edgelist,data.tsv}
# Output:
#   4.3_Niche_Interactions/correlations/*.tsv (Correlations.tsv and the *zscore_thres* tables are large and not committed)
#   4.3_Niche_Interactions/figures/correlations/*.pdf (paper: correlation_dist.pdf = Fig. 4B, cumulative_distribution_cor.pdf)
#   4.3_Niche_Interactions/figures/legends/legend_metacommunity.pdf (legend of Fig. 4B and of every
#     figure coloured with meta_palette, also those of src/sparcc_threshold_ncor.R and src/sparcc_nsamples_ncor.R)
#
# Figures are saved without legend; legends go to separate PDFs (src/legend_utils.R):
# legend_metacommunity_paired.pdf, legend_treehole_paired.pdf, cumulative_distribution_holes_legend.pdf.
# Run from anywhere inside the repository: Rscript src/sparcc_correlation_distribution.R
###############################################################################
library(here)
source(here("src", "legend_utils.R"))  # save_without_legend(), save_metacommunity_legend()
sparcc.dir <- here("4.3_Niche_Interactions", "SparCC")
cor.dir <- here("4.3_Niche_Interactions", "correlations")
fig.dir <- here("4.3_Niche_Interactions", "figures", "correlations")
dir.create(cor.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)

data_meta <- data.frame()
for (i in c(paste0("Meta_" ,c(1,2,3,4,5,6,7,8, 10))))
{
  Sparcc_cor <- read.table(file.path(sparcc.dir, i, "correlations.edgelist"), header = T, sep = "\t")
  Abd <- read.table(file.path(sparcc.dir, i, "data.tsv"), header = T, row.names =  1, sep = "\t")
  Abd <- Abd/rowSums(Abd)
  pearson_cor <- cor(t(Abd))
  idx <- which(upper.tri(pearson_cor), arr.ind = TRUE)
  row_names <- rownames(pearson_cor)[idx[, "row"]]
  col_names  <- colnames(pearson_cor)[idx[, "col"]]
  pearson_values <- pearson_cor[idx]
  
  
  pearson_cor <- pearson_cor[upper.tri(pearson_cor)]
  
  hist(pearson_cor)
  hist(Sparcc_cor$cor)
  
  
  df_plot <- data.frame(
    value = c(pearson_cor, Sparcc_cor$cor),
    method = c(rep("Pearson", length(pearson_cor)), 
               rep("SparCC", length(Sparcc_cor$cor))),
    MetaCom = i,
    taxa_1 = c(row_names, Sparcc_cor$Sp1), 
    taxa_2 = c(col_names, Sparcc_cor$Sp2)
  )
  
  data_meta <- rbind(data_meta, df_plot)
}

data_meta <- data_meta %>% filter(!is.na(value))
write.table(data_meta, file = file.path(cor.dir, "Correlations.tsv"), sep = "\t", col.names = T, row.names = F)

meta_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                  "#6D3F5B", "#1E213D", "#763C28")

# Legend shared by every figure coloured with meta_palette (Fig. 4B and the threshold figures);
# those figures are saved without legend (legend_file = NULL below)
leg.dir <- here("4.3_Niche_Interactions", "figures", "legends")
dir.create(leg.dir, recursive = TRUE, showWarnings = FALSE)
save_metacommunity_legend(file.path(leg.dir, "legend_metacommunity.pdf"), palette = meta_palette)

########### Correlation coefficient vs number of correlations
library(parallel)

data_meta <- read.table(file.path(cor.dir, "Correlations.tsv"), header = T, sep = "\t")

data_meta_SparCC <- data_meta %>% filter(method == "SparCC")
nbins <- floor(nrow(data_meta_SparCC) / 100)
inter <- max(abs(data_meta_SparCC$value)) / nbins
cut_points <- (0:nbins) * inter
abs_values <- abs(data_meta_SparCC$value) # Extract the vector so it is lighter to pass to the cores
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c("abs_values"))
n_counts <- parSapply(cl, cut_points, function(x) {
  sum(abs_values >= x)
})

stopCluster(cl)

summary_df <- data.frame(
  min = cut_points,
  nCors = n_counts
)

abs_values <- abs(data_meta_SparCC$value)
ecdf_fun <- ecdf(abs_values)

cut_points <- seq(0, max(abs_values), length.out = 100)

summary_df <- data.frame(
  threshold = cut_points,
  prob = 1 - ecdf_fun(cut_points) 
)
ggplot(summary_df) +
  geom_point( aes(x = threshold, y = prob)) +  # "+" was missing (exploratory plot, not saved)
  scale_y_log10()  +
  scale_x_log10()


  

#### Per metacommunity #####
  summary_all = data.frame()
  data_meta_sparCC_zscore <- data.frame()
  for (i in paste0("Meta_" ,c(1,2,3,4,5,6,7,8, 10)))
  {
    data_meta_SparCC_filter <- data_meta_SparCC %>% filter(MetaCom == i)
    
    data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
    
    nbins <- (max(abs(data_meta_SparCC_filter$zscore)) - min(abs(data_meta_SparCC_filter$zscore)))/min(abs(data_meta_SparCC_filter$zscore))
    inter <- max(abs(data_meta_SparCC_filter$zscore)) / nbins
    cut_points <- (0:nbins) * inter
    abs_values <- abs(data_meta_SparCC_filter$zscore) # Extract the vector so it is lighter to pass to the cores
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("abs_values"))
    n_counts <- parSapply(cl, cut_points, function(x) {
      sum(abs_values >= x)
    })
    summary_tmp <- data.frame(
      min = cut_points,
      nCors = n_counts,
      MetaCom = i
    )
    summary_tmp$nCors = 1 - summary_tmp$nCors/summary_tmp$nCors[1]  # Normalise and make it cumulative
    
    summary_all = rbind(summary_all, summary_tmp)
    data_meta_sparCC_zscore = rbind(data_meta_sparCC_zscore, data_meta_SparCC_filter)
    print(i)
  }
  
  
  
  summary_all$MetaCom <- factor(summary_all$MetaCom, 
                               levels = str_sort(unique(summary_all$MetaCom), numeric = TRUE))
  
  
  
  
  p <- ggplot(summary_all) +
    geom_point(aes(x = min, y = nCors, color = MetaCom)) +
    scale_color_manual(values = meta_palette) +  # Colour palette
    labs(
      x = "Threshold ",
      y = "Normalized proportion of correlations",
      color = "Metacommunity"
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 14), 
      axis.title.y = element_text(size = 14),
      legend.title = element_text(size = 15), # Larger "Metacommunity" title
      legend.text = element_text(size = 13),  # Larger names (Meta_1, Meta_2...)
      legend.key.size = unit(1.5, "lines")    # Larger colour dots in the legend
    )
  save_without_legend(p, file.path(fig.dir, "cumulative_distribution.pdf"), width = 15, height = 10,
                      legend_file = NULL)  # legends/legend_metacommunity.pdf

  
  
  
  
  
  
  
  
  
  
  
  
  library(RColorBrewer)
  
  p <- ggplot(data_meta_sparCC_zscore, aes(x = zscore, fill = MetaCom)) +
    geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
    theme_minimal() +            # clean theme
    scale_fill_brewer(palette = "Paired") +
    labs(x = "Correlation index",
         y = "Probability Density",
         fill = "Meta Community")     # Legend label
  save_without_legend(p, file.path(fig.dir, "zscore_of_Correlations.pdf"), width = 15, height = 10,
                      legend_file = file.path(fig.dir, "legend_metacommunity_paired.pdf"))
  
  

  ##### Reconstructed (Nov 2026): the code that wrote these three files was lost from this
  ##### script. It mirrors the tree-hole block at the end of the script (|z-score| >= 5).
  ##### Correlation_zscore_thres.tsv holds ALL the SparCC scores with their z-score (not only
  ##### those above the threshold); it is the input of src/networks_M7.R.
  data_meta_SparCC_filtered <- data_meta_sparCC_zscore %>% filter(abs(zscore) >= 5 )
  pos <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value > 0))$MetaCom)
  neg <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value < 0))$MetaCom)
  table <- data.frame(MetaCom = names(pos), Positive = as.vector(pos), Negative = as.vector(neg))
  write.table(table, file = file.path(cor.dir, "Correlation_Meta_Com_zscore_thres.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)
  p <- ggplot(data_meta_SparCC_filtered, aes(x = value, fill = MetaCom)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    scale_fill_brewer(palette = "Paired") +
    labs(x = "Correlation index",
         y = "Probability Density",
         fill = "Meta Community")
  save_without_legend(p, file.path(fig.dir, "Correlations_zscore_thres.pdf"), width = 15, height = 10,
                      legend_file = file.path(fig.dir, "legend_metacommunity_paired.pdf"))
  write.table(data_meta_sparCC_zscore, file = file.path(cor.dir, "Correlation_zscore_thres.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)
  # Kept for the |z| >= 2..10 counts below: the next loop rebuilds data_meta_sparCC_zscore without
  # the zscore column (the original counts were made on these data; MetaCom written as Meta_1).
  data_meta_sparCC_zscore_z <- data_meta_sparCC_zscore

  ##### Instead of by z-score, I will do it by correlation value
  
  summary_all = data.frame()
  data_meta_sparCC_zscore <- data.frame()
  for (i in paste0("Meta_" ,c(1,2,3,4,5,6,7,8, 10)))
  {
    data_meta_SparCC_filter <- data_meta_SparCC %>% filter(MetaCom == i)
    
    nbins <- (max(abs(data_meta_SparCC_filter$value)) - min(abs(data_meta_SparCC_filter$value)))/min(abs(data_meta_SparCC_filter$value))
    inter <- max(abs(data_meta_SparCC_filter$value)) / nbins
    cut_points <- (0:nbins) * inter
    abs_values <- abs(data_meta_SparCC_filter$value) # Extract the vector so it is lighter to pass to the cores
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("abs_values"))
    n_counts <- parSapply(cl, cut_points, function(x) {
      sum(abs_values >= x)
    })
    summary_tmp <- data.frame(
      min = cut_points,
      nCors = n_counts,
      MetaCom = i
    )
    summary_tmp$nCors = 1 - summary_tmp$nCors/summary_tmp$nCors[1]  # Normalise and make it cumulative
    
    summary_all = rbind(summary_all, summary_tmp)
    data_meta_sparCC_zscore = rbind(data_meta_sparCC_zscore, data_meta_SparCC_filter)
    print(i)
  }
  
  
  summary_all$MetaCom <- str_replace(summary_all$MetaCom, "_", " ")
  summary_all$MetaCom <- factor(summary_all$MetaCom, 
                               levels = str_sort(unique(summary_all$MetaCom), numeric = TRUE))
  
  
  
  
  p <- ggplot(summary_all) +
    geom_point(aes(x = min, y = nCors, color = MetaCom)) +
    scale_color_manual(values = meta_palette) +  # Palette
    labs(
      x = "Correlation Threshold",
      y = "Cumulative correlation distribution",
      color = "Metacommunity"
    ) +
    geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", linewidth = 2) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 25), 
      axis.title.y = element_text(size = 25),
      axis.text.x = element_text(size = 20),  # Larger X-axis numbers
      axis.text.y = element_text(size = 20),  # And larger Y-axis numbers
      legend.title = element_text(size = 25), 
      legend.text = element_text(size = 20),  
      legend.key.size = unit(1.5, "lines"),
      legend.position = "none"
    ) +
    # This enlarges the points in the legend
    guides(color = guide_legend(override.aes = list(size = 8)))
  save_without_legend(p, file.path(fig.dir, "cumulative_distribution_cor.pdf"), width = 15, height = 10,
                      legend_file = NULL)  # legends/legend_metacommunity.pdf
  
  
  
  data_meta_sparCC_zscore$MetaCom <- str_replace(data_meta_sparCC_zscore$MetaCom, "_", " ")
  data_meta_sparCC_zscore$MetaCom <- factor(data_meta_sparCC_zscore$MetaCom, 
                               levels = str_sort(unique(data_meta_sparCC_zscore$MetaCom), numeric = TRUE))
  
  data_meta_sparCC_zscore <- data_meta_sparCC_zscore %>%
    mutate(MetaCom = fct_reorder(MetaCom, value, .fun = median))
  
  p <- ggplot(data_meta_sparCC_zscore, aes(x = value, fill = MetaCom)) +
    geom_density(alpha = 0.7) +  
    geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", linewidth = 2) +
    geom_vline(xintercept = -0.3, linetype = "dashed", color = "red", linewidth = 2) +
    # The trick: force the reading order, ignoring the internal order
    scale_fill_manual(
      values = meta_palette,
      breaks = str_sort(unique(data_meta_sparCC_zscore$MetaCom), numeric = TRUE)
    ) +  
    labs(x = "SparCC Correlation",
         y = "Probability Density",
         fill = "Metacommunity") + 
    theme_bw() + 
    theme(
      axis.title.x = element_text(size = 30), 
      axis.title.y = element_text(size = 30),
      axis.text.x = element_text(size = 25),  
      axis.text.y = element_text(size = 25),  
      legend.title = element_text(size = 25), 
      legend.text = element_text(size = 20),  
      legend.key.size = unit(1.5, "lines"),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 8)))
  save_without_legend(p, file.path(fig.dir, "correlation_dist.pdf"), width = 15, height = 10,
                      legend_file = NULL)  # legends/legend_metacommunity.pdf
  
  
  
  
  positive <- data.frame()
  negative <- data.frame()
  for (i in 2:10)
  {
    data_meta_SparCC_filtered <- data_meta_sparCC_zscore_z %>% filter(abs(zscore) >= i )  # was data_meta_sparCC_zscore, which has no zscore here
    
    pos <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value > 0))$MetaCom)
    
    positive[names(pos), paste0(">",i)] = pos
    
    neg <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value < 0))$MetaCom)
    
    negative[names(neg), paste0(">",i)] = neg
    
    pr <- data_meta_SparCC_filtered %>% filter(method == "SparCC" & abs(value) > 0)
    pr2 <- data_meta_SparCC_filtered %>% filter(method == "SparCC" & abs(value) < 0)
    
    library(RColorBrewer)
    
    p <- ggplot(pr, aes(x = value, fill = MetaCom)) +
      geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
      theme_minimal() +            # clean theme
      scale_fill_brewer(palette = "Paired") +
      labs(x = "Correlation index",
           y = "Probability Density",
           fill = "Meta Community")     # Legend label
    save_without_legend(p, file.path(fig.dir, paste0("Correlations_zscore_thres_", i, ".pdf")), width = 15, height = 10,
                        legend_file = file.path(fig.dir, "legend_metacommunity_paired.pdf"))
    
  }
  positive[is.na(positive)] = 0
  negative[is.na(negative)] = 0
  
  write.table(positive, file = file.path(cor.dir, "Positive_correlations.tsv"), quote = F, col.names = T, sep = "\t", row.names = T)
  write.table(negative, file = file.path(cor.dir, "Negative_correlations.tsv"), quote = F, col.names = T, sep = "\t", row.names = T)
  
  
#### Table and plot of the SparCC correlations for the single tree:#######
data_meta$MetaCom <- as.factor(data_meta$MetaCom)
# > 0.25



pos <- table((data_meta %>% filter(method == "SparCC" & value > 0.25))$MetaCom)
neg <- table((data_meta %>% filter(method == "SparCC" & value < -0.25))$MetaCom)

table <- data.frame(MetaCom = names(pos), Positive = as.vector(pos), Negative = as.vector(neg))

write.table(table, file = file.path(cor.dir, "Correlation_Meta_Com_025.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)



pr <- data_meta %>% filter(method == "SparCC" & abs(value) > 0.25)
pr2 <- data_meta %>% filter(method == "SparCC" & abs(value) < 0.25)

library(RColorBrewer)

p <- ggplot(pr, aes(x = value, fill = MetaCom)) +
  geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
  theme_minimal() +            # clean theme
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Correlation index",
       y = "Probability Density",
       fill = "Meta Community")     # Legend label
save_without_legend(p, file.path(fig.dir, "Correlations_025.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_metacommunity_paired.pdf"))

# > 0.5

pos <- table((data_meta %>% filter(method == "SparCC" & value > 0.5  ))$MetaCom)
neg <- table((data_meta %>% filter(method == "SparCC" & value < -0.5  ))$MetaCom)

table <- data.frame(MetaCom = names(pos), Positive = as.vector(pos), Negative = as.vector(neg))

write.table(table, file = file.path(cor.dir, "Correlation_Meta_Com_05.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)



pr <- data_meta %>% filter(method == "SparCC" & abs(value) > 0.5)
pr2 <- data_meta %>% filter(method == "SparCC" & abs(value) < 0.5)

library(RColorBrewer)

p <- ggplot(pr, aes(x = value, fill = MetaCom)) +
  geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
  theme_minimal() +            # clean theme
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Correlation index",
       y = "Probability Density",
       fill = "Meta Community")     # Legend label
save_without_legend(p, file.path(fig.dir, "Correlations_05.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_metacommunity_paired.pdf"))



####################FOR THE Holes################


library(tidyverse)
library(parallel)

data_meta <- data.frame()
for (i in c(paste0("Hole_" ,c(1,2,3,4,5,6,7,8, 9, 10))))
{
  Sparcc_cor <- read.table(file.path(sparcc.dir, i, "correlations.edgelist"), header = T, sep = "\t")
  Abd <- read.table(file.path(sparcc.dir, i, "data.tsv"), header = T, row.names =  1, sep = "\t")
  Abd <- Abd/rowSums(Abd)
  pearson_cor <- cor(t(Abd))
  idx <- which(upper.tri(pearson_cor), arr.ind = TRUE)
  row_names <- rownames(pearson_cor)[idx[, "row"]]
  col_names  <- colnames(pearson_cor)[idx[, "col"]]
  pearson_values <- pearson_cor[idx]
  
  
  pearson_cor <- pearson_cor[upper.tri(pearson_cor)]
  
  hist(pearson_cor)
  hist(Sparcc_cor$cor)
  
  
  df_plot <- data.frame(
    value = c(pearson_cor, Sparcc_cor$cor),
    method = c(rep("Pearson", length(pearson_cor)), 
               rep("SparCC", length(Sparcc_cor$cor))),
    MetaCom = i,
    taxa_1 = c(row_names, Sparcc_cor$Sp1), 
    taxa_2 = c(col_names, Sparcc_cor$Sp2)
  )
  
  data_meta <- rbind(data_meta, df_plot)
}

data_meta <- data_meta %>% filter(!is.na(value))
write.table(data_meta, file = file.path(cor.dir, "Correlations_holes.tsv"), sep = "\t", col.names = T, row.names = F)


data_meta_SparCC <- data_meta %>% filter(method == "SparCC")
nbins <- floor(nrow(data_meta_SparCC) / 100)
inter <- max(abs(data_meta_SparCC$value)) / nbins
cut_points <- (0:nbins) * inter
abs_values <- abs(data_meta_SparCC$value) # Extract the vector so it is lighter to pass to the cores
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, varlist = c("abs_values"))
n_counts <- parSapply(cl, cut_points, function(x) {
  sum(abs_values >= x)
})

stopCluster(cl)

summary_df <- data.frame(
  min = cut_points,
  nCors = n_counts
)

abs_values <- abs(data_meta_SparCC$value)
ecdf_fun <- ecdf(abs_values)

cut_points <- seq(0, max(abs_values), length.out = 100)

summary_df <- data.frame(
  threshold = cut_points,
  prob = 1 - ecdf_fun(cut_points) 
)
ggplot(summary_df) +
  geom_point( aes(x = threshold, y = prob)) +  # "+" was missing (exploratory plot, not saved)
scale_y_log10()  +
  scale_x_log10()




#### Per metacommunity #####
summary_all = data.frame()
data_meta_sparCC_zscore <- data.frame()
for (i in c(paste0("Hole_" ,c(1,2,3,4,5,6,7,8, 9, 10))))
{
  data_meta_SparCC_filter <- data_meta_SparCC %>% filter(MetaCom == i)
  
  data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
  
  nbins <- (max(abs(data_meta_SparCC_filter$zscore)) - min(abs(data_meta_SparCC_filter$zscore)))/min(abs(data_meta_SparCC_filter$zscore))
  inter <- max(abs(data_meta_SparCC_filter$zscore)) / nbins
  cut_points <- (0:nbins) * inter
  abs_values <- abs(data_meta_SparCC_filter$zscore) # Extract the vector so it is lighter to pass to the cores
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  clusterExport(cl, varlist = c("abs_values"))
  n_counts <- parSapply(cl, cut_points, function(x) {
    sum(abs_values >= x)
  })
  summary_tmp <- data.frame(
    min = cut_points,
    nCors = n_counts,
    MetaCom = i
  )
  summary_tmp$nCors = 1 - summary_tmp$nCors/summary_tmp$nCors[1]  # Normalise and make it cumulative
  
  summary_all = rbind(summary_all, summary_tmp)
  data_meta_sparCC_zscore = rbind(data_meta_sparCC_zscore, data_meta_SparCC_filter)
  print(i)
}

p <- ggplot(summary_all) +
  geom_point( aes(x = min, y = nCors, color = MetaCom))  
save_without_legend(p, file.path(fig.dir, "cumulative_distribution_holes.pdf"), width = 15, height = 10)

library(RColorBrewer)

p <- ggplot(data_meta_sparCC_zscore, aes(x = zscore, fill = MetaCom)) +
  geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
  theme_minimal() +            # clean theme
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Correlation index",
       y = "Probability Density",
       fill = "Meta Community")     # Legend label
save_without_legend(p, file.path(fig.dir, "zscore_of_Correlations_holes.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_treehole_paired.pdf"))



data_meta_SparCC_filtered <- data_meta_sparCC_zscore %>% filter(abs(zscore) >= 5 )

pos <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value > 0))$MetaCom)
neg <- table((data_meta_SparCC_filtered %>% filter(method == "SparCC" & value < 0))$MetaCom)

table <- data.frame(MetaCom = names(pos), Positive = as.vector(pos), Negative = as.vector(neg))

write.table(table, file = file.path(cor.dir, "Correlation_Meta_Com_zscore_thres_holes.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)



pr <- data_meta_SparCC_filtered %>% filter(method == "SparCC" & abs(value) > 0)
pr2 <- data_meta_SparCC_filtered %>% filter(method == "SparCC" & abs(value) < 0)

library(RColorBrewer)

p <- ggplot(pr, aes(x = value, fill = MetaCom)) +
  geom_density(alpha = 0.5) +  # alpha sets the transparency to show overlaps
  theme_minimal() +            # clean theme
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Correlation index",
       y = "Probability Density",
       fill = "Meta Community")     # Legend label
save_without_legend(p, file.path(fig.dir, "Correlations_zscore_thres_holes.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_treehole_paired.pdf"))

write.table(data_meta_sparCC_zscore, file = file.path(cor.dir, "Correlation_zscore_thres_holes.tsv"), quote = F, col.names = T, sep = "\t", row.names = F)
