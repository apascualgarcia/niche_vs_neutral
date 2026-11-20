###############################################################################
# Number of SparCC links as a function of the threshold
#
# First by z-score, then by the SparCC score itself. The second block overwrites
# the PDFs of the first: Threshold_ncor.pdf (score threshold, red line at 0.3) is
# the one used in the paper.
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlations.tsv
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output:
#   4.3_Niche_Interactions/figures/thresholds/{Threshold_ncor,ncor_Threshold_}.pdf
#   4.3_Niche_Interactions/figures/thresholds/legend_metacommunity_hue.pdf (also used by src/sparcc_nsamples_ncor.R)
# Figures are saved without legend (src/legend_utils.R); the legend of Threshold_ncor.pdf is
# figures/legends/legend_metacommunity.pdf, written by src/sparcc_correlation_distribution.R.
# Run from anywhere inside the repository: Rscript src/sparcc_threshold_ncor.R
###############################################################################
library(here)
source(here("src", "legend_utils.R"))  # save_without_legend(), save_metacommunity_legend()
fig.dir <- here("4.3_Niche_Interactions", "figures", "thresholds")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)

### HERE I plot the number of correlations against the threshold, to cut at the right number.

cors <- read.table(here("4.3_Niche_Interactions", "correlations", "Correlations.tsv"), header = T, sep = "\t")
meta <- read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, sep = "\t", row.names = 1)
data <- t(read.table(file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"), header = T, row.names = 1, sep = "\t"))

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

meta_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                  "#6D3F5B", "#1E213D", "#763C28")


#### By zscore
summary_df <- data.frame()

ncor_thres <- data.frame()

for (i in c(1,2,3,4,5,6,7,8,10))
{
  data_meta_SparCC_filter <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i))
  data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
  
  met_red <- meta %>% filter(Dirichlet == i)
  
  inter <- seq(2,max(abs(data_meta_SparCC_filter$zscore)), length = 1000)
  pass = TRUE
  for (j in inter)
{  
    
    cors_red <- data_meta_SparCC_filter %>% filter(abs(zscore) > j)
  summary_df <- rbind(summary_df, data.frame(Metacom = i, 
                                   Mean = mean(abs(cors_red$value)), 
                                   nsamples = nrow(met_red), 
                                   n_cor = nrow(cors_red),
                                   se = sd(abs(cors_red$value)), 
                                   Mean_taxon = mean(rowSums(data_rel[rownames(met_red),]>0)),
                                   sd_taxon = sd(rowSums(data_rel[rownames(met_red),]>0)),
                                   Mean_possible_associations = mean(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                   sd_possible_associations = sd(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                   Threshold = j
                 ) 
                 )
  
  
  if (pass & nrow(cors_red) <= 1000)
  {
  ncor_thres <- rbind(ncor_thres, data.frame(Metacom = i,
                                                  Threshold = j,
                                                  n_cor = nrow(cors_red)
                                                  
  ))
  pass = FALSE
  }
}
}




## Threshold
p <- ggplot(summary_df, aes(x = Threshold, y = n_cor, color = factor(Metacom))) + 
  geom_point() +
  geom_line() +
  geom_hline(yintercept=1000) +
  scale_y_log10()
save_without_legend(p, file.path(fig.dir, "Threshold_ncor.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_metacommunity_hue.pdf"))


p <- ggplot(summary_df, aes(x = n_cor, y = Threshold, color = factor(Metacom))) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept=1000) +
  scale_x_log10()
save_without_legend(p, file.path(fig.dir, "ncor_threshold.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_metacommunity_hue.pdf"))




## By correlation
summary_df <- data.frame()

ncor_thres <- data.frame()

for (i in c(1,2,3,4,5,6,7,8,10))
{
  data_meta_SparCC_filter <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i))
 
  met_red <- meta %>% filter(Dirichlet == i)
  
  inter <- seq(0,max(abs(data_meta_SparCC_filter$value)), length = 1000)
  pass = TRUE
  for (j in inter)
  {  
    
    cors_red <- data_meta_SparCC_filter %>% filter(abs(value) > j)
    summary_df <- rbind(summary_df, data.frame(Metacom = i, 
                                     Mean = mean(abs(cors_red$value)), 
                                     nsamples = nrow(met_red), 
                                     n_cor = nrow(cors_red),
                                     se = sd(abs(cors_red$value)), 
                                     Mean_taxon = mean(rowSums(data_rel[rownames(met_red),]>0)),
                                     sd_taxon = sd(rowSums(data_rel[rownames(met_red),]>0)),
                                     Mean_possible_associations = mean(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                     sd_possible_associations = sd(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                     Threshold = j
    ) 
    )
  }
}




## Threshold
p <- ggplot(summary_df, aes(x = Threshold, y = n_cor/2, color = factor(Metacom))) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept=0.3, color = "red", linetype = "dashed", linewidth = 2) +
  scale_y_log10() +
  labs(x = "Threshold Correlation Weight of SparCC",
       y = "Number of SparCC Correlations",
       color = "Metacommunity") + 
  scale_color_manual(
    values = meta_palette
  ) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),  
    axis.text.y = element_text(size = 20),  
    legend.title = element_text(size = 25), 
    legend.text = element_text(size = 20),  
    legend.key.size = unit(1.5, "lines"),
    legend.position = "none"
  ) +
  guides(color = guide_legend(override.aes = list(size = 8)))
save_without_legend(p, file.path(fig.dir, "Threshold_ncor.pdf"), width = 15, height = 10,
                    legend_file = NULL)  # legends/legend_metacommunity.pdf (src/sparcc_correlation_distribution.R)


p <- ggplot(summary_df, aes(x = n_cor, y = Threshold, color = factor(Metacom))) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept=1000) +
  scale_x_log10()
save_without_legend(p, file.path(fig.dir, "ncor_threshold.pdf"), width = 15, height = 10,
                    legend_file = file.path(fig.dir, "legend_metacommunity_hue.pdf"))

