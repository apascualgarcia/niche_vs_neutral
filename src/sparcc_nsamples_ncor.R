###############################################################################
# Number of SparCC links vs number of samples of each metacommunity
#
# Links with |z| > 2..10 and with |r| > 0.3 (nsamples_ncor_T03.pdf, used in the paper).
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlations.tsv
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output:
#   4.3_Niche_Interactions/figures/thresholds/nsamples_ncor*.pdf, Threshold_ncor_zscore.pdf
#
# Threshold_ncor_zscore.pdf was called Threshold_ncor.pdf, like the figure of
# src/sparcc_threshold_ncor.R; renamed so that it does not overwrite it.
# Figures are saved without legend (src/legend_utils.R). Their legends are shared with other
# scripts: figures/legends/legend_metacommunity.pdf (nsamples_ncor_T03.pdf) and
# figures/thresholds/legend_metacommunity_hue.pdf (the rest), so this script writes none.
# Run from anywhere inside the repository: Rscript src/sparcc_nsamples_ncor.R
###############################################################################
library(here)
source(here("src", "legend_utils.R"))  # save_without_legend(), save_metacommunity_legend()
fig.dir <- here("4.3_Niche_Interactions", "figures", "thresholds")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)

cors <- read.table(here("4.3_Niche_Interactions", "correlations", "Correlations.tsv"), header = T, sep = "\t")
meta <- read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, sep = "\t", row.names = 1)
data <- t(read.table(file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"), header = T, row.names = 1, sep = "\t"))

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

meta_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                  "#6D3F5B", "#1E213D", "#763C28")

summary_df <- data.frame()
for (i in c(1,2,3,4,5,6,7,8,10))
{
  data_meta_SparCC_filter <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i))
  data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
  
  met_red <- meta %>% filter(Dirichlet == i)
  
  
  for (j in 2:10)
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
  
}
}
## Associations strength
ggplot(summary_df, aes(x = nsamples, y = Mean, color = factor(Metacom), group = factor(Threshold))) + 
  geom_point(size = 6) + 
  geom_line() +
  geom_errorbar(aes(x=nsamples, ymin= Mean-se, ymax=Mean+se), width=0.4)

## nAssociations
p <- ggplot(summary_df, aes(x = nsamples, y = n_cor, color = factor(Metacom), group = factor(Threshold))) + 
  geom_point(size = 6) +
  geom_line()
save_without_legend(p, file.path(fig.dir, "nsamples_ncor.pdf"), width = 15, height = 10,
                    legend_file = NULL)  # legend_metacommunity_hue.pdf (src/sparcc_threshold_ncor.R)


##Check roughly where to set the threshold
data_meta_SparCC_filter_met7 <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",7))

a <- hist(data_meta_SparCC_filter_met7$value, breaks = 1000)
view(data.frame(counts = a$counts, minds = a$mids))
## At about 0.3


summary_cor <- data.frame()
for (i in c(1,2,3,4,5,6,7,8,10))
{
  data_meta_SparCC_filter <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i))
  data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
  
  met_red <- meta %>% filter(Dirichlet == i)
  
  cors_red <- data_meta_SparCC_filter %>% filter(abs(value) > 0.3)
summary_cor <- rbind(summary_cor, data.frame(Metacom = i, 
                                 Mean = mean(abs(cors_red$value)), 
                                 nsamples = nrow(met_red), 
                                 n_cor = nrow(cors_red),
                                 se = sd(abs(cors_red$value)), 
                                 Mean_taxon = mean(rowSums(data_rel[rownames(met_red),]>0)),
                                 sd_taxon = sd(rowSums(data_rel[rownames(met_red),]>0)),
                                 Mean_possible_associations = mean(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                 sd_possible_associations = sd(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                 Threshold = 0.3
) 
)
}

## By correlation
res <- summary_cor %>% filter(Threshold == 0.3)
p <- ggplot(res, aes(x = nsamples, y = n_cor, color = factor(Metacom), group = factor(Threshold))) + 
  geom_point(size = 10) +  
  labs(x = "Number of samples",
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
  guides(fill = guide_legend(override.aes = list(size = 8)))
save_without_legend(p, file.path(fig.dir, "nsamples_ncor_T03.pdf"), width = 15, height = 10,
                    legend_file = NULL)  # legends/legend_metacommunity.pdf (src/sparcc_correlation_distribution.R)



## By z-score
res <- summary_df %>% filter(Threshold == 3)
p <- ggplot(res, aes(x = nsamples, y = n_cor, color = factor(Metacom), group = factor(Threshold))) + 
  geom_point(size = 6) 
save_without_legend(p, file.path(fig.dir, "nsamples_ncor_T3.pdf"), width = 15, height = 10,
                    legend_file = NULL)  # legend_metacommunity_hue.pdf (src/sparcc_threshold_ncor.R)



## NTaxons
ggplot(summary_df, aes(x = nsamples, y = Mean_taxon, color = factor(Metacom))) + 
  geom_point(size = 6) + 
  geom_errorbar(aes(x=nsamples, ymin= Mean_taxon-sd_taxon, ymax=Mean_taxon+sd_taxon), width=0.4)
##  POssible associations
ggplot(summary_df, aes(x = nsamples, y = Mean_possible_associations, color = factor(Metacom))) + 
  geom_point(size = 6) + 
  geom_errorbar(aes(x=nsamples, ymin= Mean_possible_associations-sd_possible_associations, ymax=Mean_possible_associations+sd_possible_associations), width=0.4)

## Real correlations vs possible associations
ggplot(summary_df, aes(x = n_cor, y = Mean_possible_associations, color = factor(Metacom))) + 
  geom_point(size = 6) + 
  geom_errorbar(aes(x=n_cor, ymin= Mean_possible_associations-sd_possible_associations, ymax=Mean_possible_associations+sd_possible_associations), width=0.4)
## Threshold
p <- ggplot(summary_df, aes(x = Threshold, y = n_cor, color = factor(Metacom))) + 
  geom_point(size = 6) +
  geom_line() +
  scale_y_log10()
save_without_legend(p, file.path(fig.dir, "Threshold_ncor_zscore.pdf"), width = 15, height = 10,
                    legend_file = NULL)  # legend_metacommunity_hue.pdf (src/sparcc_threshold_ncor.R)

##Histogram of associations of 2
hist((cors %>% filter(method == "Pearson" & MetaCom == paste0("Meta_",2)))$value)




