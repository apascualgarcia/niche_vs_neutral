###############################################################################
# Tree-hole characteristic ASVs of M7 by abundance ranking
#
# Ranks the ASVs by summed relative abundance in TH1 and TH7 and selects as
# characteristic those whose rank differs by > 1000 or whose abundance differs by
# > 0.05 (ranking_TH column), next to the IndVal result (indval_TH column).
# The original script continued with exploratory sections for M10, M3 and M2 and a
# copy of src/indval_treeholes_M7.R; they are in
# src/legacy/interactions/ranking_treeholes_all_metacommunities.R.
#
# Input:
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
#   4.3_Niche_Interactions/indicator_species/indval_M7.csv
# Output:
#   4.3_Niche_Interactions/indicator_species/characteristic_ASVs_M7.tsv
#   4.3_Niche_Interactions/figures/indicator_species/*.pdf (DifRankingvsDifAbun_thres.pdf in the paper)
#   Figures are saved without legend; legends: Ranking_TH_legend.pdf (Ranking_TH1/TH7),
#   Ranking_Ind_legend.pdf (RankingvRanking, DifRankingvsDifAbun), DifRankingvsDifAbun_thres_legend.pdf
# Run from anywhere inside the repository: Rscript src/ranking_treeholes_M7.R
###############################################################################
library(here)
source(here("src", "legend_utils.R"))  # save_without_legend(), save_metacommunity_legend()
ind.dir <- here("4.3_Niche_Interactions", "indicator_species")
fig.dir <- here("4.3_Niche_Interactions", "figures", "indicator_species")
dir.create(ind.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(indicspecies)
library(dplyr)
library(vegan)
library(gplots)
library(tidyverse)

data <- t(read.table(file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"), header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, row.names = 1, sep = "\t")

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

data <- data[rownames(meta)[rownames(meta) %in% rownames(data)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

##### MET 7 ########


indval_res <- read.csv(file.path(ind.dir, "indval_M7.csv"))
meta_M7 <- meta %>% filter(Dirichlet == 7)

data_M7 <- data_rel[rownames(meta_M7),]
data_M7 <- data_M7[,colSums(data_M7 > 0) == nrow(data_M7)]

## MY INDVAL

meta_TH1 <- meta_M7 %>% filter(Location == 1)
meta_TH7 <- meta_M7 %>% filter(Location == 7)

data_TH1 <- data_rel[rownames(meta_TH1),]
sums_1 <- colSums(data_TH1)
# sums_1 <- sums_1/sum(sums_1)


data_TH7 <- data_rel[rownames(meta_TH7),]
sums_7 <- colSums(data_TH7)
# sums_7 <- sums_7/sum(sums_7)

total <- data.frame(row.names = names(sums_1), TH1 = sums_1, TH7 = sums_7, Ranking_1 = rank(-sums_1), Ranking_7 = rank(-sums_7))

total[indval_res$X,"indval_TH"] <- indval_res$index
total$indval_TH[is.na(total$indval_TH)] = 0

p <- ggplot(total, aes(x = Ranking_1, shape = factor(indval_TH))) + 
  geom_point(aes(y = TH1, color = "TH1"), alpha = 0.75) + 
  geom_point(aes(y = TH7, color = "TH7"), alpha = 0.75) +
  scale_color_manual(
    name = "Tree hole", 
    values = c("TH1" = "blue", "TH7" = "red") )
save_without_legend(p, file.path(fig.dir, "Ranking_TH1.pdf"), width = 15, height = 12,
                    legend_file = file.path(fig.dir, "Ranking_TH_legend.pdf"))

p <- ggplot(total, aes(x = Ranking_7, shape = factor(indval_TH))) + 
  geom_point(aes(y = TH1, color = "TH1"), alpha = 0.75) + 
  geom_point(aes(y = TH7, color = "TH7"), alpha = 0.75) +
  scale_color_manual(
    name = "Tree hole", 
    values = c("TH1" = "blue", "TH7" = "red") )
save_without_legend(p, file.path(fig.dir, "Ranking_TH7.pdf"), width = 15, height = 12,
                    legend_file = file.path(fig.dir, "Ranking_TH_legend.pdf"))

p <- ggplot(total, aes(x = Ranking_1, y = Ranking_7, shape = factor(indval_TH))) + 
  geom_point() 
save_without_legend(p, file.path(fig.dir, "RankingvRanking.pdf"), width = 15, height = 12,
                    legend_file = file.path(fig.dir, "Ranking_Ind_legend.pdf"))


total$Dif_Rank <- (total$Ranking_1 - total$Ranking_7)
total$Dif_Abun <- (total$TH1 - total$TH7)

p <- ggplot(total, aes(x = Dif_Rank, y = Dif_Abun, shape = factor(indval_TH))) + 
  geom_point() 
save_without_legend(p, file.path(fig.dir, "DifRankingvsDifAbun.pdf"), width = 15, height = 12,
                    legend_file = file.path(fig.dir, "Ranking_Ind_legend.pdf"))

total[((total$Dif_Rank > 1000) | (total$Dif_Abun < -0.05)), "ranking_TH"] = "TH7"
total[((total$Dif_Rank < -1000) | (total$Dif_Abun > 0.05)), "ranking_TH"] = "TH1"
total[is.na(total$ranking_TH), "ranking_TH"] = "None"

total$indval_TH = ifelse(total$indval_TH == 0, "None", total$indval_TH)
total$indval_TH = ifelse(total$indval_TH == 2, "TH7", total$indval_TH)
total$indval_TH = ifelse(total$indval_TH == 1, "TH1", total$indval_TH)

p <- ggplot(total, aes(x = Dif_Rank, 
                  y = Dif_Abun, 
                  shape = factor(indval_TH), 
                  color = factor(ranking_TH))) + 
  # Check whether indval_TH is not 0, for the point size
  geom_point(aes(size = indval_TH != "None")) + 
  
  # Set the sizes (FALSE = normal, TRUE = larger) and hide this legend
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 5), guide = "none") +
  
  # Change the axis labels and the legend titles
  labs(
    x = "Difference in Rankings (TH1 - TH7)",
    y = "Difference in relative abundance (TH1 - TH7)",
    color = "Selected by ranking/abundance",
    shape = "Selected by Indicator Value"
  ) +
  theme_bw() +
  
  # Set the text size and the legend box
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    legend.position = "none"
  ) +
  
  # Force the point size inside the legend only
  guides(
    color = guide_legend(override.aes = list(size = 10)),
    shape = guide_legend(override.aes = list(size = 10))
  )
save_without_legend(p, file.path(fig.dir, "DifRankingvsDifAbun_thres.pdf"), width = 15, height = 12)


total$ASVs = rownames(total)
write.table(total, file.path(ind.dir, "characteristic_ASVs_M7.tsv"), sep = "\t", quote = F, row.names = F, col.names = T )


### TOTAL
rank_thres <- data.frame(Thres = 0, n = nrow(total))
for (i in 0:max(abs(total$Dif_Rank)))
  rank_thres = rbind(rank_thres, data.frame(Thres = i, n = sum(total$Dif_Rank > i)))  

p <- ggplot(rank_thres) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
save_without_legend(p, file.path(fig.dir, "ranking_threshold_total.pdf"), width = 15, height = 12)
#### NOT TOTAL ###

total_pos <- total %>% filter(Dif_Rank >= 0)
total_neg <- total %>% filter(Dif_Rank <= 0)

rank_thres_pos <- data.frame(Thres = 0, n = sum(total$Dif_Rank >= 0))
rank_thres_neg <- data.frame(Thres = 0, n = sum(total$Dif_Rank <= 0))
for (i in 0:(max(abs(total$Dif_Rank))-1))
{rank_thres_pos = rbind(rank_thres_pos, data.frame(Thres = i, n = sum(total_pos$Dif_Rank > i)))
rank_thres_neg = rbind(rank_thres_neg, data.frame(Thres = - i, n = sum(abs(total_neg$Dif_Rank) > i))) 
}  

p <- ggplot(rank_thres_pos) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
save_without_legend(p, file.path(fig.dir, "ranking_threshold_positive.pdf"), width = 15, height = 12)
p <- ggplot(rank_thres_neg) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
save_without_legend(p, file.path(fig.dir, "ranking_threshold_negative.pdf"), width = 15, height = 12)



Rank_coarse_1 <- total %>%
  
  arrange(Ranking_1) %>%
  mutate(Ranking = (Ranking_1 - 1) %/% 5 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_1), "-", max(Ranking_1)),
    TH1 = mean(TH1),
    TH7 = mean(TH7),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_1, aes(x = Ranking)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH7, color = "red"))

Rank_coarse_7 <- total %>%
  
  arrange(Ranking_7) %>%
  mutate(Ranking = (Ranking_7 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_7), "-", max(Ranking_7)),
    TH1 = mean(TH1),
    TH7 = mean(TH7),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_7, aes(x = Ranking)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH7, color = "red"))



