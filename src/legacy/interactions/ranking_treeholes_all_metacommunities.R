library(indicspecies)
library(dplyr)
library(vegan)
library(gplots)
library(tidyverse)

data <- t(read.table(file = "../otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = "../Metadata_singleTree_wDirich-K10.tsv", header = T, row.names = 1, sep = "\t")

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

data <- data[rownames(meta)[rownames(meta) %in% rownames(data)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

##### MET 7 ########


ind_val <- read.csv("../Ind_Val_Holes_prueba/indval_M7.csv")
meta_metacom <- meta %>% filter(Dirichlet == 7)

data_metacom <- data_rel[rownames(meta_metacom),]
data_metacom <- data_metacom[,colSums(data_metacom > 0) == nrow(data_metacom)]

## MY INDVAL

meta_th1 <- meta_metacom %>% filter(Location == 1)
meta_th7 <- meta_metacom %>% filter(Location == 7)

data_th1 <- data_rel[rownames(meta_th1),]
sums_1 <- colSums(data_th1)
# sums_1 <- sums_1/sum(sums_1)


data_th7 <- data_rel[rownames(meta_th7),]
sums_7 <- colSums(data_th7)
# sums_7 <- sums_7/sum(sums_7)

total <- data.frame(row.names = names(sums_1), TH1 = sums_1, TH7 = sums_7, Ranking_1 = rank(-sums_1), Ranking_7 = rank(-sums_7))

total[ind_val$X,"indval_TH"] <- ind_val$index
total$indval_TH[is.na(total$indval_TH)] = 0

pdf("Ranking_TH1.pdf", height = 12, width = 15)
ggplot(total, aes(x = Ranking_1, shape = factor(indval_TH))) + 
  geom_point(aes(y = TH1, color = "TH1"), alpha = 0.75) + 
  geom_point(aes(y = TH7, color = "TH7"), alpha = 0.75) +
  scale_color_manual(
    name = "Variable legend", 
    values = c("TH1" = "blue", "TH7" = "red") )
dev.off()

pdf("Ranking_TH7.pdf", height = 12, width = 15)
ggplot(total, aes(x = Ranking_7, shape = factor(indval_TH))) + 
  geom_point(aes(y = TH1, color = "TH1"), alpha = 0.75) + 
  geom_point(aes(y = TH7, color = "TH7"), alpha = 0.75) +
  scale_color_manual(
    name = "Variable legend", 
    values = c("TH1" = "blue", "TH7" = "red") )
dev.off()

pdf("RankingvRanking.pdf", height = 12, width = 15)
ggplot(total, aes(x = Ranking_1, y = Ranking_7, shape = factor(indval_TH))) + 
  geom_point() 
dev.off()


total$Dif_Rank <- (total$Ranking_1 - total$Ranking_7)
total$Dif_Abun <- (total$TH1 - total$TH7)

pdf("DifRankingvsDifAbun.pdf", height = 12, width = 15)
ggplot(total, aes(x = Dif_Rank, y = Dif_Abun, shape = factor(indval_TH))) + 
  geom_point() 
dev.off()

total[((total$Dif_Rank > 1000) | (total$Dif_Abun < -0.05)), "ranking_TH"] = "TH7"
total[((total$Dif_Rank < -1000) | (total$Dif_Abun > 0.05)), "ranking_TH"] = "TH1"
total[is.na(total$ranking_TH), "ranking_TH"] = "None"

total$indval_TH = ifelse(total$indval_TH == 0, "None", total$indval_TH)
total$indval_TH = ifelse(total$indval_TH == 2, "TH7", total$indval_TH)
total$indval_TH = ifelse(total$indval_TH == 1, "TH1", total$indval_TH)

pdf("DifRankingvsDifAbun_thres.pdf", height = 12, width = 15)
ggplot(total, aes(x = Dif_Rank, 
                  y = Dif_Abun, 
                  shape = factor(indval_TH), 
                  color = factor(ranking_TH))) + 
  # Check whether indval_TH is not "None", for the point size
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
  
dev.off()


total$ASVs = rownames(total)
write.table(total, "characteristic_ASVs_M7.tsv", sep = "\t", quote = F, row.names = F, col.names = T )


### TOTAL
rank_thres <- data.frame(Thres = 0, n = nrow(total))
for (i in 0:max(abs(total$Dif_Rank)))
  rank_thres = rbind(rank_thres, data.frame(Thres = i, n = sum(total$Dif_Rank > i)))  

pdf("ranking_threshold_total.pdf", height = 12, width = 15)
ggplot(rank_thres) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
dev.off()
#### NOT TOTAL ###

total_pos <- total %>% filter(Dif_Rank >= 0)
total_neg <- total %>% filter(Dif_Rank <= 0)

rank_thres_pos <- data.frame(Thres = 0, n = sum(total$Dif_Rank >= 0))
rank_thres_neg <- data.frame(Thres = 0, n = sum(total$Dif_Rank <= 0))
for (i in 0:(max(abs(total$Dif_Rank))-1))
{rank_thres_pos = rbind(rank_thres_pos, data.frame(Thres = i, n = sum(total_pos$Dif_Rank > i)))
rank_thres_neg = rbind(rank_thres_neg, data.frame(Thres = - i, n = sum(abs(total_neg$Dif_Rank) > i))) 
}  

pdf("ranking_threshold_positive.pdf", height = 12, width = 15)
ggplot(rank_thres_pos) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
dev.off()
pdf("ranking_threshold_negative.pdf", height = 12, width = 15)
ggplot(rank_thres_neg) + geom_point(aes(x = Thres, y = n)) + scale_y_log10()
dev.off()



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



######### MET 10 ##########

meta_metacom <- meta %>% filter(Dirichlet == 10)

data_metacom <- data_rel[rownames(meta_metacom),]
data_metacom <- data_metacom[,colSums(data_metacom > 0) == nrow(data_metacom)]

## MY INDVAL

meta_th8 <- meta_metacom %>% filter(Location == 8)
meta_th9 <- meta_metacom %>% filter(Location == 9)

data_th8 <- data_rel[rownames(meta_th8),]
sums_8 <- colSums(data_th8)

data_th9 <- data_rel[rownames(meta_th9),]
sums_9 <- colSums(data_th9)

total <- data.frame(row.names = names(sums_8), TH8 = sums_8, TH9 = sums_9, Ranking_8 = rank(-sums_8), Ranking_9 = rank(-sums_9))


ggplot(total, aes(x = Ranking_8)) + 
  geom_point(aes(y = TH8, color = "blue")) + 
  geom_point(aes(y = TH9, color = "red"))

ggplot(total, aes(x = Ranking_9)) + 
  geom_point(aes(y = TH8, color = "blue")) + 
  geom_point(aes(y = TH9, color = "red"))

ggplot(total, aes(x = Ranking_8, y = Ranking_9)) + 
  geom_point() 


Rank_coarse_8 <- total %>%
  
  arrange(Ranking_8) %>%
  mutate(Ranking = (Ranking_8 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_8), "-", max(Ranking_8)),
    TH8 = mean(TH8),
    TH9 = mean(TH9),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_8, aes(x = Ranking)) + 
  geom_point(aes(y = TH8, color = "blue")) + 
  geom_point(aes(y = TH9, color = "red"))

Rank_coarse_9 <- total %>%
  
  arrange(Ranking_9) %>%
  mutate(Ranking = (Ranking_9 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_9), "-", max(Ranking_9)),
    TH8 = mean(TH8),
    TH9 = mean(TH9),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_9, aes(x = Ranking)) + 
  geom_point(aes(y = TH8, color = "blue")) + 
  geom_point(aes(y = TH9, color = "red"))




######### MET 3 ##########

meta_metacom <- meta %>% filter(Dirichlet == 3)

data_metacom <- data_rel[rownames(meta_metacom),]
data_metacom <- data_metacom[,colSums(data_metacom > 0) == nrow(data_metacom)]

## MY INDVAL

meta_th1 <- meta_metacom %>% filter(Location == 1)
meta_th2 <- meta_metacom %>% filter(Location == 2)

data_th1 <- data_rel[rownames(meta_th1),]
sums_1 <- colSums(data_th1)

data_th2 <- data_rel[rownames(meta_th2),]
sums_2 <- colSums(data_th2)

total <- data.frame(row.names = names(sums_1), TH1 = sums_1, TH2 = sums_2, Ranking_1 = rank(-sums_1), Ranking_2 = rank(-sums_2))


ggplot(total, aes(x = Ranking_1)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH2, color = "red"))

ggplot(total, aes(x = Ranking_2)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH2, color = "red"))

ggplot(total, aes(x = Ranking_1, y = Ranking_2)) + 
  geom_point() 


Rank_coarse_1 <- total %>%
  
  arrange(Ranking_1) %>%
  mutate(Ranking = (Ranking_1 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_1), "-", max(Ranking_1)),
    TH1 = mean(TH1),
    TH2 = mean(TH2),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_1, aes(x = Ranking)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH2, color = "red"))

Rank_coarse_2 <- total %>%
  
  arrange(Ranking_2) %>%
  mutate(Ranking = (Ranking_2 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_2), "-", max(Ranking_2)),
    TH1 = mean(TH1),
    TH2 = mean(TH2),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_2, aes(x = Ranking)) + 
  geom_point(aes(y = TH1, color = "blue")) + 
  geom_point(aes(y = TH2, color = "red"))



######### MET 2 ##########

meta_metacom <- meta %>% filter(Dirichlet == 2)

data_metacom <- data_rel[rownames(meta_metacom),]
data_metacom <- data_metacom[,colSums(data_metacom > 0) == nrow(data_metacom)]

## MY INDVAL

meta_th10 <- meta_metacom %>% filter(Location == 10)
meta_th5 <- meta_metacom %>% filter(Location == 5)

data_th10 <- data_rel[rownames(meta_th10),]
sums_10 <- colSums(data_th10)

data_th5 <- data_rel[rownames(meta_th5),]
sums_5 <- colSums(data_th5)

total <- data.frame(row.names = names(sums_10), TH10 = sums_10, TH5 = sums_5, Ranking_10 = rank(-sums_10), Ranking_5 = rank(-sums_5))


ggplot(total, aes(x = Ranking_10)) + 
  geom_point(aes(y = TH10, color = "blue")) + 
  geom_point(aes(y = TH5, color = "red"))

ggplot(total, aes(x = Ranking_5)) + 
  geom_point(aes(y = TH10, color = "blue")) + 
  geom_point(aes(y = TH5, color = "red"))


ggplot(total, aes(x = Ranking_5, y = Ranking_10)) + 
  geom_point() 



Rank_coarse_5 <- total %>%
  
  arrange(Ranking_5) %>%
  mutate(Ranking = (Ranking_5 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_5), "-", max(Ranking_5)),
    TH5 = mean(TH5),
    TH10 = mean(TH10),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_5, aes(x = Ranking)) + 
  geom_point(aes(y = TH5, color = "blue")) + 
  geom_point(aes(y = TH10, color = "red"))

Rank_coarse_10 <- total %>%
  
  arrange(Ranking_10) %>%
  mutate(Ranking = (Ranking_10 - 1) %/% 10 + 1) %>%
  group_by(Ranking) %>%
  summarise(
    range = paste0("Top ", min(Ranking_10), "-", max(Ranking_10)),
    TH5 = mean(TH5),
    TH10 = mean(TH10),
    Elements = n() # Useful to check whether the last one has fewer than 10
  )

ggplot(Rank_coarse_10, aes(x = Ranking)) + 
  geom_point(aes(y = TH5, color = "blue")) + 
  geom_point(aes(y = TH10, color = "red"))























#IndVal

indval  <- multipatt(data_metacom, meta_metacom$Location, func = "IndVal.g", duleg = TRUE)

a <- indval$sign

b <- a[a$p.value <= 0.05,]

rowSums(data_metacom[,rownames(b)])/rowSums(data_metacom)

b$proper <- ifelse(b$index == 1, 1, 7)

a$Location <- b[rownames(a),"proper"]

write.csv(a,"indval_M7.csv", quote = F, row.names = T, col.names = T)

# Simper

simper <- simper(data_metacom, meta_metacom$Location)


# Wilcox 


results <- data.frame(Species = colnames(data_metacom), P_value = NA)

# Loop to run the Wilcoxon test for each species
for (i in 1:ncol(data_metacom)) {
  test <- wilcox.test(data_metacom[, i] ~ meta_metacom$Location)
  results$P_value[i] <- test$p.value
}

# Adjust p-values for multiple comparisons
results$P_adjusted <- p.adjust(results$P_value, method = "fdr")

# See which species differ significantly (adjusted p < 0.05)
significant_species <- subset(results, P_adjusted < 0.05)


significant_species$Species %in% rownames(b)
rownames(b) %in% significant_species$Species 


hclust.av = function(x){hclust(x,method = "average")}

rownames(data_metacom) = paste0( meta_metacom$Location, rownames(data_metacom))

heatmap.2(log(data_metacom),
          density.info = "none",
          trace = "none",
          hclustfun = hclust.av, 
          
          # --- KEY SETTINGS (thinner, large labels) ---
          keysize = 0.5,                # Smaller multiplier
          key.title = NA,               
          key.xlab = "Overlap",         
          key.par = list(cex.lab = 2.0,  # VERY large "Overlap" label
                         cex.axis = 1.5, # Large scale numbers
                         mar = c(6, 2, 2, 2)), # Larger bottom margin so the label fits
          
          # --- HEATMAP SETTINGS (larger labels) ---
          cexRow = 1.8,                 # Large row names
          cexCol = 1.8,                 # Large column names
          margins = c(12, 12),          # Larger margins so the large text is not cut off
          
          # --- CELL BORDERS ---
          colsep = 1:ncol(data_metacom), 
          rowsep = 1:nrow(data_metacom),
          sepcolor = "white",          
          sepwidth = c(0.02, 0.02)      # Slightly more visible borders
)
