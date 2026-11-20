library(tidyverse)

library(robCompositions)
library(philentropy)


Shadow_Richness <- function(abun)
{
  
  effective_ASV = data.frame()
  abun[is.na(abun)] = 0
  abun_log = abun*log(abun)
  abun_log[is.na(abun_log)] = 0
  effective_ASV = exp(-rowSums(abun_log))
  
  
  mean_Seff = mean(effective_ASV)
  Error_Seff = sd(effective_ASV)/sqrt(length(effective_ASV))
  
  abun.m = as.matrix(abun)
  Sp_R = apply(abun.m, 1 , function(row) sum(row != 0))
  Species_Richness = mean(apply(abun.m, 1 , function(row) sum(row != 0)))
  eSpecies_Richness = sd(apply(abun.m, 1 , function(row) sum(row != 0)))/sqrt(ncol(abun))
  
  SR <- mean((Sp_R-effective_ASV)/Sp_R)
  eSR <- sd((Sp_R-effective_ASV)/Sp_R)
  
  return(SR)
}



ST = read.table("otu_table_treeholes_dada2_renamed_clean_red.tsv", header = T, sep = "\t")
ST_metadata = read.table("Metadata_singleTree_wDirich-K10.tsv", header = T, sep = "\t")

ST_metadata_filter <- ST_metadata[ST_metadata$Location %in% c(1:10),]


Samples <- colnames(ST)[colnames(ST) %in% c("Taxon",ST_metadata_filter$ID.Sample)]

ST_metadata_filter <- ST_metadata_filter[ST_metadata_filter$ID.Sample %in% Samples,]

ST_filter <- ST[,c("Taxon",ST_metadata_filter$ID.Sample)]


rownames(ST_filter) <- ST_filter$Taxon


# First, the beta and alpha diversity with all the samples

dist_JSD <- sqrt(JSD(as.matrix(t(ST_filter[,-1] / rowSums(ST_filter[,-1]))), unit = "log"))
edist_JSD = sd(dist_JSD)/sqrt(length(dist_JSD))
# dist_JSD = mean(dist_JSD)




########## HERE I SUBTRACT THE within from ALL ###########
dataset <- t(ST_filter[,-1])

metadata <- ST_metadata_filter
data_prob <- dataset / rowSums(dataset)

full_dist_matrix <- dist_JSD
rownames(full_dist_matrix) <- rownames(dataset)
colnames(full_dist_matrix) <- rownames(dataset)

get_mean_dist_between_groups <- function(dist_mat, meta) {
  groups <- unique(meta$Dirichlet)
  n_groups <- length(groups)
  
  # Empty matrix for the final results (group x group)
  res_matrix <- matrix(NA, nrow = n_groups, ncol = n_groups,
                       dimnames = list(groups, groups))
  between_group_dists <- c()
  for (i in groups) {
    for (j in groups) {
      # Skip pairs of the same group
      if (i == j) next
      
      # Get the sample IDs of each group
      ids_i <- meta$ID.Sample[meta$Dirichlet == i]
      ids_j <- meta$ID.Sample[meta$Dirichlet == j]
      
      # Extract the sub-block of the distance matrix
      sub_block <- dist_mat[ids_i, ids_j]
      
      # Mean of all the distances in that block
      res_matrix[i, j] <- mean(sub_block)
      between_group_dists = c(between_group_dists,sub_block)
    }
  }
  
  return(mean(between_group_dists))
}

dist_JSD <- get_mean_dist_between_groups(full_dist_matrix, metadata)



div_location <- data.frame(Location = "ALL", Beta = dist_JSD, Depth = 1, Volume = 1,
                Alpha = Shadow_Richness(as.matrix(ST_filter[,-1] / rowSums(ST_filter[,-1]))))

div_meta <- data.frame(Location = "ALL", Beta = dist_JSD, Depth = 1, Volume = 1,
                      Alpha = Shadow_Richness(as.matrix(ST_filter[,-1] / rowSums(ST_filter[,-1]))))

# Then, the beta and alpha diversity for each group

for (Meta in unique(ST_metadata_filter$Location))
{
  
  dist_JSD <- sqrt(JSD(as.matrix(t(ST_filter[,(ST_metadata_filter %>% filter(Location == Meta))$ID.Sample]) / rowSums(t(ST_filter[,(ST_metadata_filter %>% filter(Location == Meta))$ID.Sample]))), unit = "log"))
  
  Beta <- mean(dist_JSD)
  Alpha <- Shadow_Richness(as.matrix(t(ST_filter[,(ST_metadata_filter %>% filter(Location == Meta))$ID.Sample]) / rowSums(t(ST_filter[,(ST_metadata_filter %>% filter(Location == Meta))$ID.Sample]))))
  
  div_location[nrow(div_location) + 1 ,"Location"] = Meta
  div_location[nrow(div_location) ,"Beta"] = Beta
  div_location[nrow(div_location) ,"Alpha"] = Alpha
  div_location[nrow(div_location) ,"Depth"] = mean((ST_metadata_filter %>% filter(Location == Meta))$Depth)
  div_location[nrow(div_location) ,"Volume"] = mean((ST_metadata_filter %>% filter(Location == Meta))$Volume)
}


# Then, the beta and alpha diversity for each group

for (Meta in unique(ST_metadata_filter$Dirichlet))
{
  
  dist_JSD <- sqrt(JSD(as.matrix(t(ST_filter[,(ST_metadata_filter %>% filter(Dirichlet == Meta))$ID.Sample]) / rowSums(t(ST_filter[,(ST_metadata_filter %>% filter(Dirichlet == Meta))$ID.Sample]))), unit = "log"))
  
  Beta <- mean(dist_JSD)
  Alpha <- Shadow_Richness(as.matrix(t(ST_filter[,(ST_metadata_filter %>% filter(Dirichlet == Meta))$ID.Sample]) / rowSums(t(ST_filter[,(ST_metadata_filter %>% filter(Dirichlet == Meta))$ID.Sample]))))
  
  div_meta[nrow(div_meta) + 1 ,"Location"] = Meta
  div_meta[nrow(div_meta) ,"Beta"] = Beta
  div_meta[nrow(div_meta) ,"Alpha"] = Alpha
  div_meta[nrow(div_meta) ,"Depth"] = mean((ST_metadata_filter %>% filter(Dirichlet == Meta))$Depth)
  div_meta[nrow(div_meta) ,"Volume"] = mean((ST_metadata_filter %>% filter(Dirichlet == Meta))$Volume)
}


ggplot(div_meta) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Depth))

ggplot(div_meta) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Volume))


div_meta_notall <- div_meta %>% filter(Location != "ALL")

ggplot(div_meta_notall) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Depth))

ggplot(div_meta_notall) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Volume))



ggplot(div_location) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Depth))

ggplot(div_location) +
  geom_point(aes(x = Beta, y = Alpha, color = Location, size = Volume))

