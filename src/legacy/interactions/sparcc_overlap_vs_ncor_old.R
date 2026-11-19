library(tidyverse)
library(gplots)
library(viridisLite)

cluster_com <- read.table("../6.Correlation_Distribution/Correlations.tsv", header = T, sep = "\t") 
cluster_com <- cluster_com %>% filter(method == "SparCC")

thres_all <- read.table("../8.Mean_Correlations_vs_samples/Thres_cors.tsv", header = T, sep = "\t")

target_ncors <- unique(thres_all$n_cor)

df_all_overlap <- data.frame()
df_all_counts <- data.frame()
for (nas in target_ncors)
{
thres <- thres_all %>% filter(n_cor == nas)
        
        
all <- list()
ASVs <- c()
for (c in as.character(unique(thres$Metacom)))
{
        cluster_com_met <- cluster_com %>% filter(MetaCom == paste0("Meta_",c))
        cluster_com_met$zscore = (cluster_com_met$value- mean(cluster_com_met$value))/ sd(cluster_com_met$value)
        
        cluster_com_met <- cluster_com_met %>% filter(round(abs(zscore),digits = 6) >= round(thres[thres$Metacom == c, "Threshold"], digits = 6))
        
        all[[c]] <- cluster_com_met
        ASVs <- c(ASVs, cluster_com_met$taxa_1, cluster_com_met$taxa_2)
}


ASVs <- unique(ASVs)


mat <- matrix(0, ncol = length(ASVs), nrow = length(ASVs))
colnames(mat) = rownames(mat) = ASVs

networks = list()


for (c in as.character(unique(thres$Metacom)))
{
        networks[[c]] = mat
        for (i in 1:nrow(all[[c]]))
                networks[[c]][all[[c]]$taxa_1[i],all[[c]]$taxa_2[i]] = all[[c]]$value[i]   
}

networks_coex = networks

for (c in names(networks))
        networks_coex[[c]] = ifelse(networks[[c]] == 0, 0, 1)


# networks_coex = networks_coex[names(networks)[names(networks) != "Flashweave"]]


cvc <- matrix(0, ncol = length(names(networks_coex)), nrow = length(names(networks_coex)))
ncons <- matrix(0, ncol = length(names(networks_coex)), nrow = length(names(networks_coex)))
colnames(cvc) = rownames(cvc) = names(networks_coex)
colnames(ncons) = rownames(ncons) = names(networks_coex)
for (i in names(networks_coex))
        for (j in names(networks_coex))
        {
                ncons[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/2
                cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/min(c(sum(abs(networks_coex[[i]])), sum(abs(networks_coex[[j]]))))                                      #Szymkiewicz-Simpson
                # cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/(sum(abs(networks_coex[[i]])) + sum(abs(networks_coex[[j]])) -sum(networks_coex[[i]]*networks_coex[[j]]) ) #Jaccard
        }

colnames(cvc) = rownames(cvc) = paste0("Meta_", rownames(cvc))
colnames(ncons) = rownames(ncons) = paste0("Meta_", rownames(ncons))

ut <- upper.tri(cvc)

# 2. Build the data.frame by mapping the positions with row() and col()
df_over <- data.frame(
        Meta_1 = rownames(cvc)[row(cvc)[ut]],
        Meta_2 = colnames(cvc)[col(cvc)[ut]],
        overlap = cvc[ut],
        ncors = nas,
        stringsAsFactors = FALSE
)

df_all_overlap <- rbind(df_all_overlap, df_over)


ut <- upper.tri(ncons)

# 2. Build the data.frame by mapping the positions with row() and col()
df_counts <- data.frame(
        Meta_1 = rownames(ncons)[row(ncons)[ut]],
        Meta_2 = colnames(ncons)[col(ncons)[ut]],
        overlap = ncons[ut],
        ncors = nas,
        stringsAsFactors = FALSE
)

df_all_counts <- rbind(df_all_counts, df_counts)

}
