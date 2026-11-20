library(tidyverse)
library(gplots)
library(viridisLite)


all <- list()
ASVs <- c()

## THE NETWORK OF MET7
cluster_com <- read.table(paste0("../Networks/network_7.tsv"), header = F, sep = "\t") 
cluster_com$Com = "Meta_7"
all[["Meta_7"]] <- cluster_com
ASVs <- c(ASVs, cluster_com$V1, cluster_com$V2)
cluster_com$Com = "Meta_7_pos"
cluster_com <- cluster_com %>% filter(V3 > 0)
all[["Meta_7_pos"]] <- cluster_com

## THE NETWORK OF THE MET 7 BUT ONLY WITH THE INDVALUE ASVS
cluster_com <- read.table(paste0("../Networks/network_7_Indval.tsv"), header = F, sep = "\t") 
cluster_com$Com = "Meta_7_indVal"
all[["Meta_7_indVal"]] <- cluster_com
ASVs <- c(ASVs, cluster_com$V1, cluster_com$V2)
cluster_com$Com = "Meta_7_indVal_pos"
cluster_com <- cluster_com %>% filter(V3 > 0)
all[["Meta_7_indVal_pos"]] <- cluster_com

## THE NETWORK OF THE HOLES
cluster_com <- read.table(paste0("../Networks/network_holes.tsv"), header = F, sep = "\t") 
cluster_com$Com = "Hole"
all[["Hole"]] <- cluster_com
ASVs <- c(ASVs, cluster_com$V1, cluster_com$V2)
cluster_com$V4 <- str_remove(cluster_com$V4, "-.*")
cluster_com_7 = cluster_com %>% filter(V4 %in% c("Hole_7", "Hole_1 + Hole_7"))
cluster_com_1 = cluster_com %>% filter(V4 %in% c("Hole_1", "Hole_1 + Hole_7"))

all[["Hole_7"]] <- cluster_com_1
all[["Hole_1"]] <- cluster_com_7

## THE NETWORK OF THE HOLES WITH THE INDVALUE ASVS
characteristic_ASVs <- read.table("../../5.Real_data/Single_tree/My_Ind_val/characteristic_ASVs_M7.tsv", header = T, sep = "\t")
characteristic_ASVs <- characteristic_ASVs %>% filter(indval_TH != 0 | ranking_TH != 0)


cluster_com <- read.table(paste0("../Networks/network_holes.tsv"), header = F, sep = "\t") 

cluster_com <- cluster_com %>% filter((V1 %in% characteristic_ASVs$ASVs) & (V2 %in% characteristic_ASVs$ASVs))


cluster_com$Com = "Hole"
all[["Hole_IndVal"]] <- cluster_com
ASVs <- c(ASVs, cluster_com$V1, cluster_com$V2)
cluster_com$V4 <- str_remove(cluster_com$V4, "-.*")
cluster_com_7 = cluster_com %>% filter(V4 %in% c("Hole_7", "Hole_1 + Hole_7"))
cluster_com_1 = cluster_com %>% filter(V4 %in% c("Hole_1", "Hole_1 + Hole_7"))

all[["Hole_7_IndVal"]] <- cluster_com_1
all[["Hole_1_IndVal"]] <- cluster_com_7






ASVs <- unique(ASVs)


mat <- matrix(0, ncol = length(ASVs), nrow = length(ASVs))
colnames(mat) = rownames(mat) = ASVs

networks = list()


for (c in names(all))
{
        networks[[c]] = mat
        for (i in 1:nrow(all[[c]]))
                networks[[c]][all[[c]]$V1[i],all[[c]]$V2[i]] = all[[c]]$V3[i]   
}

networks_coex = networks

for (c in names(networks))
        networks_coex[[c]] = ifelse(networks[[c]] == 0, 0, 1)


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



my_lmat <- rbind(c(4, 3), c(2, 1))

# 2. lwid sets the column widths.
# The first value (0.8) makes the key and the left dendrogram very thin.
# The second value (4) gives all the space to the heatmap.
my_lwid <- c(0.8, 4)

# 3. lhei sets the row heights.
my_lhei <- c(1.2, 4)

hclust.av = function(x){hclust(x,method = "average")}


threshold <- 0.4 
text_color_matrix <- ifelse(cvc > threshold, "black", "white")

hm_temp <- heatmap.2(cvc,
                     trace = "none",
                     hclustfun = hclust.av)

# 3. THE KEY STEP: reorder the colour matrix
# with the row and column indices computed by hm_temp
reordered_color_matrix <- text_color_matrix[hm_temp$rowInd, hm_temp$colInd]
pdf("Heatmap_overlap.pdf", height = 10, width = 15)
# 4. Now draw the final heatmap with the already ordered matrix
heatmap.2(cvc,
          col = inferno(100),            
          density.info = "none",
          trace = "none",
          hclustfun = hclust.av, 
          
          # Pass the data frame of numbers and the reordered colour matrix
          cellnote = as.matrix(ncons), 
          notecol = reordered_color_matrix, 
          notecex = 1.2,                    
          
          # --- KEY SETTINGS ---
          keysize = 0.5,                
          lmat = my_lmat,               
          lwid = my_lwid,               
          lhei = my_lhei,               
          key.title = NA,               
          key.xlab = "Overlap",         
          key.par = list(cex.lab = 2.0,  
                         cex.axis = 1.5, 
                         mar = c(6, 2, 2, 2)), 
          
          # --- HEATMAP SETTINGS ---
          cexRow = 1.8,                 
          cexCol = 1.8,                 
          margins = c(12, 12),          
          
          # --- CELL BORDERS ---
          colsep = 1:ncol(cvc), 
          rowsep = 1:nrow(cvc),
          sepcolor = "white",           
          sepwidth = c(0.02, 0.02)      
)
dev.off()

