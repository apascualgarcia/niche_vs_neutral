library(tidyverse)
library(gplots)
library(viridisLite)

cluster_com <- read.table("../6.Correlation_Distribution/Correlations.tsv", header = T, sep = "\t") 
cluster_com <- cluster_com %>% filter(method == "SparCC")

thres_all <- read.table("../8.Mean_Correlations_vs_samples/Thres_cors.tsv", header = T, sep = "\t")



all <- list()
ASVs <- c()
for (c in c("1","2","3","4","5","6","7","8","10"))
{
        cluster_com_met <- cluster_com %>% filter(MetaCom == paste0("Meta_",c))
        cluster_com_met$zscore = (cluster_com_met$value- mean(cluster_com_met$value))/ sd(cluster_com_met$value)
        
        cluster_com_met <- cluster_com_met %>% filter(abs(zscore) >= thres[thres$Metacom == c, "Threshold"])
        
        all[[c]] <- cluster_com_met
        ASVs <- c(ASVs, cluster_com_met$taxa_1, cluster_com_met$taxa_2)
}


ASVs <- unique(ASVs)


mat <- matrix(0, ncol = length(ASVs), nrow = length(ASVs))
colnames(mat) = rownames(mat) = ASVs

networks = list()


for (c in c("1","2","3","4","5","6","7","8","10"))
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


pdf("Heatmap_overlap.pdf", height = 10, width = 15)
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



