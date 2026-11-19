###############################################################################
# Overlap between the SparCC networks of the metacommunities (Fig. 4A)
#
# Szymkiewicz-Simpson overlap of the links (all significant links, p <= 0.001)
# between every pair of metacommunity networks; cells show the links in common.
# Same for the tree-hole networks.
#
# Input:
#   4.3_Niche_Interactions/SparCC/{Meta_*,Hole_*,Total}/correlations.edgelist
# Output:
#   4.3_Niche_Interactions/figures/correlations/Heatmap_overlap.pdf (Fig. 4A)
#   4.3_Niche_Interactions/figures/correlations/Heatmap_overlap_Holes.pdf
# Run from anywhere inside the repository: Rscript src/sparcc_network_overlap.R
###############################################################################
library(here)
sparcc.dir <- here("4.3_Niche_Interactions", "SparCC")
fig.dir <- here("4.3_Niche_Interactions", "figures", "correlations")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)
library(gplots)
library(viridisLite)


all <- list()
ASVs <- c()
for (c in c("1","2","3","4","5","6","7","8","10"))
{
        cluster_com <- read.table(file.path(sparcc.dir, paste0("Meta_", c), "correlations.edgelist"), header = T, sep = "\t") 
        cluster_com$Com = c
        all[[c]] <- cluster_com
        ASVs <- c(ASVs, cluster_com$Sp1, cluster_com$Sp2)
}

cluster_com <- read.table(file.path(sparcc.dir, "Total", "correlations.edgelist"), header = T, sep = "\t") 
cluster_com$Com = "Total"
all[["Total"]] <- cluster_com
ASVs <- c(ASVs, cluster_com$Sp1, cluster_com$Sp2)


ASVs <- unique(ASVs)


mat <- matrix(0, ncol = length(ASVs), nrow = length(ASVs))
colnames(mat) = rownames(mat) = ASVs

networks = list()


for (c in c("1","2","3","4","5","6","7","8","10"))
{
        networks[[c]] = mat
        for (i in 1:nrow(all[[c]]))
                networks[[c]][all[[c]]$Sp1[i],all[[c]]$Sp2[i]] = all[[c]]$cor[i]   
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
                cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/min(c(sum(abs(networks_coex[[i]])), sum(abs(networks_coex[[j]]))))                                      #Szymkiewicz–Simpson
                # cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/(sum(abs(networks_coex[[i]])) + sum(abs(networks_coex[[j]])) -sum(networks_coex[[i]]*networks_coex[[j]]) ) #Jaccard
        }

colnames(cvc) = rownames(cvc) = paste0("Meta_", rownames(cvc))
colnames(ncons) = rownames(ncons) = paste0("Meta_", rownames(ncons))

library(viridisLite)


hm_lmat <- rbind(c(4, 3), c(2, 1))

# 2. lwid sets the column widths.
# The first value (0.8) makes the key and the left dendrogram very thin.
# The second value (4) gives all the space to the heatmap.
hm_lwid <- c(0.8, 4)

# 3. lhei sets the row heights.
hm_lhei <- c(1.2, 4)

hclust.av = function(x){hclust(x,method = "average")}


threshold <- 0.4 
text_colors <- ifelse(cvc > threshold, "black", "white")

hm_temp <- heatmap.2(cvc,
                     trace = "none",
                     hclustfun = hclust.av)

# 3. THE KEY STEP: reorder the colour matrix
# with the row and column indices computed by hm_temp
text_colors_reordered <- text_colors[hm_temp$rowInd, hm_temp$colInd]
pdf(file.path(fig.dir, "Heatmap_overlap.pdf"), height = 10, width = 15)
heatmap.2(cvc,
          col = inferno(100),            
          density.info = "none",
          trace = "none",
          hclustfun = hclust.av, 
          
          # Numbers inside the cells
          cellnote = as.matrix(ncons), 
          notecol = text_colors_reordered, 
          notecex = 1.8, # <-- Increase for larger numbers inside the cells
          
          # --- KEY SETTINGS ---
          keysize = 0.5,                
          lmat = hm_lmat,               
          lwid = hm_lwid,               
          lhei = hm_lhei,               
          key.title = NA,               
          key.xlab = "Overlap",         
          key.par = list(cex.lab = 2.5,  # <-- Legend title
                         cex.axis = 2.0, # <-- Legend axis
                         mar = c(6, 2, 8, 2)), 
          
          # --- HEATMAP SETTINGS (ROWS AND COLUMNS) ---
          cexRow = 2.5, # <-- Increase for the row names
          cexCol = 2.5, # <-- Increase for the column names
          margins = c(15, 15), # <-- The margin may need to grow if the font is very large
          
          # --- CELL BORDERS ---
          colsep = 1:ncol(cvc), 
          rowsep = 1:nrow(cvc),
          sepcolor = "white",           
          sepwidth = c(0.02, 0.02)      
)
dev.off()




all <- list()
ASVs <- c()
for (c in c("1","2","3","4","5","6","7","8","9", "10"))
{
        cluster_com <- read.table(file.path(sparcc.dir, paste0("Hole_", c), "correlations.edgelist"), header = T, sep = "\t") 
        cluster_com$Com = c
        all[[c]] <- cluster_com
        ASVs <- c(ASVs, cluster_com$Sp1, cluster_com$Sp2)
}

ASVs <- c(ASVs, cluster_com$Sp1, cluster_com$Sp2)


ASVs <- unique(ASVs)


mat <- matrix(0, ncol = length(ASVs), nrow = length(ASVs))
colnames(mat) = rownames(mat) = ASVs

networks = list()


for (c in c("1","2","3","4","5","6","7","8","9", "10"))
{
        networks[[c]] = mat
        for (i in 1:nrow(all[[c]]))
                networks[[c]][all[[c]]$Sp1[i],all[[c]]$Sp2[i]] = all[[c]]$cor[i]   
}

networks_coex = networks

for (c in names(networks))
        networks_coex[[c]] = ifelse(networks[[c]] == 0, 0, 1)




cvc <- matrix(0, ncol = length(names(networks_coex)), nrow = length(names(networks_coex)))
ncons <- matrix(0, ncol = length(names(networks_coex)), nrow = length(names(networks_coex)))

colnames(cvc) = rownames(cvc) = names(networks_coex)
colnames(ncons) = rownames(ncons) = names(networks_coex)
for (i in names(networks_coex))
        for (j in names(networks_coex)){
                ncons[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/2
                cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/min(c(sum(abs(networks_coex[[i]])), sum(abs(networks_coex[[j]]))))            }                          #Szymkiewicz–Simpson
# cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/(sum(abs(networks_coex[[i]])) + sum(abs(networks_coex[[j]])) -sum(networks_coex[[i]]*networks_coex[[j]]) ) #Jaccard


colnames(cvc) = rownames(cvc) = paste0("Hole_", rownames(cvc))
colnames(ncons) = rownames(ncons) = paste0("Hole_", rownames(ncons))

library(viridisLite)


hm_lmat <- rbind(c(4, 3), c(2, 1))

# 2. lwid sets the column widths.
# The first value (0.8) makes the key and the left dendrogram very thin.
# The second value (4) gives all the space to the heatmap.
hm_lwid <- c(0.8, 4)

# 3. lhei sets the row heights.
hm_lhei <- c(1.2, 4)

threshold <- 0.4 
text_colors <- ifelse(cvc > threshold, "black", "white")

hclust.av = function(x){hclust(x,method = "average")}

hm_temp <- heatmap.2(cvc,
                     trace = "none",
                     hclustfun = hclust.av)

# 3. THE KEY STEP: reorder the colour matrix
# with the row and column indices computed by hm_temp
text_colors_reordered <- text_colors[hm_temp$rowInd, hm_temp$colInd]
pdf(file.path(fig.dir, "Heatmap_overlap_Holes.pdf"), height = 10, width = 15)
# 4. Now draw the final heatmap with the already ordered matrix
heatmap.2(cvc,
          col = inferno(100),            
          density.info = "none",
          trace = "none",
          hclustfun = hclust.av, 
          
          # Pass the data frame of numbers and the reordered colour matrix
          cellnote = as.matrix(ncons), 
          notecol = text_colors_reordered, 
          notecex = 1.2,                    
          
          # --- KEY SETTINGS ---
          keysize = 0.5,                
          lmat = hm_lmat,               
          lwid = hm_lwid,               
          lhei = hm_lhei,               
          key.title = NA,               
          key.xlab = "Overlap",         
          key.par = list(cex.lab = 2.0,  
                         cex.axis = 1.5, 
                         mar = c(6, 2, 8, 2)),    
          
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
