###############################################################################
# Tree hole(s) in which each ASV of the merged TH1-TH7 network has links
#
# Adds to the functionInk partition the tree holes of the links of each ASV, used
# in Cytoscape to split the network (Suppl. Fig., panel A).
#
# Input:
#   4.3_Niche_Interactions/networks_M7/network_holes.tsv
#   4.3_Niche_Interactions/networks_M7/functionInk/network_holes/Partition-NL_Average_StopStep-580_network_holes.tsv
# Output:
#   4.3_Niche_Interactions/networks_M7/functionInk/network_holes/Clustering_holes.tsv
# Run from anywhere inside the repository: Rscript src/network_treeholes_membership.R
###############################################################################
library(here)
net.dir <- here("4.3_Niche_Interactions", "networks_M7")
fink.dir <- here("4.3_Niche_Interactions", "networks_M7", "functionInk", "network_holes")
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)
### I want to create this script in order to simplify the representation of the holes. So I will check if the different Species have a connection from the hole 1 or hole 7 and give a factor. So, in cytoscape I can split the network.

clustering <- read.table(file.path(fink.dir, "Partition-NL_Average_StopStep-580_network_holes.tsv"), sep = "\t")
network <- read.table(file.path(net.dir, "network_holes.tsv"), header = F, sep = "\t")
## A couple of links appear in different communities with different signs. I originally removed them in network_holes_no_sign.tsv, but I will not take them into account here. They are only 2 nodes, so I will separate them by hand in the figure.
network$V4 <- str_remove(network$V4, "-.*")
clustering$Hole = ""
rownames(clustering) = clustering$V1

for(i in unique(network$V4))
{
  network_hole <- network %>% filter(V4 == i)
  ASVs <- unique( c(network_hole$V1, network_hole$V2))
  clustering[ASVs,"Hole"] <-  ifelse(clustering[ASVs,"Hole"] == "" , i,  paste0( clustering[ASVs,"Hole"], " | ", i))
  
  
}


colnames(clustering) = c("ASV", "Cluster", "Hole")
write.table(clustering, file.path(fink.dir, "Clustering_holes.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)


