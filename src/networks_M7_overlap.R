###############################################################################
# Overlap between the M7 network and the TH1/TH7 networks
#
# Szymkiewicz-Simpson overlap between the networks of networks_M7.R; cells show
# the number of links in common.
#   Heatmap_overlap.pdf: all the networks (M7, M7 positive links, M7 restricted
#     to the characteristic ASVs, TH1 + TH7, TH1, TH7, the same restricted to the
#     characteristic ASVs). Links are counted once per direction, so cells show
#     links / 2 (as in the original script).
#   Heatmap_overlap_reduced.pdf: only the networks restricted to the
#     characteristic ASVs (M7, TH1 + TH7, TH1, TH7).
#
# Input:
#   4.3_Niche_Interactions/networks_M7/{network_7,network_7_Indval,network_holes}.tsv
#   4.3_Niche_Interactions/indicator_species/characteristic_ASVs_M7.tsv
# Output:
#   4.3_Niche_Interactions/figures/networks_M7/Heatmap_overlap.pdf
#   4.3_Niche_Interactions/figures/networks_M7/Heatmap_overlap_reduced.pdf
#
# Known issue: the Hole_1 and Hole_7 labels are swapped in both heatmaps (see README).
# Run from anywhere inside the repository: Rscript src/networks_M7_overlap.R
###############################################################################
library(here)
net.dir <- here("4.3_Niche_Interactions", "networks_M7")
fig.dir <- here("4.3_Niche_Interactions", "figures", "networks_M7")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)
library(gplots)
library(viridisLite)


## Characteristic ASVs of TH1 / TH7 (IndVal or ranking)
char_ASVs <- read.table(here("4.3_Niche_Interactions", "indicator_species", "characteristic_ASVs_M7.tsv"), header = T, sep = "\t")
# 0/None = not characteristic (see src/networks_M7.R)
char_ASVs <- char_ASVs %>% filter(!(indval_TH %in% c("0", "None")) | !(ranking_TH %in% c("0", "None")))

read_net <- function(f) read.table(file.path(net.dir, f), header = F, sep = "\t")

## Splits the merged TH1 + TH7 network into the links of each tree hole
## (the labels are swapped, as in the original script: see README)
split_holes <- function(cluster_com) {
        cluster_com$V4 <- str_remove(cluster_com$V4, "-.*")
        cluster_com_7 = cluster_com %>% filter(V4 %in% c("Hole_7", "Hole_1 + Hole_7"))
        cluster_com_1 = cluster_com %>% filter(V4 %in% c("Hole_1", "Hole_1 + Hole_7"))
        list(Hole_7 = cluster_com_1, Hole_1 = cluster_com_7)
}

## Overlap (Szymkiewicz-Simpson) and links in common between every pair of networks
overlap_matrices <- function(all, ncons_div = 1) {
        ASVs <- unique(unlist(lapply(all, function(x) c(x$V1, x$V2))))
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
                        ncons[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/ncons_div
                        cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/min(c(sum(abs(networks_coex[[i]])), sum(abs(networks_coex[[j]]))))                                      #Szymkiewicz–Simpson
                        # cvc[i,j] = sum(networks_coex[[i]]*networks_coex[[j]])/(sum(abs(networks_coex[[i]])) + sum(abs(networks_coex[[j]])) -sum(networks_coex[[i]]*networks_coex[[j]]) ) #Jaccard
                }
        list(cvc = cvc, ncons = ncons)
}

## Heatmap of the overlap with the number of links in common in each cell
plot_overlap_heatmap <- function(cvc, ncons, file, key_mar) {
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
        pdf(file, height = 10, width = 15)
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
                                 mar = key_mar),

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
        invisible(dev.off())
}


######## Full heatmap: all the networks ########
all <- list()

## THE NETWORK OF MET7
cluster_com <- read_net("network_7.tsv")
all[["Meta_7"]] <- cluster_com
all[["Meta_7_pos"]] <- cluster_com %>% filter(V3 > 0)

## THE NETWORK OF THE MET 7 BUT ONLY WITH THE INDVALUE ASVS
cluster_com <- read_net("network_7_Indval.tsv")
all[["Meta_7_indVal"]] <- cluster_com
all[["Meta_7_indVal_pos"]] <- cluster_com %>% filter(V3 > 0)

## THE NETWORK OF THE HOLES
cluster_com <- read_net("network_holes.tsv")
all[["Hole"]] <- cluster_com
all[c("Hole_7", "Hole_1")] <- split_holes(cluster_com)[c("Hole_7", "Hole_1")]

## THE NETWORK OF THE HOLES WITH THE INDVALUE ASVS
cluster_com <- cluster_com %>% filter((V1 %in% char_ASVs$ASVs) & (V2 %in% char_ASVs$ASVs))
all[["Hole_IndVal"]] <- cluster_com
holes <- split_holes(cluster_com)
all[["Hole_7_IndVal"]] <- holes$Hole_7
all[["Hole_1_IndVal"]] <- holes$Hole_1

ov <- overlap_matrices(all, ncons_div = 2)
plot_overlap_heatmap(ov$cvc, ov$ncons, file.path(fig.dir, "Heatmap_overlap.pdf"),
                     key_mar = c(6, 2, 2, 2))


######## Reduced heatmap: only the characteristic ASVs ########
all <- list()

## THE NETWORK OF THE MET 7 BUT ONLY WITH THE INDVALUE ASVS
all[["Meta_7"]] <- read_net("network_7_Indval.tsv")

## THE NETWORK THE HOLES WITH THE INDVALUE ASVS
cluster_com <- read_net("network_holes.tsv")
cluster_com <- cluster_com %>% filter((V1 %in% char_ASVs$ASVs) & (V2 %in% char_ASVs$ASVs))
all[["Hole"]] <- cluster_com
holes <- split_holes(cluster_com)
all[["Hole_7"]] <- holes$Hole_7
all[["Hole_1"]] <- holes$Hole_1

ov <- overlap_matrices(all)
colnames(ov$cvc) = rownames(ov$cvc) = colnames(ov$ncons) = rownames(ov$ncons) = c("Meta 7", "Hole 1 + Hole 7", "Hole 7", "Hole 1")
plot_overlap_heatmap(ov$cvc, ov$ncons, file.path(fig.dir, "Heatmap_overlap_reduced.pdf"),
                     key_mar = c(6, 2, 8, 2))
