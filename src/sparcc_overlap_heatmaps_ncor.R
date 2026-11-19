###############################################################################
# Overlap heatmaps of the metacommunity networks for each number of links N
#
# Same overlap as src/sparcc_network_overlap.R, keeping the N strongest links.
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlations.tsv
#   4.3_Niche_Interactions/thresholds/Thres_cors.tsv
# Output:
#   4.3_Niche_Interactions/figures/thresholds/Heatmap_overlap_<N>.pdf
#
# Known issue: links are filtered with abs(zscore) >= Threshold, but Thres_cors.tsv
# holds thresholds of the SparCC score (see README).
# Run from anywhere inside the repository: Rscript src/sparcc_overlap_heatmaps_ncor.R
###############################################################################
library(here)
fig.dir <- here("4.3_Niche_Interactions", "figures", "thresholds")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)
library(gplots)
library(viridisLite)

# --- 1. READING AND PRE-COMPUTATION OUTSIDE THE LOOP ---
# Vectorise the z-score so it is not computed N times
cluster_com <- read_tsv(here("4.3_Niche_Interactions", "correlations", "Correlations.tsv")) %>%
        filter(method == "SparCC") %>%
        group_by(MetaCom) %>%
        mutate(zscore = scale(value)[, 1]) %>% 
        ungroup()

thres_all <- read_tsv(here("4.3_Niche_Interactions", "thresholds", "Thres_cors.tsv"))
target_ncors <- unique(thres_all$n_cor)

# --- 2. FIXED HEATMAP SETTINGS ---
hclust.av <- function(x) { hclust(x, method = "average") }

hm_lmat <- rbind(c(4, 3), c(2, 1))
hm_lwid <- c(0.8, 4)
hm_lhei <- c(1.2, 4)

# --- 3. MAIN LOOP OVER THRESHOLDS ---
for (nas in target_ncors) {
        
        # Vectorised filtering with the thresholds
        thres <- thres_all %>% 
                filter(n_cor == nas) %>% 
                mutate(MetaCom_join = paste0("Meta_", Metacom))
        
        current_data <- cluster_com %>%
                inner_join(thres %>% select(MetaCom_join, Threshold), by = c("MetaCom" = "MetaCom_join")) %>%
                filter(round(abs(zscore), 6) >= round(Threshold, 6))
        
        ASVs <- unique(c(current_data$taxa_1, current_data$taxa_2))
        n_asv <- length(ASVs)
        
        if (n_asv == 0) next # Guard against empty thresholds
        
        metacoms_nas <- paste0("Meta_", unique(thres$Metacom))
        n_meta <- length(metacoms_nas)
        
        # --- Fast construction of the co-occurrence networks ---
        networks_coex <- list()
        for (met in metacoms_nas) {
                # This is the fix: the ASVs are passed as row and column names.
                mat <- matrix(0, ncol = n_asv, nrow = n_asv, dimnames = list(ASVs, ASVs))
                
                df_met <- current_data %>% filter(MetaCom == met)
                if (nrow(df_met) > 0) {
                        mat[as.matrix(df_met[, c("taxa_1", "taxa_2")])] <- 1 
                }
                networks_coex[[met]] <- mat
        }
        
        # --- CVC and NCONS matrices ---
        cvc <- matrix(0, ncol = n_meta, nrow = n_meta, dimnames = list(metacoms_nas, metacoms_nas))
        ncons <- matrix(0, ncol = n_meta, nrow = n_meta, dimnames = list(metacoms_nas, metacoms_nas))
        
        # Pre-compute the sums to speed up Szymkiewicz–Simpson
        sums_coex <- vapply(networks_coex, sum, numeric(1))
        
        for (i in metacoms_nas) {
                for (j in metacoms_nas) {
                        intersec <- sum(networks_coex[[i]] * networks_coex[[j]])
                        ncons[i, j] <- intersec / 2
                        cvc[i, j] <- intersec / min(sums_coex[i], sums_coex[j])
                }
        }
        
        # --- HEATMAP ---
        threshold <- 0.4 
        text_colors <- ifelse(cvc > threshold, "black", "white")
        
        # Dummy run to extract the dendrogram (the price of using gplots)
        hm_temp <- heatmap.2(cvc, trace = "none", hclustfun = hclust.av)
        text_colors_reordered <- text_colors[hm_temp$rowInd, hm_temp$colInd]
        
        pdf(file.path(fig.dir, paste0("Heatmap_overlap_", nas, ".pdf")), height = 10, width = 15)
        heatmap.2(cvc,
                  col = inferno(100),            
                  density.info = "none",
                  trace = "none",
                  hclustfun = hclust.av, 
                  cellnote = as.matrix(ncons), 
                  notecol = text_colors_reordered, 
                  notecex = 1.2,                    
                  keysize = 0.5,                
                  lmat = hm_lmat,               
                  lwid = hm_lwid,               
                  lhei = hm_lhei,               
                  key.title = NA,               
                  key.xlab = "Overlap",         
                  key.par = list(cex.lab = 2.0, cex.axis = 1.5, mar = c(6, 2, 2, 2)), 
                  cexRow = 1.8,                 
                  cexCol = 1.8,                 
                  margins = c(12, 12),          
                  colsep = 1:ncol(cvc), 
                  rowsep = 1:nrow(cvc),
                  sepcolor = "white",           
                  sepwidth = c(0.02, 0.02)      
        )
        dev.off()
}