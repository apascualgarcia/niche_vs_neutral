###############################################################################
# Overlap between metacommunity networks as a function of the number of links
#
# For each N in Thres_cors.tsv, keeps the N strongest links of each network and
# computes the Szymkiewicz-Simpson overlap of every pair of networks.
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlations.tsv
#   4.3_Niche_Interactions/thresholds/Thres_cors.tsv
# Output:
#   4.3_Niche_Interactions/figures/thresholds/ncors_overlap.pdf and ncors_overlap_legend.pdf
# Run from anywhere inside the repository: Rscript src/sparcc_overlap_vs_ncor.R
###############################################################################
library(here)
source(here("src", "legend_utils.R"))  # save_without_legend(), save_metacommunity_legend()
fig.dir <- here("4.3_Niche_Interactions", "figures", "thresholds")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)
# library(gplots) 
# library(viridisLite)

meta_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                  "#6D3F5B", "#1E213D", "#763C28")


# 1. Read the data and pre-compute the z-scores for the WHOLE dataset at once
cluster_com <- read_tsv(here("4.3_Niche_Interactions", "correlations", "Correlations.tsv")) %>%
        filter(method == "SparCC") %>%
        group_by(MetaCom) %>%
        mutate(zscore = scale(value)[, 1]) %>% # scale() vectorises (value - mean) / sd
        ungroup()

thres_all <- read_tsv(here("4.3_Niche_Interactions", "thresholds", "Thres_cors.tsv"))

target_ncors <- unique(thres_all$n_cor)

# List to store the results instead of calling rbind in every iteration
results_list <- list()

for (nas in target_ncors) {
        
        # Prepare the thresholds for this 'nas'
        thres <- thres_all %>%
                filter(n_cor == nas) %>%
                mutate(MetaCom_join = paste0("Meta_", Metacom))
        
        # 2. Filter cluster_com with an inner_join to apply the thresholds in a vectorised way
        current_data <- cluster_com %>%
                inner_join(thres %>% select(MetaCom_join, Threshold), 
                           by = c("MetaCom" = "MetaCom_join")) %>%
                filter(round(abs(value), 6) >= round(Threshold, 6))
                ## Changed so that links are ranked by correlation value instead of z-score, because the data are not normal
                # filter(round(abs(zscore), 6) >= round(Threshold, 6))
        
        # 3. Get the unique ASVs to set the matrix dimensions
        ASVs <- unique(c(current_data$taxa_1, current_data$taxa_2))
        n_asv <- length(ASVs)
        
        if (n_asv == 0) next # Avoid errors when no data pass the threshold
        
        metacoms_nas <- paste0("Meta_", unique(thres$Metacom))
        
        # 4. Build the (binary) co-occurrence networks
        networks_coex <- list()
        
        for (met in metacoms_nas) {
                # Initialise a matrix of zeros
                mat <- matrix(0, ncol = n_asv, nrow = n_asv)
                dimnames(mat) <- list(ASVs, ASVs)
                
                df_met <- current_data %>% filter(MetaCom == met)
                
                if (nrow(df_met) > 0) {
                        # Vectorised: matrix indexing with 2 columns. Replaces the inner for loop.
                        # Co-occurrence networks only need 0 or 1, so the 1s are set directly.
                        idx <- as.matrix(df_met[, c("taxa_1", "taxa_2")])
                        mat[idx] <- 1 
                }
                networks_coex[[met]] <- mat
        }
        
        # 5. Compute the overlaps efficiently over the pairs of networks
        meta_names <- names(networks_coex)
        
        if (length(meta_names) >= 2) {
                # Pre-compute the sums for the Szymkiewicz–Simpson index
                sums_coex <- vapply(networks_coex, sum, numeric(1)) 
                
                # combn gives only the unique pairs (equivalent to upper.tri)
                combos <- combn(meta_names, 2, simplify = FALSE)
                
                res_nas <- lapply(combos, function(pair) {
                        i <- pair[1]
                        j <- pair[2]
                        
                        # Intersection of two binary matrices
                        intersection <- sum(networks_coex[[i]] * networks_coex[[j]])
                        
                        data.frame(
                                Meta_1 = i,
                                Meta_2 = j,
                                overlap = intersection / min(sums_coex[i], sums_coex[j]), # Szymkiewicz–Simpson
                                overlap_counts = intersection / 2,                        # Ncons
                                ncors = nas,
                                stringsAsFactors = FALSE
                        )
                })
                
                results_list[[as.character(nas)]] <- bind_rows(res_nas)
        }
}

# 6. Combine all the results at the end (avoids repeated rbind)
final_results <- bind_rows(results_list)

# To get the two data frames in their original form:
df_all_overlap <- final_results %>% select(Meta_1, Meta_2, overlap, ncors)
df_all_counts  <- final_results %>% select(Meta_1, Meta_2, overlap = overlap_counts, ncors)

df_all_overlap$Meta_1 <- str_replace(df_all_overlap$Meta_1 , "_", " ")
df_all_overlap$Meta_2 <- str_replace(df_all_overlap$Meta_2 , "_", " ")

df_all_overlap <- df_all_overlap %>% filter(ncors != 10)



combo_palette <- c(
        "#AF3427", "#A06941", "#726472", "#866B53", "#CC8319", "#984B46", 
        "#703C37", "#9C4A2C", "#8C4638", "#5E4168", "#72474A", "#B85F0F", 
        "#84283C", "#5C192D", "#882623", "#507682", "#637C64", "#AA9429", 
        "#755D56", "#4E4E47", "#7A5B3D", "#357794", "#7C8F5A", "#475887", 
        "#204978", "#4C566D", "#8F963B", "#5B5E68", "#334F59", "#5F5D4F", 
        "#A1762E", "#7A671F", "#A67514", "#45304C", "#713D41", "#4A2E32"
)

p <- ggplot(df_all_overlap, aes(x = ncors/2, 
                           y = overlap, 
                           group = paste(Meta_1, Meta_2, sep = "-"), 
                           color = paste(Meta_1, Meta_2, sep = "-"))) + 
        geom_point(size = 2) + 
        geom_line(linewidth = 2) +
        scale_y_log10() +  
        # THIS WAS THE ERROR: "values =" was missing
        scale_color_manual(values = combo_palette) +
        labs(x = "Number of SparCC Correlations",
             y = "Overlap Coefficient between networks",
             color = "Pair of networks overlaped") + 
        theme_bw() + 
        theme(
                axis.title.x = element_text(size = 25), 
                axis.title.y = element_text(size = 25),
                axis.text.x = element_text(size = 20),  
                axis.text.y = element_text(size = 20),  
                legend.title = element_text(size = 25), 
                legend.text = element_text(size = 20),  
                legend.key.size = unit(1.5, "lines"),
                legend.position = "none"
        ) +
        guides(color = guide_legend(override.aes = list(size = 8)))
save_without_legend(p, file.path(fig.dir, "ncors_overlap.pdf"), width = 15, height = 10)

