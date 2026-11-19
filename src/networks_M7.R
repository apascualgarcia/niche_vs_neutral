###############################################################################
# Networks of metacommunity 7 and of tree holes 1 and 7 (|SparCC| >= 0.3)
#
# network_7: M7 network; network_7_Indval: restricted to the tree-hole
# characteristic ASVs (IndVal or ranking); network_holes: TH1 and TH7 networks
# merged (type = tree hole(s) and sign). Format: functionInk edge list.
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlation_zscore_thres{,_holes}.tsv
#   4.3_Niche_Interactions/indicator_species/characteristic_ASVs_M7.tsv
# Output:
#   4.3_Niche_Interactions/networks_M7/{network_7,network_7_Indval,network_holes}.tsv
# Run from anywhere inside the repository: Rscript src/networks_M7.R
###############################################################################
library(here)
net.dir <- here("4.3_Niche_Interactions", "networks_M7")
dir.create(net.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)

remove_bidirectional <- function(df) {
  df %>%
    mutate(
      # pmin and pmax compare the elements pairwise and return the alphabetically smaller/larger one
      t1_temp = pmin(taxa_1, taxa_2),
      t2_temp = pmax(taxa_1, taxa_2)
    ) %>%
    # Overwrite the original columns with the standardised order
    mutate(
      taxa_1 = t1_temp,
      taxa_2 = t2_temp
    ) %>%
    # Drop the temporary columns
    select(-t1_temp, -t2_temp) %>%
    # distinct removes the duplicated rows based on the taxa and the metacommunity.
    # .keep_all = TRUE keeps the 'value' column
    distinct(taxa_1, taxa_2, MetaCom, .keep_all = TRUE)
}

cors_met <- read.table(here("4.3_Niche_Interactions", "correlations", "Correlation_zscore_thres.tsv"), header = T, sep = "\t")
cors_holes <- read.table(here("4.3_Niche_Interactions", "correlations", "Correlation_zscore_thres_holes.tsv"), header = T, sep = "\t")


cors_met_7 <- cors_met %>% filter(MetaCom == "Meta_7" & method == "SparCC" & abs(value) >= 0.3)

cors_met_7 <- remove_bidirectional(cors_met_7)

cors_hole_1 <- cors_holes %>% filter(MetaCom == "Hole_1" & method == "SparCC" & abs(value) >= 0.3)
cors_hole_1 <- remove_bidirectional(cors_hole_1)

cors_hole_7 <- cors_holes %>% filter(MetaCom == "Hole_7" & method == "SparCC" & abs(value) >= 0.3)
cors_hole_7 <- remove_bidirectional(cors_hole_7)

#We merge the Species on both holes with the same sign
df_final <- bind_rows(cors_hole_1, cors_hole_7) %>%
  # Group by the taxa and by the sign of the 'value' column
  group_by(taxa_1, taxa_2, link_sign = sign(value)) %>%
  summarise(
    # If there is more than one (same sign and taxa), join them with " + "
    # If there is only one, paste() with collapse leaves the original string unchanged
    MetaCom = paste(MetaCom, collapse = " + "),
    
    # Take the mean. For a single value, the mean is the value itself
    value = mean(value),
    
    .groups = "drop"
  ) %>%
  # Drop the helper sign column
  select(-link_sign)


## We filter for functionINK
cors_met_7$type <- ifelse(cors_met_7$value > 0, "P", "N" )
cors_met_7_network <- cors_met_7[, c("taxa_1", "taxa_2", "value", "type")]

df_final$type <- ifelse(df_final$value > 0, "P", "N" )
df_final$type<- paste(df_final$MetaCom, df_final$type, sep = "-")
df_final_network <- df_final[, c("taxa_1", "taxa_2", "value", "type")]

df_final_network <- df_final_network %>%
  # Group by the taxa and by the sign of the 'value' column
  group_by(taxa_1, taxa_2) %>%
  summarise(
    # Take the mean. For a single value, the mean is the value itself
    value = mean(abs(value)),
    # If there is more than one (taxa), join them with " + "
    # If there is only one, paste() with collapse leaves the original string unchanged
    type = paste(type, collapse = " + "),
    
    .groups = "drop"
  ) 



colnames(cors_met_7_network) = colnames(df_final_network) = c("#SpeciesA", "SpeciesB", "interaction", "type")


write.table(df_final_network, file.path(net.dir, "network_holes.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(cors_met_7_network, file.path(net.dir, "network_7.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)





#####Now, I want to filter from the network_7 only the species that appear as indicator values #####


char_ASVs <- read.table(here("4.3_Niche_Interactions", "indicator_species", "characteristic_ASVs_M7.tsv"), header = T, sep = "\t")

# "None" = not characteristic. Files written by older versions of src/ranking_treeholes_M7.R
# used 0 instead; both are accepted.
char_ASVs_filter <- char_ASVs %>% filter(!(ranking_TH %in% c("0", "None")) | !(indval_TH %in% c("0", "None")))

ASVs <- char_ASVs_filter$ASVs

cors_met_7_network_IndVal <- cors_met_7_network %>% filter(`#SpeciesA` %in% ASVs & SpeciesB %in% ASVs)

write.table(cors_met_7_network_IndVal, file.path(net.dir, "network_7_Indval.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



