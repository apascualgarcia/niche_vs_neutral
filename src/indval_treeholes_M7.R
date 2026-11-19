###############################################################################
# Indicator ASVs of tree holes 1 and 7 within metacommunity 7 (IndVal)
#
# indicspecies::multipatt (IndVal.g) on the ASVs present in all the M7 samples;
# Location = tree hole of the ASVs with p <= 0.05.
#
# Input:
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output:
#   4.3_Niche_Interactions/indicator_species/indval_M7.csv
# Run from anywhere inside the repository: Rscript src/indval_treeholes_M7.R
###############################################################################
library(here)
ind.dir <- here("4.3_Niche_Interactions", "indicator_species")
dir.create(ind.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(indicspecies)
library(dplyr)
library(vegan)
library(gplots)


data <- t(read.table(file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"), header = T, row.names = 1, sep = "\t"))
meta <- read.table(file = here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, row.names = 1, sep = "\t")

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

data <- data[rownames(meta)[rownames(meta) %in% rownames(data)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))


meta_M7 <- meta %>% filter(Dirichlet == 7)

data_M7 <- data_rel[rownames(meta_M7),]
data_M7 <- data_M7[,colSums(data_M7 > 0) == nrow(data_M7)]

#IndVal

indval  <- multipatt(data_M7, meta_M7$Location, func = "IndVal.g", duleg = TRUE)

a <- indval$sign

b <- a[a$p.value <= 0.05,]

rowSums(data_M7[,rownames(b)])/rowSums(data_M7)

b$proper <- ifelse(b$index == 1, 1, 7)

a$Location <- b[rownames(a),"proper"]

write.csv(a, file.path(ind.dir, "indval_M7.csv"), quote = F, row.names = T, col.names = T)

# Simper

simper <- simper(data_M7, meta_M7$Location)


# Wilcox 


wilcox_results <- data.frame(ASV = colnames(data_M7), p_value = NA)

# Loop to run the Wilcoxon test for each species
for (i in 1:ncol(data_M7)) {
  test <- wilcox.test(data_M7[, i] ~ meta_M7$Location)
  wilcox_results$p_value[i] <- test$p.value
}

# Adjust p-values for multiple comparisons
wilcox_results$p_adjusted <- p.adjust(wilcox_results$p_value, method = "fdr")

# See which species differ significantly (adjusted p < 0.05)
significant_ASVs <- subset(wilcox_results, p_adjusted < 0.05)


significant_ASVs$ASV %in% rownames(b)
rownames(b) %in% significant_ASVs$ASV 


hclust.av = function(x){hclust(x,method = "average")}

rownames(data_M7) = paste0( meta_M7$Location, rownames(data_M7))

heatmap.2(log(data_M7),
          density.info = "none",
          trace = "none",
          hclustfun = hclust.av, 
          
          # --- KEY SETTINGS (thinner, large labels) ---
          keysize = 0.5,                # Smaller multiplier
          key.title = NA,               
          key.xlab = "Overlap",         
          key.par = list(cex.lab = 2.0,  # VERY large "Overlap" label
                         cex.axis = 1.5, # Large scale numbers
                         mar = c(6, 2, 2, 2)), # Larger bottom margin so the label fits
          
          # --- HEATMAP SETTINGS (larger labels) ---
          cexRow = 1.8,                 # Large row names
          cexCol = 1.8,                 # Large column names
          margins = c(12, 12),          # Larger margins so the large text is not cut off
          
          # --- CELL BORDERS ---
          colsep = 1:ncol(data_M7), 
          rowsep = 1:nrow(data_M7),
          sepcolor = "white",          
          sepwidth = c(0.02, 0.02)      # Slightly more visible borders
)
