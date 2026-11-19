###############################################################################
# SparCC score threshold that keeps the N strongest links of each metacommunity
#
# N = 1, 10, ..., 10^6 (all links when N exceeds half of them).
#
# Input:
#   4.3_Niche_Interactions/correlations/Correlations.tsv
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output:
#   4.3_Niche_Interactions/thresholds/Thres_cors.tsv
# Run from anywhere inside the repository: Rscript src/sparcc_threshold_table.R
###############################################################################
library(here)
thr.dir <- here("4.3_Niche_Interactions", "thresholds")
dir.create(thr.dir, recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

library(tidyverse)

### HERE I plot the number of correlations against the threshold, to cut at the right number.

cors <- read.table(here("4.3_Niche_Interactions", "correlations", "Correlations.tsv"), header = T, sep = "\t")
meta <- read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, sep = "\t", row.names = 1)
data <- t(read.table(file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"), header = T, row.names = 1, sep = "\t"))

data_rel <- data / rowSums(data)

data_rel <- data_rel[rownames(meta)[rownames(meta) %in% rownames(data_rel)],]
meta <- meta[rownames(data_rel)[rownames(data_rel) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))

summary_df <- data.frame()
ncor_thres <- data.frame()

for (i in c(1,2,3,4,5,6,7,8,10))
{
  data_meta_SparCC_filter <- cors %>% filter(method == "SparCC" & MetaCom == paste0("Meta_",i))
  data_meta_SparCC_filter$zscore = (data_meta_SparCC_filter$value- mean(data_meta_SparCC_filter$value))/ sd(data_meta_SparCC_filter$value)
  
  met_red <- meta %>% filter(Dirichlet == i)
  
  # 1. Sort the data frame by the absolute z-score, in decreasing order
  # data_meta_SparCC_filter <- data_meta_SparCC_filter %>% 
  #   arrange(desc(abs(zscore)))
  
  #Changed to sort by the value instead of the z-score, because the data are not normally distributed
  data_meta_SparCC_filter <- data_meta_SparCC_filter %>% 
    arrange(desc(abs(value)))
  
  # 2. Define the exact numbers of correlations to explore.
  # This list can be changed or extended (e.g. 1000, 1500, 2000, etc.)
  # target_ncors <- c(1, 10, 50, 100, 200, 300, 500, 1000, 1500, 2000, 2500, 3000, 5000, 7000, 8000, 10000, 11000, 12000, 13000, 14000, 15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 1000000)
  
  ## Reduced list
  target_ncors <- c(1, 10, 100, 1000, 10000, 100000, 200000, 300000, 500000, 700000, 800000, 900000, 1000000)
  for (n in target_ncors)
  {  
    # Safety check: if more correlations are requested than exist in this MetaCom, use all of them
    if (n < nrow(data_meta_SparCC_filter)/2){
    
    # Select exactly the top 'n' correlations with the largest absolute z-score
    cors_red <- data_meta_SparCC_filter[1:n, ]
    
    # The corresponding threshold is the absolute z-score at position 'n'
    # j <- abs(cors_red$zscore[n])    
    #Same as above
    j <- abs(cors_red$value[n])
    }
    else
    {
      cors_red <- data_meta_SparCC_filter
      # j <- abs(cors_red$zscore[nrow(cors_red)])
      #Same as above
      j <- abs(cors_red$value[nrow(cors_red)])
    }
    summary_df <- rbind(summary_df, data.frame(Metacom = i, 
                                     Mean = mean(abs(cors_red$value)), 
                                     nsamples = nrow(met_red), 
                                     n_cor = nrow(cors_red), # Will be exactly n
                                     se = sd(abs(cors_red$value)), 
                                     Mean_taxon = mean(rowSums(data_rel[rownames(met_red),]>0)),
                                     sd_taxon = sd(rowSums(data_rel[rownames(met_red),]>0)),
                                     Mean_possible_associations = mean(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                     sd_possible_associations = sd(rowSums(data_rel[rownames(met_red),]>0) * (rowSums(data_rel[rownames(met_red),]>0) -1)/2),
                                     Threshold = j
    ) 
    )
    
    # Store in ncor_thres the exact threshold (z-score) computed for that number 'n' of correlations
    ncor_thres <- rbind(ncor_thres, data.frame(Metacom = i,
                                               Threshold = j,
                                               n_cor = n))
  }
}



## Threshold

write.table(ncor_thres, file.path(thr.dir, "Thres_cors.tsv"), sep = "\t", row.names = F, col.names = T)

