###############################################################################
# Split the ASV table into one table per metacommunity and per tree hole (SparCC input)
#
# One table per Dirichlet component (K = 10; components with a single sample are
# skipped), one per tree hole (all its samples, whatever their metacommunity) and
# the whole table (Total, Total_Holes). Counts are not rescaled.
#
# Input:
#   1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#   1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output:
#   4.3_Niche_Interactions/SparCC/{Meta_*,Hole_*,Total,Total_Holes}/data.tsv
# Run from anywhere inside the repository: Rscript src/sparcc_split_tables.R
###############################################################################
library(here)
sparcc.dir <- here("4.3_Niche_Interactions", "SparCC")
dir.create(file.path(sparcc.dir, "Total"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sparcc.dir, "Total_Holes"), recursive = TRUE, showWarnings = FALSE)
# Plots that the script does not save (exploratory) go to a null device under Rscript
if (!interactive()) pdf(NULL, width = 20, height = 20)

## Counts are not rescaled here.
library(tidyverse)

otu_file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv")
meta <- read.table(file = here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, row.names = 1, sep = "\t")


a <- read.table(file = otu_file, header = F, sep = "\t")
a[1,1] = "#OTU_ID"
colnames(a) = a[1,]
a = a[-1,]


a <- a[,c("#OTU_ID", rownames(meta)[rownames(meta) %in% colnames(a)])]
meta <- meta[colnames(a)[colnames(a) %in% rownames(meta)],]

meta <- meta %>% filter(Dirichlet %in% names(table(meta$Dirichlet)[table(meta$Dirichlet) > 1]))


for (i in unique(meta$Dirichlet[!is.na(meta$Dirichlet)]))
{
  data_meta <- a[,c("#OTU_ID",rownames(meta[meta$Dirichlet == i & !is.na(meta$Dirichlet), ]))]
  dir.create(file.path(sparcc.dir, paste0("Meta_", i)), showWarnings = FALSE)
  out_file = file.path(sparcc.dir, paste0("Meta_", i), "data.tsv")
  write.table(data_meta, file = out_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
}

data_meta <- a[,c("#OTU_ID",rownames(meta[!is.na(meta$Dirichlet), ]))]
out_file = file.path(sparcc.dir, "Total", "data.tsv")
write.table(data_meta, file = out_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


#### Now for treeholes

library(tidyverse)

otu_file = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv")
meta <- read.table(file = here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, row.names = 1, sep = "\t")


a <- read.table(file = otu_file, header = F, sep = "\t")
a[1,1] = "#OTU_ID"
colnames(a) = a[1,]
a = a[-1,]


a <- a[,c("#OTU_ID", rownames(meta)[rownames(meta) %in% colnames(a)])]
meta <- meta[colnames(a)[colnames(a) %in% rownames(meta)],]

meta <- meta %>% filter(Location %in% names(table(meta$Location)[table(meta$Location) > 1]))


for (i in unique(meta$Location[!is.na(meta$Location)]))
{
  data_meta <- a[,c("#OTU_ID",rownames(meta[meta$Location == i & !is.na(meta$Location), ]))]
  dir.create(file.path(sparcc.dir, paste0("Hole_", i)), showWarnings = FALSE)
  out_file = file.path(sparcc.dir, paste0("Hole_", i), "data.tsv")
  write.table(data_meta, file = out_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
}

data_meta <- a[,c("#OTU_ID",rownames(meta[!is.na(meta$Location), ]))]
out_file = file.path(sparcc.dir, "Total_Holes", "data.tsv")
write.table(data_meta, file = out_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
