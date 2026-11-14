
rm(list = ls())
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(rprojroot)
library(ALDEx2)
library(data.table)

dir.root <- rprojroot::find_root(is_rstudio_project) # need to open the Rproject

# START EDITING ---------------------

dir.pred = "4.2_Niche_PiCRUSTsc/KO_metagenome_out"
dir.kegg = "4.2_Niche_PiCRUSTsc/KEGG_pathways_out"
dir.meta = "1_samples_properties"
dir.out = "4.2_Niche_PiCRUSTsc/analysis"

file.kegg = "KEGGs_vs_samples.L3.tsv" # if you already provide kegg pathways categorize_by_function.R output (ko2kegg = F)
file.pred = "pred_metagenome_unstrat.tsv" # if you need to transform prediction to kegg pathway (ko2kegg = T)
file.meta = "Metadata_singleTree_wDirich-K10.tsv"

ko2kegg = F # if you use file.kegg = F, if you use file.pred = T

# STOP EDITING -----------------------

setwd(dir.root)

# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.

# Read data and match tables ------
# Load metadata as a tibble
# data(metadata)
setwd(dir.meta)
metadata <- read_delim(file.meta, delim = "\t", escape_double = FALSE, trim_ws = TRUE)

metadata = metadata[metadata$Type == "DNA",]
metadata = metadata[!is.na(metadata$Dirichlet),]

# Load KEGG pathway abundance
# data(kegg_abundance)
setwd(dir.root)
if(ko2kegg == T){
  setwd(dir.pred)
  kegg_abundance_all <- ko2kegg_abundance(file.pred)
}else{
  setwd(dir.kegg)
  kegg_abundance_all <- read.csv(file.kegg, sep = "\t", row.names = 1,  check.names = FALSE, stringsAsFactors = FALSE)#, col.names = T)
}

# --- match both tables
colnames(kegg_abundance_all)[1:5]
metadata$ID.Sample[1:5]

matched = match(metadata$ID.Sample, colnames(kegg_abundance_all))
which(is.na(matched)==TRUE) # none
matched = match(colnames(kegg_abundance_all),metadata$ID.Sample)
which(is.na(matched)==TRUE) # none

# consider these lines if it does not return 0 above
#metadata = metadata[!is.na(matched),]
#matched = matched[!is.na(matched)]
#kegg_abundance = kegg_abundance_all[, matched]

#DAA ----------kegg_abundance-------------
# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance_all, 
                              metadata = metadata, 
                              group = "Dirichlet", 
                              daa_method = "ALDEx2", 
                              include_effect_size = TRUE,
                              select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

# Annotate pathway results using KO to KEGG conversion
if(ko2kegg == T){
  daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)
}else{ # already annotated
  daa_annotated_sub_method_results_df <- daa_sub_method_results_df 
  daa_annotated_sub_method_results_df$pathway_name = daa_sub_method_results_df$feature
  daa_annotated_sub_method_results_df$description = daa_sub_method_results_df$feature
  #daa_annotated_sub_method_results_df$description = daa_sub_method_results_df$feature
}


# Plots ----------------------------

# --- PCA all
pca = pathway_pca(abundance = kegg_abundance_all, metadata = metadata, "Dirichlet")

print(pca)

# --- Generate pathway error bar plot
# ... select pathways
#id_sel = grep("metabolism", daa_annotated_sub_method_results_df$feature)
id_sel = c(1,8,10,18,21,24,33,52,54,63,69,82,91,98,100,102,106,109,116,121,149,151,163,180,181,187,190,198,204,209,217,234,235,236,260,268,269,273,283,285,290,293,294)
# Correspond to the following annotations (inspect the following df to subset)
df_sel = data.frame(id_sel,daa_annotated_sub_method_results_df$feature[id_sel])
# iterate after visual inspection to reduce the set
id_sel = c(234,69,235,54,190,283,1,98,294,285,100,217,21)
df_sel = data.frame(id_sel,daa_annotated_sub_method_results_df$feature[id_sel])
names_sel = as.vector(df_sel[,2])

# Please change Group to metadata$your_group_column if you are not using example dataset
p_err_all <- pathway_errorbar(abundance = kegg_abundance_all, daa_results_df = daa_annotated_sub_method_results_df[id_sel,], 
                              Group = metadata$Dirichlet, p_values_threshold = 0.001,
                              #order = "pathway_class",
                              select = NULL, 
                              ko_to_kegg = F, p_value_bar = TRUE, 
                              colors = NULL)#, x_lab = "pathway_name")
print(p_err_all)

# --- The PCA suggests Dirichlet 1 and 10 are far from each other. But previous results
# suggest that 2 vs 7 are more relevant
# .... filter selecting pairs of Dirichlet components to analyse:
dir.pairs.list = list()
dir.pairs.list[[1]] = c(1,2)
dir.pairs.list[[2]] = c(1,3)
dir.pairs.list[[3]] = c(1,6)
dir.pairs.list[[4]] = c(1,7)
dir.pairs.list[[5]] = c(2,3)
dir.pairs.list[[6]] = c(3,6)
dir.pairs.list[[7]] = c(3,7)
dir.pairs.list[[8]] = c(6,7)

colors = c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
           "#6D3F5B", "#1E213D", "#1C542D", "#763C28")

for(i in 1:length(dir.pairs.list)){
  pairs = dir.pairs.list[[i]]
  metadata_red = metadata[metadata$Dirichlet == pairs[1] | metadata$Dirichlet == pairs[2],]
  matched = match(metadata_red$ID.Sample, colnames(kegg_abundance_all))
  kegg_abundance_red = kegg_abundance_all[, matched]
  
  # .... repeat for the reduced set
  daa_results_df_red <- pathway_daa(abundance = kegg_abundance_red, 
                                    metadata = metadata_red, 
                                    group = "Dirichlet", 
                                    daa_method = "ALDEx2", 
                                    include_effect_size = TRUE,
                                    select = NULL, reference = NULL)
  # --- extract most important pathways
  quantile(daa_results_df_red$effect_size)
  hist(daa_results_df_red$effect_size, breaks = 50)
  plot(sort(daa_results_df_red$effect_size))
  hist(daa_results_df_red$log2_fold_change, breaks = 50)
  ##daa_results_df_red$rescale = (daa_results_df_red$log2_fold_change - median(daa_results_df_red$log2_fold_change))#/sd(daa_results_df_red$log2_fold_change)
  #hist(daa_results_df_red$rescale, breaks = 50)
  #daa_results_df_red_filt = daa_results_df_red[abs(daa_results_df_red$rescale)> 0.25,]
  daa_results_df_red_filt = daa_results_df_red # no filtering
  #daa_results_df_red_filt = daa_results_df_red[(abs(daa_results_df_red$log2_fold_change)< 2),] # remove artefacts
  #daa_results_df_red_filt = daa_results_df_red_filt[((abs(daa_results_df_red_filt$log2_fold_change)< 1)&(abs(daa_results_df_red_filt$log2_fold_change)>0.5)),] # remove artefacts
  
  #daa_results_df_red_filt = daa_results_df_red_filt[((daa_results_df_red_filt$effect_size)> 3 | ((daa_results_df_red_filt$effect_size) < -1.5)),]
  #daa_results_df_red_filt = daa_results_df_red_filt[((daa_results_df_red_filt$effect_size) < 0.5),]
  #daa_results_df_red_filt = daa_results_df_red[((daa_results_df_red$effect_size)> 4 | ((daa_results_df_red$effect_size) < -2)),]
  
  daa_sub_method_results_df_red_filt <- daa_results_df_red_filt[daa_results_df_red_filt$method == "ALDEx2_Welch's t test", ]
  
  
  # --- Annotate pathway results using KO to KEGG conversion
  if(ko2kegg == T){
    daa_annotated_sub_method_results_df_red_filt <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df_red_filt, ko_to_kegg = TRUE)
  }else{ # already annotated
    daa_annotated_sub_method_results_df_red_filt <- daa_sub_method_results_df_red_filt
    daa_annotated_sub_method_results_df_red_filt$pathway_name = daa_sub_method_results_df_red_filt$feature
    daa_annotated_sub_method_results_df_red_filt$description = daa_sub_method_results_df_red_filt$feature
    daa_annotated_sub_method_results_df_red_filt$description = daa_sub_method_results_df_red_filt$feature
  }
  
  # filter by type
  #id.all = grep("Metabolism", daa_annotated_sub_method_results_df_red$pathway_class)
  #id.all = c(id.all, grep("Genetic", daa_annotated_sub_method_results_df_red$pathway_class))
  #id.all = c(id.all, grep("Environmental", daa_annotated_sub_method_results_df_red$pathway_class))
  #daa_annotated_sub_method_results_df_red = daa_annotated_sub_method_results_df_red[id.all,]
  
  # re-select pathways. After selecting groups, there may be pathways lost and the id_sel selection will point
  # to the wrong rows, we need relabelling them.
  id_sel_filt = match(names_sel, daa_annotated_sub_method_results_df_red_filt$feature)
  
  thr = 0.0001 # select a threshold to represent pathways
  check = which(daa_annotated_sub_method_results_df_red_filt[id_sel_filt,"p_adjust"] < thr)
  if(length(check) == 0){next} # skip representation
  
  p_err_red <- pathway_errorbar(abundance = kegg_abundance_red, daa_results_df = daa_annotated_sub_method_results_df_red_filt[id_sel_filt,], 
                                Group = metadata_red$Dirichlet,
                                p_values_threshold = thr,
                                legend_title = "Metacommunity",
                                legend_ncol = 2,
                                #order = "pathway_class", 
                                select = NULL, ko_to_kegg = FALSE, p_value_bar = TRUE, 
                                colors = colors[c(pairs[1],pairs[2])])#, x_lab = "pathway_name")
  setwd(dir.root)
  setwd(dir.out)
  file_out = paste0("Plot_ExtErrBar_Met",pairs[1],"vsMet",pairs[2],"_selpaths.pdf")
  pdf(file_out, height = 5, width = 10)
  print(p_err_red)
  dev.off()
}