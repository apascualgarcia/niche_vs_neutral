###################################
#  split_otu_table_Dirichlet.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-03-21
# Script Name: split_otu_table_Dirichlet.R 
# Script Description: This script takes the
#   ASV table (after filtering rare ASVs and low
#   abundant samples) and splits it according to
#   the output given by the optimal division of 
#   data into Dirichlet components given by the
#   method by Holmes et al. It also returns transposed
#   and rarefied matrices, including the original
#

rm(list = ls())
# SET WORKING DIRECTORY -----------------------------
this.dir = here::here() # project root (works both in RStudio and with Rscript)
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of input data 
dirMeta=paste(this.dir,"/1_samples_properties/",sep="") # Dir of dada2 output
dirDirich=paste(this.dir,"/2_Dirichlet/",sep="") # Dir of Dirichlet analysis

args = commandArgs(trailingOnly = TRUE)
Ncomp = if (length(args) > 0) as.integer(args[1]) else 10  # Number of Dirichlet components (Rscript split_otu_table_Dirichlet.R 7)
dirOutput=paste(this.dir,"/2_Dirichlet/otu_tables_K",Ncomp,"/",sep="") # Dir of output data
dir.create(dirOutput, showWarnings = FALSE)
file.otu = "otu_table_treeholes_dada2_renamed_clean_red.tsv"
file.otu.label = "otu_table_treeholes_dada2_renamed_clean_red" # file name without extension
file.meta = "Metadata_singleTree.tsv"
file.meta.label = "Metadata_singleTree" # file name without extension
file.dirich = paste0("otu_table_treeholes_dada2_Dirichlet_K",Ncomp,".z")

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("stringr","phyloseq") # list of packages to load
n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}


# READ FILES  ------------------------------------

setwd(dirInput)
otu.table = read.table(file = file.otu, header = T, sep = "\t", row.names = 1)
head(otu.table)[1:5,1:5]
otu.pseq = otu_table(as.matrix(otu.table), taxa_are_rows = TRUE)


setwd(dirMeta)
df.meta = read.table(file = file.meta, header = T, sep = "\t")
head(df.meta)[1:5,1:5]

setwd(dirDirich)
df.dirich = read.table(file = file.dirich, sep = ",",skip = 1, row.names = 1)
colnames(df.dirich) = c(unlist(strsplit(paste0("C",seq(1,Ncomp),collapse = "_"),
                                        split = "_", fixed = T)))
head(df.dirich)[1:5,1:Ncomp]

# SPLIT DATA ---------------

matched = match(rownames(df.dirich), df.meta$ID.Sample) # look for those samples classified by a Dirichlet component
df.meta$Dirichlet = 0
head(df.meta)[1:5,]
df.meta$Dirichlet[matched] = apply(df.dirich, MARGIN = 1, # identify the component, assign an integer
                          FUN = function(x){which(x > 0)})
head(df.meta)[1:5,]
df.meta$Dirichlet[df.meta$Dirichlet == 0] = NA  # identify the remainder as NA

#comp = colnames(df.dirich)[1] # debug
setwd(dirOutput)
for(comp in colnames(df.dirich)){
  cat(" -------------  \n")
  cat(paste0("Extracting component ",comp,"\n"))
  id.comp = which(df.dirich[, comp] > 0)
  samples = rownames(df.dirich[id.comp, ])
  matched = match(samples, colnames(otu.table))
  otu.table.tmp = as.matrix(otu.table[, matched, drop = FALSE]) # drop = FALSE keeps the ASV names when a component has a single sample
  otu.table.tmp = otu.table.tmp[rowSums(otu.table.tmp)> 0, , drop = FALSE]
  otu.table.tmp.t = t(otu.table.tmp)
  otu.pseq.tmp = otu_table(as.matrix(otu.table.tmp), taxa_are_rows = TRUE)
  set.seed(1) # Seed 1 is the one used for the committed results: rngseed = TRUE in phyloseq meant set.seed(1)
  otu.pseq.tmp.rar =  rarefy_even_depth(otu.pseq.tmp, 
                                         sample.size = min(sample_sums(otu.pseq.tmp)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  cat(paste0("Dimensions rarefied table ",dim(otu.pseq.tmp.rar@.Data)))
  cat(paste0("\n Number ASVs with zero counts: ",length(which(rowSums(otu.pseq.tmp.rar@.Data) == 0))),"\n")
  
  file.out = paste0(file.otu.label,"_Dirich",comp,".tsv")
  write.table(otu.table.tmp, file = file.out, quote = F, sep = "\t")
  file.out = paste0(file.otu.label,"_Dirich",comp,"_t.tsv")
  write.table(otu.table.tmp.t, file = file.out, quote = F, sep = "\t")
  file.out = paste0(file.otu.label,"_Dirich",comp,"_rar.csv")
  # NMGS input: first column must be labelled "OTUs"
  write.table(data.frame(OTUs = rownames(otu.pseq.tmp.rar@.Data), otu.pseq.tmp.rar@.Data, check.names = F),
              file = file.out, quote = F, sep = ",", row.names = F)
}

# --- Rarefy original table
set.seed(1)
otu.pseq.rar =  rarefy_even_depth(otu.pseq, 
                                      sample.size = min(sample_sums(otu.pseq)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
setwd(dirInput)
file.out = paste0(file.otu.label,"_rar.csv")
# NMGS input: first column must be labelled "OTUs"
write.table(data.frame(OTUs = rownames(otu.pseq.rar@.Data), otu.pseq.rar@.Data, check.names = F),
            file = file.out, quote = F, sep = ",", row.names = F)

# --- Write new metadata table
setwd(dirMeta)
file.out = paste0(file.meta.label,"_wDirich-K",Ncomp,".tsv")
write.table(df.meta,file = file.out, quote = F, sep = "\t", row.names = F)
