###################################
#     filter_rare_ASVs.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c) Alberto Pascual-García, 2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-03-21
# Script Name: filter_rare_ASVs.R  
# Script Description: This script simply takes the ASVs
# table and it filters any rare ASVs and samples with
# low number of reads.
# Notes:
#
#
rm(list = ls())

# SET WORKING DIRECTORY -----------------------------
this.dir = here::here() # project root (works both in RStudio and with Rscript)

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of dada2 output 
dirMeta=paste(this.dir,"/1_samples_properties/",sep="") # Dir of dada2 output
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the filtered matrix will be saved


# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("tidyverse", "stringr", "ggplot2") # list of packages to load
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


# READ IN DATA ------------------------------------------------------------

setwd(dirInput)
#-- sequence table from DADA2
seqtab_treeholes<-read.table('otu_table_treeholes_dada2_renamed_clean.txt', 
                           sep = "\t", header = T)
rownames(seqtab_treeholes) = seqtab_treeholes$ASV_ID
seqtab_treeholes = subset(seqtab_treeholes, select = -c(ASV_ID))
head(seqtab_treeholes)[1:10,1:10]
seqtab_treeholes = t(seqtab_treeholes)
dim(seqtab_treeholes)              

#-- Taxonomy
taxa_wsp<-read.csv('taxonomy.list')

#-- Metadata
setwd(dirMeta)
samdf_treeholes<-read.table('Metadata_singleTree.tsv', 
                          sep = "\t", header = T)
head(samdf_treeholes)[1:10,1:10]

taxa_wsp<-read.csv(file.path(dirInput,'taxonomy.list'))

nrow(seqtab_treeholes)
ncol(seqtab_treeholes)
nrow(samdf_treeholes)
quantile(rowSums(seqtab_treeholes),q = c(0,0.25,0.5,0.75,1)) 
plot(sort(rowSums(seqtab_treeholes))) # there is a change in the slope for samples < 5000 and those > 5000
# 0%   25%   50%   75%  100% 
# 0  5907  8275 11046 24314 

# FILTER SEQUENCES --------------------------------------------------------

#remove communities with less than 5K sequences (all ASV columns are counted)
seqtab_treeholes<-seqtab_treeholes[which(rowSums(seqtab_treeholes)>5000),]
nrow(seqtab_treeholes)

#remove OTUs with less than 100 reads across all samples
seqtab_treeholes<-seqtab_treeholes[,colSums(seqtab_treeholes)>=100] # (-which() would drop all columns if none matched)
ncol(seqtab_treeholes)

#remove OTUs occurring in less than 10 samples
seqtab_treeholes<-seqtab_treeholes[,colSums(seqtab_treeholes>0)>=10]
ncol(seqtab_treeholes)

#remove CONTROL communities
seqtab_treeholes <- seqtab_treeholes[!grepl("CONTROL", rownames(seqtab_treeholes)),]

# FILTER THE SAMPLE METADATA AND TAXA TABLE ACCORDING TO REMAINING COMMUNITIES AND ASVS, respectively -------------------------------------------------------------------------

#list communities in samdf which are still in seqtab_treeholes after filtering
samples_to_keep<-which(samdf_treeholes$filename%in%rownames(seqtab_treeholes))
#filter them out
samdf_treeholes<-samdf_treeholes[samples_to_keep,]

#list sequences/ASVs in samdf which are still in seqtab_treeholes after filtering
sequences_to_keep<-which(taxa_wsp$X%in%colnames(seqtab_treeholes))
#filter them out
taxa_wsp<-taxa_wsp[sequences_to_keep,]


# Write outputs ---------------
# ------- Otu tables
setwd(dirInput)
seqtab_treeholes_t = t(seqtab_treeholes)
write.table(seqtab_treeholes, file = "otu_table_treeholes_dada2_renamed_clean_red_t.tsv",
            sep = "\t", col.names = NA, row.names = T, quote = F)
# The first column is labelled "Taxon": with an empty label the Dirichlet
# multinomial mixture code (DirichletMixtureGHPFit, input = the .csv) silently
# drops one sample. R scripts reading these tables rename it as needed.
seqtab_treeholes_t = data.frame(Taxon = rownames(seqtab_treeholes_t), seqtab_treeholes_t, check.names = F)
write.table(seqtab_treeholes_t, file = "otu_table_treeholes_dada2_renamed_clean_red.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(seqtab_treeholes_t, file = "otu_table_treeholes_dada2_renamed_clean_red.csv",
            sep = ",", row.names = F, quote = F)

# ------ Metadata and taxa


# CONTINUE HERE ***************** 22fEB24

