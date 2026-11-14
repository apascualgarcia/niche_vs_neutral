#####MAKE PHYLOGENETIC TREE OF THE ASVs#####

rm(list=ls()) # clean/reset environment 

here::here()

library(stringi)
library(stringr)
library(phyloseq)
library(ggplot2)
library(ape)
#library(spider) # not used
library(Biostrings)
library(DECIPHER)

# --- Set directories
this.dir = here::here() # project root (works both in RStudio and with Rscript)

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of the configuration file
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the output of the simulation will be saved

# READ IN DATA ------------------------------------------------------------

setwd(dirInput)

###IMPORT AS PHYLOSEQ
mapfile = "Metadata_singleTree_phyloseq.tsv"
otufile = "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv"
repset = "repset_treeholes_dada2_phyloseq.txt"

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile, refseqfilename = repset)
physeq

##REMOVAL OF MITOCHONDRIA AND CHLOROPLAST
Dt <- physeq %>% subset_taxa(Family!= "Mitochondria")
Dt <- Dt %>% subset_taxa(Order!= "Chloroplast")
Dt

#Dt %>%
#  refseq() %>%
#  Biostrings::writeXStringSet("repset_treeholes_dada2_phyloseq_clean.fasta", append=FALSE,
#                              compress=FALSE, compression_level=NA, format="fasta")

##ALIGN ALL SEQUENCES TO THE SAME LENGTH

xx = refseq(Dt)
aligned_repset = subseq(xx, start = 1, end = 252)
aligned_repset_tree = AlignSeqs(aligned_repset)
aligned_repset_tree = as.DNAbin(aligned_repset_tree)
aligned_repset_tree

#1 COMPUTE PHYLOGENETIC DISTANCE MATRIX
dist_matrix_treeholes = dist.dna(aligned_repset_tree)

#2 COMPUTE NEIGHBOUR JOINING TREE
nj_tree_stegen = nj(dist_matrix_treeholes)
nj_tree_stegen

#3 SAVE TREE FILE
write.tree(nj_tree_stegen, file = file.path(this.dir, "3.2_Neutral_Stegen", "nj_tree_stegen.tre"))
