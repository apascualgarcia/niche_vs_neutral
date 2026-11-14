####FORMAT REPRESENTATIVE 16S SEQUENCE FILE FOR PHYLOSEQ####

rm(list = ls())

here::here()

library(stringi)
library(stringr)
library(Biostrings)
library(phyloseq) # import_qiime, refseq, subset_taxa

# --- Set directories
this.dir = here::here() # project root (works both in RStudio and with Rscript)

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of the configuration file
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the output of the simulation will be saved

# READ IN DATA ------------------------------------------------------------

setwd(dirInput)

repset_fasta = read.table("repset_treeholes_dada2.fasta", sep = "\t")
head(repset_fasta)

repset_fasta[c("ASV", "Taxonomy")] <- str_split_fixed(repset_fasta$V1, ' ', 2)
repset_fasta = subset(repset_fasta, select = -V1)

write.table(repset_fasta$ASV, file = "repset_treeholes_dada2_phyloseq.txt",
            sep = "\t", col.names = F, row.names = F, quote = F)


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

Dt %>%
  refseq() %>%
  Biostrings::writeXStringSet("repset_treeholes_dada2_phyloseq_clean.fasta", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")


