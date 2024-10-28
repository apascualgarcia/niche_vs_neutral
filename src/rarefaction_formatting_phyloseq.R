rm(list=ls()) # clean/reset environment 

here::here()

library(stringi)
library(stringr)
library(phyloseq)
library(ggplot2)


# --- Set directories
this.dir = strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] 

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of the configuration file
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the output will be saved


# READ IN DATA ------------------------------------------------------------

setwd(dirInput)

#sequence table from DADA2 
seqtab_treeholes <- read.table('otu_table_treeholes_dada2_renamed_clean_red.tsv', 
                               sep = "\t", header = T)
rownames(seqtab_treeholes) = seqtab_treeholes$X
head(seqtab_treeholes)[1:10,1:10]

#ASV taxonomy
taxa_wsp <- read.csv('taxonomy.list')

taxa_wsp[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'X')] <- str_split_fixed(taxa_wsp$taxonomy, ';', 7)
head(taxa_wsp)

#Format for phyloseq import
taxa_wsp$'Consensus Lineage' = paste(paste("k__", taxa_wsp$Kingdom, ";"),paste("p__", taxa_wsp$Phylum, ";"), paste("c__", taxa_wsp$Class, ";"),
                                     paste("o__", taxa_wsp$Order, ";"), paste("f__", taxa_wsp$Family, ";"), paste("g__", taxa_wsp$Genus, ";"), 
                                     rep("s__", nrow(taxa_wsp)))
taxa_wsp$`Consensus Lineage` = gsub(" ;", ";", taxa_wsp$`Consensus Lineage`)
taxa_wsp$`Consensus Lineage` = gsub("_ ", "_", taxa_wsp$`Consensus Lineage`)
taxa_wsp_phyloseq = taxa_wsp[,8:9]
rownames(taxa_wsp_phyloseq) = taxa_wsp_phyloseq$X
head(taxa_wsp_phyloseq)

otu_table_treeholes_phyloseq = merge(seqtab_treeholes, taxa_wsp_phyloseq, by = "X")
otu_table_treeholes_phyloseq[1:10,1:10]

names(otu_table_treeholes_phyloseq)[names(otu_table_treeholes_phyloseq) == "X"] = "#OTU ID"

#Read metadata
samdf_treeholes <- read.table('Metadata_singleTree.csv', 
                              sep = "\t", header = T)
samdf_treeholes = subset(samdf_treeholes, select = -c(UID, X))
names(samdf_treeholes)[names(samdf_treeholes) == "ID.OTU"] = "#SampleID"
samdf_treeholes[1:10,1:10]

write.table(samdf_treeholes, "Metadata_singleTree_phyloseq.tsv", sep = "\t", quote = F, row.names = F)
write.table(otu_table_treeholes_phyloseq, "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv", sep = "\t", quote = F, row.names = F)

###IMPORT AS PHYLOSEQ
mapfile = "Metadata_singleTree_phyloseq.tsv"
otufile = "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv"

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)
physeq

##RAREFACTION
#Rarefaction plot
source("Rfunctions.R")
rarefaction.plot = ggrare(physeq, step = 200, se = FALSE) + xlim(c(0,max(sample_sums(physeq)))) + ylim(0,2000) + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y=element_text(size=20), legend.position = "null",
                                                                                                                axis.title.x = element_text(size = 20, margin = margin(r = 15)), axis.title.y = element_text(size =20, margin = margin(r = 10))) + xlab("Number of reads")

tiff("Rarefaction_plot_treeholes.tiff", units="cm", width=15, height=15, res=300)
rarefaction.plot
dev.off()

set.seed(123)
Dt = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                       rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
Dt

##REMOVAL OF MITOCHONDRIA AND CHLOROPLAST
Dt <- Dt %>% subset_taxa(Family!= "Mitochondria")
Dt <- Dt %>% subset_taxa(Order!= "Chloroplast")
Dt

# CONTINUE WITH RAREFIED READ COUNTS ***************** 