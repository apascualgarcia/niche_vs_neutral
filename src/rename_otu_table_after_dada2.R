rm(list = ls())

here::here()

library(stringi)

# --- Set directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of the configuration file
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the output of the simulation will be saved

# READ IN DATA ------------------------------------------------------------

setwd(dirInput)

##LOAD SEQUENCING FILE NAMES
filenames = list.files(dirInput, pattern = "R1_001.fastq.gz")
filenames

##SPLIT INTO THE OLD AND NEW NAMES
library(stringr)

filenames.split = as.data.frame(do.call("rbind", str_split(filenames, "_", n = 2)))
filenames.split$V2 = gsub(".fastq.gz", "", filenames.split$V2)
head(filenames.split)

##LOAD OTU TABLE AND CHANGE NAMES
otutable = read.csv(paste(dirInput,"otu_table_treeholes_dada2.txt", sep =""), sep = "\t")
head(otutable)

otutable.temp = otutable[,c(2:(ncol(otutable)-1))]
otutable.temp = otutable.temp[,colnames(otutable.temp) %in% filenames.split[,2]]
filenames.split = filenames.split[filenames.split$V2 %in% colnames(otutable.temp),]
otutable.temp = otutable.temp[,match(filenames.split[,2], colnames(otutable.temp))]
colnames(otutable.temp) = gsub("-", ".", filenames.split[,1])
head(otutable.temp)

##REASSEMBLE OTU TABLE WITH NEW NAMES
otutable.renamed = as.data.frame(cbind(otutable$X.ASV_ID, otutable.temp, otutable$taxonomy))
names(otutable.renamed)[names(otutable.renamed) == "otutable$taxonomy"] = "taxonomy"
names(otutable.renamed)[names(otutable.renamed) == "otutable$X.ASV_ID"] = "ASV_ID"
head(otutable.renamed)

##SAVE AS TAB-DELIMITED
write.table(otutable.renamed, paste(dirInput,"otu_table_treeholes_dada2_renamed_tab.txt", sep =""), sep = "\t", row.names = F, quote = F)
