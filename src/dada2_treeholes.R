rm(list = ls())

here::here()

library(stringi)

# --- Set directories
this.dir = here::here() # project root (works both in RStudio and with Rscript)

dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirInput=paste(this.dir,"/1_dada2/",sep="") # Dir of the configuration file
dirOutput=paste(this.dir,"/1_dada2/",sep="") # Dir where the output of the simulation will be saved

# READ IN DATA ------------------------------------------------------------

setwd(dirInput)

#' Load DADA2 and required packages

library(dada2); packageVersion("dada2") # the dada2 pipeline
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(dplyr); packageVersion("dplyr") # for manipulating data
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots


# Set up pathway to idemp (demultiplexing tool) and test
# ---- EDIT: external tools and reference (not in the repo) ----
# idemp (github.com/yhwu/idemp) demultiplexes the raw reads with a barcode table:
#   idemp -b barcodes.txt -I1 <I1.fastq.gz> -R1 <R1.fastq.gz> -R2 <R2.fastq.gz> -o 1_dada2/
# which writes one Undetermined_S0_L001_R{1,2}_001.fastq.gz_<sample> file per sample (renamed below).
idemp <- "idemp" # CHANGE ME: full path to the idemp binary if it is not in your PATH
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "cutadapt" # CHANGE ME: full path to cutadapt if it is not in your PATH
# dada2-formatted SILVA training set used by assignTaxonomy() (version not recorded; the
# taxonomy names, e.g. Bacteroidota, imply SILVA >= 138; https://zenodo.org/records/4587955)
silva.ref <- file.path(dirInput, "silva_nr99_v138.1_train_set.fa.gz") # CHANGE ME
system2(cutadapt, args = "--version") # Check by running shell command from R

# Set path to shared data folder and contents
data.fp <- dirInput

# List all files in shared folder to check path
list.files(data.fp)

# Set file paths for fastqs
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

project.fp <- dirInput 

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(project.fp, "01_preprocess")
demultiplex.fp <- file.path(preprocess.fp, "demultiplexed_2")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 


# Rename files - use gsub to get names in order!
R1_names <- gsub(paste0(project.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                 list.files(project.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(project.fp, pattern="R1", full.names = TRUE), 
            paste0(project.fp, "/R1_", R1_names))

R2_names <- gsub(paste0(project.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                 list.files(project.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(project.fp, pattern="R2", full.names = TRUE),
            paste0(project.fp, "/R2_", R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(project.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(project.fp, pattern="R2_", full.names = TRUE))

# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

xx = gsub("R2", "R1", fnRs)
fnFs = fnFs[fnFs %in% xx]
length(fnFs)
length(fnRs)

yy = gsub("R1", "R2", fnFs)
fnRs = fnRs[fnRs %in% yy]
length(fnFs)
length(fnRs)

dd = gsub("R2", "R1", fnRs.filtN)
fnFs.filtN = fnFs.filtN[fnFs.filtN %in% dd]
length(fnFs.filtN)
length(fnRs.filtN)

zz = gsub("R1", "R2", fnFs.filtN)
fnRs.filtN = fnRs.filtN[fnRs.filtN %in% zz]
length(fnFs.filtN)
length(fnRs.filtN)

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = T) 
# CHANGE multithread to FALSE on Windows (here and elsewhere in the program)


# Set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME # this is 515f
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME # this is 806Br

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# Create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

# Run Cutadapt
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
## should all be zero!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[4]]))



# Put filtered reads into separate sub-directories for big data workflow
dir.create(filter.fp)
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 
dir.create(subF.fp)
dir.create(subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs))
file.rename(from = fnFs.cut, to = fnFs.Q)
file.rename(from = fnRs.cut, to = fnRs.Q)

# File parsing; create file names and make sure that forward and reverse files match
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered") # ...
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#CHECK FOR QUALITY
if( length(fastqFs) <= 20) {
  plotQualityProfile(paste0(subF.fp, "/", fastqFs))
  plotQualityProfile(paste0(subR.fp, "/", fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(fastqFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs[rand_samples]))
}

ggsave(plot = fwd_qual_plots, filename = paste0(filter.fp, "/fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots, filename = paste0(filter.fp, "/rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")


filt_out <- filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
                          rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                          truncLen=c(120,220), maxEE=c(2,2), truncQ=2, maxN=0, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=TRUE)

# look at how many reads were kept
head(filt_out)


# summary of samples in filt_out by percentage
filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))


# figure out which samples, if any, have been filtered out
remaining_samplesF <-  fastqFs[rand_samples][
  which(fastqFs[rand_samples] %in% list.files(filtpathF))] # keep only samples that haven't been filtered out
remaining_samplesR <-  fastqRs[rand_samples][
  which(fastqRs[rand_samples] %in% list.files(filtpathR))] # keep only samples that haven't been filtered out
fwd_qual_plots_filt <- plotQualityProfile(paste0(filtpathF, "/", remaining_samplesF))
rev_qual_plots_filt <- plotQualityProfile(paste0(filtpathR, "/", remaining_samplesR))


ggsave(plot = fwd_qual_plots_filt, filename = paste0(filter.fp, "/fwd_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots_filt, filename = paste0(filter.fp, "/rev_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")


# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- substring(basename(filtFs), regexpr("_", basename(filtFs)) + 1) # doesn't drop fastq.gz
sample.names <- gsub(".fastq.gz", "", sample.names)
sample.namesR <- substring(basename(filtRs), regexpr("_", basename(filtRs)) + 1) # doesn't drop fastq.gz
sample.namesR <- gsub(".fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names


set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates (Notes: randomize default is FALSE)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)

errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

#ggsave(plot = errF_plot, filename = paste0(filter.fp, "/errF_plot.png"), 
#       width = 10, height = 10, dpi = "retina")
#ggsave(plot = errR_plot, filename = paste0(filter.fp, "/errR_plot.png"), 
#       width = 10, height = 10, dpi = "retina")


# Dereplicate forward reads
derepF.p <- derepFastq(filtFs)
names(derepF.p) <- sample.names
# Infer sequences for forward reads
dadaF.p <- dada(derepF.p, err = errF, multithread = TRUE, pool = TRUE)
names(dadaF.p) <- sample.names

# Dereplicate reverse reads
derepR.p <- derepFastq(filtRs)
names(derepR.p) <- sample.names
# Infer sequences for reverse reads
dadaR.p <- dada(derepR.p, err = errR, multithread = TRUE, pool = TRUE)
names(dadaR.p) <- sample.names

# Merge reads together
mergers <- mergePairs(dadaF.p, derepF.p, dadaR.p, derepR.p)

seqtab <- makeSequenceTable(mergers)

# Save table as an r data object file
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds"))


# Read in RDS 
st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Print percentage of our seqences that were not chimeric.
100*sum(seqtab.nochim)/sum(seqtab)
## [1] 98.33841

# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, silva.ref, tryRC = TRUE,
                      multithread=TRUE)

# Write results to disk
saveRDS(seqtab.nochim, paste0(table.fp, "/seqtab_final.rds"))
saveRDS(tax, paste0(table.fp, "/tax_final.rds"))


# Flip table
seqtab.t <- as.data.frame(t(seqtab.nochim))

# Pull out ASV repset
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- dplyr::mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
taxonomy <- as.data.frame(tax)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
taxonomy_for_mctoolsr <- unite_(taxonomy, "taxonomy", 
                                c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"),
                                sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the taxonomy dataframe for the writeRepSetFasta function
taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  dplyr::select(ASV, name) %>%
  dplyr::rename(seq = ASV)

# write fasta file
writeRepSetFasta(taxonomy_for_fasta, paste0(table.fp, "/repset.fasta"))

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ASV <- NULL 

# Set name of table in mctoolsr format and save
out_fp <- paste0(table.fp, "/seqtab_wTax_mctoolsr.txt")
names(seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)
