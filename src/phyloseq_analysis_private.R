##################################################
# phyloseq_analysis.R
##################################################
# In this script I perform the typical pipeline of phyloseq analysis, to
# study phylogenetic distributions of taxa in the different samples, alpha diversity
# quantification, and dimensionality reduction.
# This is an ongoing modification of Cordero's analysis. I haven't done anything yet, 
# I just introduced the path of the data.
###################################################
# Madrid, Oct 18th, 2023
# Spanish National Centre for Biotechnology
# apascualgarcia.github.io
###################################################


library(phyloseq)
library(reshape2)
#library(RDPutils) # Notused
library(ggplot2)

# START EDITING --------
dir.otu = "../1_dada2" # set directories relative to src dir
dir.samples = "../1_samples_properties"
fileOTU = "otu_table_treeholes_dada2_clean.txt"
fileSample = "Metadata_AllSamples_4cytos.tsv"
fileTaxonomy = "taxonomy.list" 

# --- Minimum number of reads to keep a sample/taxa
Nprune.sample = 1000
Nprune.taxa = 10
# STOP EDITING --------

# Read files ----------
# --- Set environment 
current_path = rstudioapi::getActiveDocumentContext()$path 
dir.src = dirname((current_path))
setwd(dir.src)

# --- Load OTUs
setwd(dir.otu)
otu.in=read.csv(fileOTU, sep = "\t", row.names = 1)
#colSums(otu.in)
head(otu.in)[1:5,1:5]
otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)
setwd(dir.src)

# --- Load samples metadata
setwd(dir.samples)
sample_metadata = import_qiime_sample_data(fileSample)
setwd(dir.src)

# --- Load taxonomy and clean it 
setwd(dir.otu)
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("Kingdom","Phylum",
                                 "Class","Order","Family","Genus","taxa_id")
rownames(taxa.in)=taxa.in$taxa_id
taxa.in=subset(taxa.in,select=-c(taxa_id))
tax.pseq = tax_table(as.matrix(taxa.in))
setwd(dir.src)

# Preprocessing (controls) ------------------
# --- Control that there are no empty rows/columns 
which(taxa_sums(otu.pseq) == 0) # 0
otu.pseq = prune_taxa(taxa_sums(otu.pseq) > 0, otu.pseq) # delete these guys
any(taxa_sums(otu.pseq) == 0) # double check (FALSE)
which(sample_sums(otu.pseq) == 0)
otu.pseq = prune_samples(sample_sums(otu.pseq) > 0, otu.pseq)
any(sample_sums(otu.pseq) == 0) # double check (FALSE)

# --- Check the distribution of samples and taxa
readsumsdf = data.frame(nreads = sort(taxa_sums(otu.pseq), TRUE), 
                        sorted = 1:ntaxa(otu.pseq),type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(otu.pseq), TRUE), 
                                          sorted = 1:nsamples(otu.pseq), type = "Samples"))

# ....... plot distribution
setwd(dir.samples)
file.out = "Plot_Number_reads.pdf"
pdf(file = file.out, width = 10)
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p = p + ggtitle(title) + scale_y_log10() + 
  facet_wrap(~type, 1, scales = "free") + 
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 25)) 
print(p)
dev.off()
min(sample_sums(otu.pseq)) # The minimum sample sum is 1

# --- Remove samples with less than 1K reads and taxa with less than 10
otu.pseq = prune_taxa(taxa_sums(otu.pseq) > Nprune.taxa, otu.pseq) # delete these guys
otu.pseq = prune_samples(sample_sums(otu.pseq) > Nprune.sample, otu.pseq)
min(sample_sums(otu.pseq))
min(taxa_sums(otu.pseq))

# Finally, create a single phyloseq object
treeholes = merge_phyloseq(otu.pseq, sample_metadata, tax.pseq) 

# Biodiversity analysis ----------------
# Now we want to quantify biodiversity more systematically
# --- Rarefy
# Quick check to see how much things change
set.seed(30052018) # Today's date 30/05/2018. Stored for reproducibility

# Select the analysis to be done, WARNING: after the final choice of rarefaction level you should repeat this step
# now skip samples that you may not want to study
particles.enke=subset_samples(particles, Substrate != "AgaroseChitosan") # Exclude chitosan to work on Tim Enke paper
particles.enke=subset_samples(particles.enke, Substrate != "Chitosan") # Exclude chitosan
particles.enke=subset_samples(particles.enke, Media != "Seawater") 
particles.nosea = subset_samples(particles, Media != "Seawater") # exclude seawater from the analysis
particles.nomix=subset_samples(particles, Substrate != "AgaroseChitosan") # Exclude mixed substrates
particles.nomix=subset_samples(particles.nomix, Substrate != "AgaroseAlginate") # Exclude mixed substrates
particles.nomix=subset_samples(particles.nomix, Substrate != "AgaroseCarageenan") # Exclude mixed substrates
particles.nomix.nosea=subset_samples(particles.nomix, Media != "Seawater") 
particles.nomix.nobeads=subset_samples(particles.nomix, Media != "Beads") 

# and select the final set
particles.select=particles.nomix.nosea # Determine the set for subsequent analysis
ExtraLabel=" Beads; "#; Media = Seawater" # to add in the title, double check the code below of the plot 
Merged="RepMerged" # should the replicas be merged? one of "RepMerged" or "NotMerged"
Analysis=paste("Substrate-Beads",Merged,sep="-") # To add in the name of the files

# Rarefaction, determine levels
seqSize=c(1000,2500) #c(50,100,500,1000,2500,5000,7500,10000,25000,50000)

dir="/home/apascual/Nextcloud/Research/Projects/Cordero/16sRNA/Biodiv_vs_raref"
setwd(dir)
substrate_list=levels(particles.select@sam_data$Substrate)
for(i in 1:length(seqSize)){
  sampling=seqSize[i]
  particles.select.rar = rarefy_even_depth(particles.select, sample.size = sampling)
  if(Merged=="RepMerged"){ # Merge the replicas
    # First prepare a new field combining Media, Substrates and Time, we will merge samples using this factor
    sample_data(particles.select.rar)$MediaSubTime <- as.factor(paste0(sample_data(particles.select.rar)$Media, 
                                                                       sample_data(particles.select.rar)$Substrate,sample_data(particles.select.rar)$Time))
    particles.select.rar.RepMerged =merge_samples(particles.select.rar, "MediaSubTime")
    # Repair the levels used in the new data
    sample_data(particles.select.rar.RepMerged)$Substrate <- as.factor(levels(sample_data(particles.select.rar)$Substrate)[get_variable(particles.select.rar.RepMerged, 
                                                                                                                                        "Substrate")]) # Repair levels
    sample_data(particles.select.rar.RepMerged)$Media <- as.factor(levels(sample_data(particles.select.rar)$Media)[get_variable(particles.select.rar.RepMerged, 
                                                                                                                                "Media")]) # Repair levels
    sample_data(particles.select.rar.RepMerged)$Time <- as.factor(levels(sample_data(particles.select.rar)$Time)[get_variable(particles.select.rar.RepMerged, 
                                                                                                                              "Time")])
    particles.select.rar=particles.select.rar.RepMerged
    idx=sort(as.numeric(levels(particles.select.rar@sam_data$Time)),index.return=TRUE)$ix # reorder levels
    particles.select.rar@sam_data$Time=factor(particles.select.rar@sam_data$Time,
                                              levels(particles.select.rar@sam_data$Time)[idx]) # reorder level
  }
  richness.dd=estimate_richness(particles.select.rar)
  titleOut=paste("Table_Richness_Raref",sampling,".tsv",sep="") # For the output files
  write.table(richness.dd,file=titleOut)
  titlePlot=paste("Plot_RichnessBy",Analysis,"_Raref",sampling,".pdf",sep="") # For the plot title
  pdf(file=titlePlot,width=15)
  print(plot_richness(particles.select.rar,x="Substrate",
                      measures=c("Observed", "Chao1", "Shannon"),color="Time",
                      title=paste("Samples size = ",sampling,ExtraLabel,sep="")))
  dev.off()
  
  for(j in 1:length(substrate_list)){
    Sub=substrate_list[j]
    titlePlot2=paste("Plot_RichnessBy",Analysis,"_",Sub,"_Raref",sampling,".pdf",sep="")
    particles.select.rar.sub = subset_samples(particles.select.rar, Substrate == Sub) # levels are already ordered
    # but now I don't want to use plot_richness cause I need full control of the plot so I have to extract the data and reorder the dataframe (not the levels) manually
    time.factor=particles.select.rar.sub@sam_data$Time
    time.numeric=as.numeric(levels(time.factor))[time.factor]
    vecColor=time.factor
    #idx=sort(as.numeric(levels(particles.select.rar.sub@sam_data$Time)),index.return=TRUE)$ix 
    #vecColor=factor(particles.select.rar.sub@sam_data$Time,
    #                                         levels(particles.select.rar.sub@sam_data$Time)[idx]) # reorder levels
    #particles.select.rar.sub@sam_data$Time=factor(particles.select.rar.sub@sam_data$Time,
    #                                              levels(particles.select.rar.sub@sam_data$Time)[idx]) # reorder levels
    #particles.select.rar.sub@sam_data$Time=as.numeric(levels(particles.select.rar.sub@sam_data$Time)[idx])
    richness.dd=estimate_richness(particles.select.rar.sub) # I am gonna build manually a df to use geom_smooth
    df.tmp=as.data.frame(cbind(as.matrix(time.numeric),richness.dd$Shannon))
    idx=sort(df.tmp$V1,index.return=TRUE)$ix
    df.tmp=df.tmp[idx,]
    vecColor=vecColor[idx]
    colnames(df.tmp)=c("Time","Shannon")
    # To increase control add here some aesthetics, note I keep only Shannon diversity
    #p=plot_richness(particles.select.rar.sub,x="Time",
    #measures=c("Shannon"),color="Time", #"Time",vecColor,
    #title=paste("Substrate = ",Sub,ExtraLabel,sep="")) # "; Samples size = ",sampling
    p=ggplot(df.tmp,aes(x=Time,y=Shannon))+ggtitle(paste(ExtraLabel,Sub,sep="")) #,color=vecColor))
    p=p+geom_point(aes(color=vecColor), size=15,alpha=0.7)+labs(y = "Shannon diversity")
    p=p+geom_smooth(method="lm",show.legend = FALSE,se=FALSE) # does not work if data are factors
    p=p+ theme(plot.title=element_text(hjust = 0.5,size=30), axis.text = element_text(size=40),axis.title = element_text(size=30),
               legend.position = "none",strip.text.x=element_blank(),strip.background = element_blank()) # Remove strip text and legend 
    #p
    pdf(file=titlePlot2,width=10)
    # print(plot_richness(particles.select.rar.sub,x="Time",
    #                     measures=c("Observed", "Chao1", "Shannon"),color="Replica",
    #                     title=paste("Substrate = ",Sub,"; Samples size = ",sampling,sep=""))) # With replicas (old code)
    print(p)
    dev.off()
  }
}

# Representations for accepted rarefied level -------------------------
# --- Set an output directory
dir="/home/apascual/Nextcloud/Research/Projects/Cordero/16sRNA/Biodiv_vs_raref"
setwd(dir)
# .... Create a palette for the plots
library(RColorBrewer) #
qual_col_pals = brewer.pal.info[(brewer.pal.info$category == 'qual'),] # Get qualitative palettes
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[1]=col_vector[8] # I rearrange the colors to have grey for Others
col_vector[8]=col_vector[60] # and then I assigned to the two guys with grey positions
col_vector[16]=col_vector[61] # other colors
col_vector[6]=col_vector[59] # And to one of the pink ones, to much pink around

# --- Create some plots for the minimum acceptable level of rarefaction
#otu.pseq.rar = rarefy_even_depth(otu.pseq, sample.size = 1000) # 1018 is the minimum sampling sites, and the alpha diversity pattern is already there
#particles.rar=merge_phyloseq(otu.pseq.rar,sample_metadata,tax.pseq) 
sampling=1000 # select the size here, 1018 is the minimum sampling sites, and the alpha diversity pattern is already there
particles.rar = rarefy_even_depth(particles.select, sample.size = sampling)
# **** Warning! Note that in the following results depend on the choice of set made above, and may be even some lines
#               that will not work for some choices

# ..... Repeat the plot to see if the OTU distro significantly changes
otu.pseq.rar=particles.rar@otu_table

readsumsdf = data.frame(nreads = sort(taxa_sums(otu.pseq.rar), TRUE), 
                        sorted = 1:ntaxa(otu.pseq.rar),type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(otu.pseq.rar), TRUE), 
                                          sorted = 1:nsamples(otu.pseq.rar), type = "Samples"))
title = "Total number of reads - rarefied"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# --- Identify the top 19 genus and families for each substrate (the 20th will be the remainder)
substrates=levels(particles.rar@sam_data$Substrate)
topESVs=matrix()
particles.rar.bySub.list=list() # Separate the samples by substrate
#particles.rar.sea=prune_samples(particles.rar@sam_data$Media=="Seawater",particles.rar)
#particles.rar.beads=prune_samples(particles.rar@sam_data$Media=="Beads",particles.rar)
particles.tmp=particles.rar #.beads # First look by substrate for beads
for(i in 1:length(substrates)){ # for each substrate
  Substrate=substrates[i]
  particles.rar.bySub.list[[i]]=prune_samples(particles.tmp@sam_data$Substrate==Substrate,particles.tmp) # get the subset of samples
  topESVs = c(topESVs,names(sort(taxa_sums(particles.rar.bySub.list[[i]]), TRUE)[1:19])) # look for the 19 most abundant species for that sample and add it to the previous
}

# ..... Now look for the common set 
topESVs=unique(topESVs) # clean up those species that are repeated
topESVs=topESVs[!is.na(topESVs)]
taxTop = cbind(tax_table(particles.rar), family = NA, genus=NA)
taxTop[topESVs, "family"] <- as(tax_table(particles.rar)[topESVs, "Family"], 
                                "character")
taxTop[topESVs, "genus"] <- as(tax_table(particles.rar)[topESVs, "Genus"], 
                               "character")
taxTop[which(is.na(taxTop))]="_Other_"
# ..... Add the genus to the object
tax_table(particles.rar) <- tax_table(taxTop)

# --- Merge the samples by substrate
# First prepare a new field combining Media and substrates, we will merge samples using this factor
sample_data(particles.rar)$MediaSub <- as.factor(paste0(sample_data(particles.rar)$Media, sample_data(particles.rar)$Substrate))
particles.rar.bySub =merge_samples(particles.rar, "MediaSub")
sample_data(particles.rar.bySub)$MediaSub <- levels(sample_data(particles.rar)$MediaSub)[get_variable(particles.rar.bySub, 
                                                                                                      "MediaSub")]
sample_data(particles.rar.bySub)$Substrate <- levels(sample_data(particles.rar)$Substrate)[get_variable(particles.rar.bySub, 
                                                                                                        "Substrate")] # Reorder levels
sample_data(particles.rar.bySub)$Media <- levels(sample_data(particles.rar)$Media)[get_variable(particles.rar.bySub, 
                                                                                                "Media")] # Reorder levels
#sample_data(particles.rar.bySub)$Substrate <- levels(sample_data(particles.rar)$Substrate) # Add the names back
particles.rar.bySub = transform_sample_counts(particles.rar.bySub, function(x) 100 * x/sum(x))

# .... abundances per substrate
#particles.rar.bySub.19 = prune_taxa(topESVs, particles.rar.bySub) # if you want to skip the non-top species
particles.rar.bySub.19 = particles.rar.bySub  # If you want to keep the remainder in grey (only those in the tax_table list different than NA will be coloured alredy)

# .... First plot
plotTitle=paste("BarPlot_diversity-",Analysis,"_ByGenus.pdf",sep="")
pdf(file=plotTitle,width=20,height = 10)
title = " "
#plot_bar(particles.rar.bySub.19,"Substrate",fill="genus",title=NULL)+
#  coord_flip()+ylab("Percentage of Sequences")+ scale_fill_manual(values = col_vector)
plot_bar(particles.rar.bySub.19,x ="MediaSub",fill="genus",title=NULL)+
  coord_flip()+ylab("Percentage of Sequences")+xlab("Media/Substrate")+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4)+
  scale_fill_manual(values = col_vector)+theme(axis.ticks.x = )
dev.off()

# --- Now, we do the same representing also the time
# ..... first, prepare a dummy variable
sample_data(particles.rar)$MediaSubTime <- paste0(sample_data(particles.rar)$MediaSub, sample_data(particles.rar)$Time)
sample_data(particles.rar)$SubTime <- paste0(sample_data(particles.rar)$Substrate, sample_data(particles.rar)$Time)
particles.rar.bySubTime=merge_samples(particles.rar, "SubTime") # And proceed as above, first merge
particles.rar.byMediaSubTime=merge_samples(particles.rar, "MediaSubTime") # And proceed as above, first merge

# ..... repair factors
sample_data(particles.rar.bySubTime)$Substrate <- as.factor(levels(sample_data(particles.rar)$Substrate)[get_variable(particles.rar.bySubTime, 
                                                                                                                      "Substrate")])
sample_data(particles.rar.bySubTime)$Time <- as.factor(levels(sample_data(particles.rar)$Time)[get_variable(particles.rar.bySubTime, 
                                                                                                            "Time")])
sample_data(particles.rar.byMediaSubTime)$MediaSub <- as.factor(levels(sample_data(particles.rar)$MediaSub)[get_variable(particles.rar.byMediaSubTime, 
                                                                                                                         "MediaSub")])
sample_data(particles.rar.byMediaSubTime)$Substrate <- as.factor(levels(sample_data(particles.rar)$Substrate)[get_variable(particles.rar.byMediaSubTime, 
                                                                                                                           "Substrate")])
sample_data(particles.rar.byMediaSubTime)$Media <- as.factor(levels(sample_data(particles.rar)$Media)[get_variable(particles.rar.byMediaSubTime, 
                                                                                                                   "Media")])
sample_data(particles.rar.byMediaSubTime)$Time <- as.factor(levels(sample_data(particles.rar)$Time)[get_variable(particles.rar.byMediaSubTime, 
                                                                                                                 "Time")])
# ..... reorder the levels for time in both particles.rar, particles.rar.bySubTime and particles.rar.byMediaSubTime
particles.rar@sam_data$Time=as.factor(particles.rar@sam_data$Time)
idx=sort(as.numeric(levels(particles.rar@sam_data$Time)),index.return=TRUE)$ix # reorder levels
particles.rar@sam_data$Time=factor(particles.rar@sam_data$Time,
                                   levels(particles.rar@sam_data$Time)[idx]) # reorder levels


particles.rar.bySubTime@sam_data$Time=as.factor(particles.rar.bySubTime@sam_data$Time)
idx=sort(as.numeric(levels(particles.rar.bySubTime@sam_data$Time)),index.return=TRUE)$ix # reorder levels
particles.rar.bySubTime@sam_data$Time=factor(particles.rar.bySubTime@sam_data$Time,
                                             levels(particles.rar.bySubTime@sam_data$Time)[idx]) # reorder levels

particles.rar.byMediaSubTime@sam_data$Time=as.factor(particles.rar.byMediaSubTime@sam_data$Time)
idx=sort(as.numeric(levels(particles.rar.byMediaSubTime@sam_data$Time)),index.return=TRUE)$ix # reorder levels
particles.rar.byMediaSubTime@sam_data$Time=factor(particles.rar.byMediaSubTime@sam_data$Time,
                                                  levels(particles.rar.byMediaSubTime@sam_data$Time)[idx]) # reorder levels

# .... transform to proportions for SubTime set
particles.rar.bySubTime.freq= transform_sample_counts(particles.rar.bySubTime, function(x) 100 * x/sum(x))
particles.rar.byMediaSubTime.freq= transform_sample_counts(particles.rar.byMediaSubTime, function(x) 100 * x/sum(x))

#particles.rar.bySubTime.freq.19 = prune_taxa(topESVs, particles.rar.bySubTime) # Only if you want to prune remainder taxa
particles.rar.bySubTime.freq.19 =particles.rar.bySubTime.freq # Otherwise
particles.rar.byMediaSubTime.freq.19 =particles.rar.byMediaSubTime.freq # Otherwise

# .... Finally plot
plotTitle=paste("BarPlot_diversity-",Analysis,"_MediaSubAndTime_ByGenus.pdf",sep="")
pdf(file=plotTitle,width = 22,height = 10)
label.substrate=as.character(levels(particles.rar.byMediaSubTime.freq.19@sam_data$Substrate))
title = " "
p = plot_bar(particles.rar.byMediaSubTime.freq.19, "Time", fill = "genus", title = title) + coord_flip() + 
  labs(colour = "family")+scale_fill_manual(values = col_vector)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~MediaSub,nrow=2,ncol=3)+ylab("Percentage of Sequences")
#p + facet_wrap(~MediaSub,nrow=2,ncol=3,label_value(label.substrate))+ylab("Percentage of Sequences") # gives Error in match.arg(scales, c("fixed", "free_x", "free_y", "free")) : 'arg' must be NULL or a character vector
dev.off()

# --- Perform some ordinantion analysis 
plotTitle=paste("PlotPCoA_BrayCurtis-",Analysis,"_SubstratesAndTime.pdf",sep="")
pdf(file=plotTitle,width=20,height=5)
title="PCoA of Bray-Curtis distance"
particles_ord = ordinate(particles.rar, "PCoA", "bray")
p = plot_ordination(particles.rar, particles_ord, color = "Time")
p = p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~Substrate,nrow=1)
dev.off()

plotTitle=paste("PlotPCoA_BrayCurtis-",Analysis,"_TimeAndSubstrate.pdf",sep="")
pdf(file=plotTitle,width=15,height=6)
title="PCoA of Bray-Curtis distance"
particles_ord = ordinate(particles.rar, "PCoA", "bray")
p = plot_ordination(particles.rar, particles_ord, color = "Substrate")
p = p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~Time,nrow=2)
dev.off()


# --- Perform some ordinantion analysis with the replicas merged
plotTitle=paste("PlotPCoA_BrayCurtis-",Analysis,"_Merged_SubstratesAndTime.pdf",sep="")
pdf(file=plotTitle,width=20,height=5)
title="PCoA of Bray-Curtis distance"
particles_ord = ordinate(particles.rar.bySubTime, "PCoA", "bray")
p = plot_ordination(particles.rar.bySubTime, particles_ord, color = "Time")
p = p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~Substrate,nrow=1)
dev.off()

# --- Perform some ordinantion analysis with the replicas merged
plotTitle=paste("PlotPCoA_BrayCurtis-",Analysis,"_Merged_MediaAndSubstrateByTime.pdf",sep="")
pdf(file=plotTitle,width=15,height=6)
title="PCoA of Bray-Curtis distance"
particles_ord = ordinate(particles.rar.byMediaSubTime, "PCoA", "bray")
p = plot_ordination(particles.rar.byMediaSubTime, particles_ord, color = "Substrate",shape = "Media")
p = p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~Time,nrow=2)
dev.off()

# --- Perform some ordinantion analysis with the replicas merged
plotTitle=paste("PlotPCoA_BrayCurtis-",Analysis,"_Merged_TimeAndSubstrate.pdf",sep="")
pdf(file=plotTitle,width=15,height=6)
title="PCoA of Bray-Curtis distance"
particles_ord = ordinate(particles.rar.byMediaSubTime, "PCoA", "bray")
p = plot_ordination(particles.rar.byMediaSubTime, particles_ord, color = "Time")
p = p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)+ 
  theme(axis.text = element_text(size=16),axis.title = element_text(size=20),
        legend.text = element_text(size=14),legend.title = element_text(size=20),legend.title.align =0.4,
        strip.text = element_text(size=14))
p + facet_wrap(~MediaSub,nrow=2)
dev.off()

# --- Create new vectors with ecological preferences

Ntimes=length(sample_data(particles.rar)$Time)

for(i in 1:Ntimes){
  if((sample_data(particles.rar)$Time[i] == 0)||(sample_data(particles.rar)$Time[i] == 12)){
    sample_data(particles.rar)$EcoTime[i]="Origin" 
  }else if((sample_data(particles.rar)$Time[i] == 24)||(sample_data(particles.rar)$Time[i] == 48)||
           (sample_data(particles.rar)$Time[i] == 60)){
    sample_data(particles.rar)$EcoTime[i]="Early" 
  }else if((sample_data(particles.rar)$Time[i] == 72)||(sample_data(particles.rar)$Time[i] == 84)){
    sample_data(particles.rar)$EcoTime[i]="Mid" 
  }else{
    sample_data(particles.rar)$EcoTime[i]="Late" 
  }
}

sample_data(particles.rar)$SubEcoTime <- paste0(sample_data(particles.rar)$Substrate, sample_data(particles.rar)$EcoTime)
particles.rar.byEcoTime=merge_samples(particles.rar, "EcoTime") # And proceed as above, first merge

particles.rar.Origin=prune_samples(sample_data(particles.rar)$EcoTime == "Origin", particles.rar)
particles.rar.Early=prune_samples(sample_data(particles.rar)$EcoTime == "Early", particles.rar)
particles.rar.Mid=prune_samples(sample_data(particles.rar)$EcoTime == "Mid", particles.rar)
particles.rar.Late=prune_samples(sample_data(particles.rar)$EcoTime == "Late", particles.rar)

# --- Perform some ordinantion analysis now with the samples independently, let's see if it is consistent with Pankaj prediction

# .... Origin times
plotTitle="PlotPCoA_BrayCurtis_SubstratesAndOrigin.pdf"
pdf(file=plotTitle,width=12,height=10)
title="PCoA of Bray-Curtis distance, \"Origin\" samples "
particles_ord = ordinate(particles.rar.Origin, "PCoA", "bray")
p = plot_ordination(particles.rar.Origin, particles_ord, color = "Substrate")
p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)
dev.off()

# .... Early times
plotTitle="PlotPCoA_BrayCurtis_SubstratesAndEarly.pdf"
pdf(file=plotTitle,width=12,height=10)
title="PCoA of Bray-Curtis distance, \"Early\" samples "
particles_ord = ordinate(particles.rar.Early, "PCoA", "bray")
p = plot_ordination(particles.rar.Early, particles_ord, color = "Substrate")
p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)
dev.off()

# .... Mid times
plotTitle="PlotPCoA_BrayCurtis_SubstratesAndMid.pdf"
pdf(file=plotTitle,width=12,height=10)
title="PCoA of Bray-Curtis distance, \"Mid\" samples "
particles_ord = ordinate(particles.rar.Mid, "PCoA", "bray")
p = plot_ordination(particles.rar.Mid, particles_ord, color = "Substrate")
p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)
dev.off()


# .... Later times
plotTitle="PlotPCoA_BrayCurtis_SubstratesAndLate.pdf"
pdf(file=plotTitle,width=12,height=10)
title="PCoA of Bray-Curtis distance, \"Late\" samples "
particles_ord = ordinate(particles.rar.Late, "PCoA", "bray")
p = plot_ordination(particles.rar.Late, particles_ord, color = "Substrate")
p + geom_point(size = 4, alpha = 0.7) + ggtitle(title)
dev.off()
