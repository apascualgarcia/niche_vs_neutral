###################################
#  taxonomy_betadiversity_envfit.R
##################################
# Community composition of the tree holes (K = 10 metacommunities):
#  * barplots of phyla/families/genera per tree hole, ordered by time
#    (rarefied data);
#  * robust Aitchison PCoA (unrarefied data, mitochondria/chloroplasts removed)
#    and envfit of the environmental variables on pairs of PCoA axes;
#  * boxplots of tree-hole variables per metacommunity and scatterplots of the
#    PCoA axes against depth, depth:surface and volume.
# Based on niche_envir_K10_git.R (J. Ramoneda, 2026). The previous version, which
# applied the CLR across samples instead of within samples, is in
# src/legacy/taxonomy_betadiversity_envfit_rarefiedCLR.R.
#
# Input:  1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv
#         1_dada2/otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv
#         1_dada2/taxonomy.list
#         1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output: 1_samples_properties/Metadata_singleTree_wDirich-K10_phyloseq.tsv
#         4.1_Niche_Envir/figures/ (figures, envfit tables, PCoA scores)
###################################

rm(list=ls()) # clean/reset environment 

library(stringi)
library(stringr)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(here)

# READ IN DATA ------------------------------------------------------------

#sequence table from DADA2 
seqtab_treeholes <- read.table(here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"),
                               sep = "\t", header = T)
names(seqtab_treeholes)[1] = "X" # first column is labelled "Taxon" (needed by the Dirichlet code)
rownames(seqtab_treeholes) = seqtab_treeholes$X
head(seqtab_treeholes)[1:10,1:10]

#ASV taxonomy
taxa_wsp <- read.csv(here("1_dada2", "taxonomy.list"))

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
samdf_treeholes <- read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"), header = T, sep = "\t")
samdf_treeholes <- samdf_treeholes[, !grepl("^X(\\.\\d+)?$", names(samdf_treeholes))] # empty column from the trailing tab
head(samdf_treeholes)

names(samdf_treeholes)[names(samdf_treeholes) == "ID.Sample"] = "#SampleID"
samdf_treeholes[1:10,1:10]

write.table(samdf_treeholes, here("1_samples_properties", "Metadata_singleTree_wDirich-K10_phyloseq.tsv"), sep = "\t", quote = F, row.names = F)
#write.table(otu_table_treeholes_phyloseq, here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv"), sep = "\t", quote = F, row.names = F)

###IMPORT AS PHYLOSEQ
mapfile = here("1_samples_properties", "Metadata_singleTree_wDirich-K10_phyloseq.tsv")
otufile = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv")

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)
physeq

###############################################################
# REMOVE NON-BACTERIAL SEQUENCES FIRST (before rarefying)
###############################################################

Dt.unrare <- physeq %>% subset_taxa(Family != "Mitochondria")
Dt.unrare <- Dt.unrare %>% subset_taxa(Order != "Chloroplast")
Dt.unrare

###############################################################
# RAREFIED VERSION (used for barplots and boxplots only)
###############################################################

# With rngseed = TRUE phyloseq calls set.seed(1) internally (overriding any previous
# set.seed), so the published figures used seed 1
set.seed(1)
Dt <- rarefy_even_depth(Dt.unrare, sample.size = min(sample_sums(Dt.unrare)),
                        rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
Dt

##################################################################################################################
###COMMUNITY COMPOSITION PLOT DISPLAYING THE TREEHOLES GROUPED BY METACOMMUNITY AND THE SAMPLES ORDERED BY TIME###
##################################################################################################################

## Phylum level
glom <- tax_glom(Dt, taxrank = 'Phylum')
glom 

plotdata <- psmelt(glom) 
plotdata$Phylum <- as.character(plotdata$Phylum)
plotdata = plotdata[complete.cases(plotdata$Phylum),]
length(unique(plotdata$Phylum))
head(plotdata)

xx = aggregate(plotdata$Abundance, list(Phylum = plotdata$Phylum), sum)
names(xx)[names(xx) == "x"] = "Abundance"
xx$Abundance.rel = xx$Abundance / sum(xx$Abundance)
xx = subset(xx, Abundance.rel > 0.01) #phyla with total abundance <1% will be relabelled
head(xx)

binned.plot.all = plotdata[plotdata$Phylum %in% xx$Phylum,]
binned.plot.other = plotdata[!plotdata$Phylum %in% xx$Phylum,]
binned.plot.other$Phylum = rep("< 1% abund.", nrow(binned.plot.other))
binned.plot.all = as.data.frame(rbind(binned.plot.all, binned.plot.other))
head(binned.plot.all)

binned.plot.all$Phylum <- factor(binned.plot.all$Phylum,                                    
                                 levels = c(xx$Phylum, "< 1% abund."))

myCol = c("#75151E", "#C7B446", "#424632", "#C1876B", "#A18594", "#35682D", "#686C5E", "#C6A664", "#0E294B", "#606E8C", "black")

binned.plot.all = binned.plot.all[order(binned.plot.all$R.date),]
binned.plot.all$Sample = as.factor(binned.plot.all$Sample)
binned.plot.all = binned.plot.all[!binned.plot.all$Location == "CONTROL",]

p = ggplot(binned.plot.all, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + scale_fill_discrete(type = myCol) + xlab("")  + ylab("Relative abundance") + labs(fill = "Phylum") +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size =14), legend.position = "top",
        axis.title.x = element_text(size =14, face = "bold", vjust = -0.1), axis.title.y = element_text(size =16, vjust = 2, face = "bold"), strip.text = element_text(size = 12),  
        legend.text = element_text(size=10),legend.title = element_text(size=14, face = "bold"), axis.ticks.y =  element_blank(), strip.background =element_rect(fill="white"),
        plot.title = element_text(size=12), strip.text.y = element_text(size = 8, colour = "black")) +
  scale_y_continuous(expand=c(0,0))
p 

f = p + facet_grid(~ factor(Location, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")), scales = "free_x")
f

tiff(here("4.1_Niche_Envir", "figures", "Barplot_composition_phyla_treeholes_xtime_new.tiff"), units="cm", width=24, height=12, res=300)
f
dev.off()

pdf(here("4.1_Niche_Envir", "figures", "Barplot_composition_phyla_treeholes_xtime_new.pdf"), width=24/2.54, height=12/2.54)
f
dev.off()


##Family level
glom.f <- tax_glom(Dt, taxrank = 'Family')
glom.f 

plotdata <- psmelt(glom.f) 
plotdata$Family <- as.character(plotdata$Family)
plotdata = plotdata[complete.cases(plotdata$Family),]
length(unique(plotdata$Family))
head(plotdata)

xx = aggregate(plotdata$Abundance, list(Family = plotdata$Family), sum)
names(xx)[names(xx) == "x"] = "Abundance"
xx$Abundance.rel = xx$Abundance / sum(xx$Abundance)
xx = subset(xx, Abundance.rel > 0.01) #families with total abundance <1% will be relabelled
head(xx)

binned.plot.all = plotdata[plotdata$Family %in% xx$Family,]
binned.plot.other = plotdata[!plotdata$Family %in% xx$Family,]
binned.plot.all = binned.plot.all[binned.plot.all$Family != "NA",]
head(binned.plot.all)

binned.plot.other$Family = rep("< 1% abund.", nrow(binned.plot.other))
binned.plot.all = as.data.frame(rbind(binned.plot.all, binned.plot.other))
head(binned.plot.all)

binned.plot.all$Family <- factor(binned.plot.all$Family,                                    
                                 levels = c(xx$Family, "< 1% abund."))

library(RColorBrewer)

palette_24 <- brewer.pal(11, "Paired")
palette_24 <- c(palette_24, brewer.pal(11, "Set3"))

myCol = c(palette_24, "black")

binned.plot.all = binned.plot.all[order(binned.plot.all$R.date),]
binned.plot.all$Sample = as.factor(binned.plot.all$Sample)
binned.plot.all = binned.plot.all[!binned.plot.all$Location == "CONTROL",]

p = ggplot(binned.plot.all, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + scale_fill_discrete(type = myCol) + xlab("Samples (sorted by sampling date)")  + ylab("Relative abundance") + labs(fill = "Family") +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size =10), legend.position = "top",
        axis.title.x = element_text(size =12, face = "bold", vjust = -0.1), axis.title.y = element_text(size =12, vjust = 3, face = "bold"), strip.text = element_text(size = 10),  
        legend.text = element_text(size=8),legend.title = element_text(size=10), axis.ticks.y =  element_blank(), strip.background =element_rect(fill="white"),
        plot.title = element_text(size=12), strip.text.y = element_text(size = 8, colour = "black")) +
  scale_y_continuous(expand=c(0,0))
p 

f = p + facet_grid(~ factor(Location, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")), scales = "free_x")
f

tiff(here("4.1_Niche_Envir", "figures", "Barplot_composition_family_treeholes_xtime_new.tiff"), units="cm", width=24, height=12, res=300)
f
dev.off()

pdf(here("4.1_Niche_Envir", "figures", "Barplot_composition_family_treeholes_xtime_new.pdf"), width=24/2.54, height=12/2.54)
f
dev.off()


##Genus level
glom.g <- tax_glom(Dt, taxrank = 'Genus')
glom.g 

plotdata <- psmelt(glom.g) 
plotdata$Genus <- as.character(plotdata$Genus)
plotdata = plotdata[complete.cases(plotdata$Genus),]
length(unique(plotdata$Genus))
head(plotdata)

xx = aggregate(plotdata$Abundance, list(Genus = plotdata$Genus), sum)
names(xx)[names(xx) == "x"] = "Abundance"
xx$Abundance.rel = xx$Abundance / sum(xx$Abundance)
xx = subset(xx, Abundance.rel > 0.001) #genera with total abundance <0.1% will be relabelled
head(xx)

binned.plot.all = plotdata[plotdata$Genus %in% xx$Genus,]
binned.plot.other = plotdata[!plotdata$Genus %in% xx$Genus,]
binned.plot.all = binned.plot.all[binned.plot.all$Genus != "NA",]
head(binned.plot.all)

binned.plot.other$Genus = rep("< 0.1% abund.", nrow(binned.plot.other))
binned.plot.all = as.data.frame(rbind(binned.plot.all, binned.plot.other))
head(binned.plot.all)

binned.plot.all$Genus <- factor(binned.plot.all$Genus,                                    
                                levels = c(xx$Genus, "< 0.1% abund."))

library(colorspace)

palette_93 <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2"),
  brewer.pal(9, "Set1"),
  brewer.pal(11, "Spectral"),
  brewer.pal(8, "Pastel2"),
  rainbow(37, start = 0.1, end = 0.9)
)

myCol = c(palette_93[1:92], "black")

binned.plot.all = binned.plot.all[order(binned.plot.all$R.date),]
binned.plot.all$Sample = as.factor(binned.plot.all$Sample)
binned.plot.all = binned.plot.all[!binned.plot.all$Location == "CONTROL",]

p = ggplot(binned.plot.all, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + theme_bw() + scale_fill_discrete(type = myCol) + xlab("")  + ylab("Relative abundance") + labs(fill = "Genus") +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size =10), legend.position = "none",
        axis.title.x = element_text(size =12, face = "bold", vjust = -0.1), axis.title.y = element_text(size =12, vjust = 3, face = "bold"), strip.text = element_text(size = 10),  
        legend.text = element_text(size=8),legend.title = element_text(size=10), axis.ticks.y =  element_blank(), strip.background =element_rect(fill="white"),
        plot.title = element_text(size=12), strip.text.y = element_text(size = 8, colour = "black")) +
  scale_y_continuous(expand=c(0,0))
p 

f = p + facet_grid(~ factor(Location, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")), scales = "free_x")
f

tiff(here("4.1_Niche_Envir", "figures", "Barplot_composition_genus_nolegend_treeholes_xtime_new.tiff"), units="cm", width=24, height=12, res=300)
f
dev.off()

pdf(here("4.1_Niche_Envir", "figures", "Barplot_composition_genus_nolegend_treeholes_xtime_new.pdf"), width=24/2.54, height=12/2.54)
f
dev.off()


##########################################################################################
### LINEAR ASSOCIATIONS BETWEEN COMMUNITY COMPOSITION AND ENVIRONMENTAL FACTORS (ENVFIT)
##########################################################################################

library(vegan)
library(compositions)

# Shared helper for saving figures
save_tiff_pdf <- function(plot_obj, filename, width, height) {
  tiff(here("4.1_Niche_Envir", "figures", paste0(filename, ".tiff")), units = "cm", width = width, height = height, res = 300)
  print(plot_obj)
  dev.off()
  
  pdf(here("4.1_Niche_Envir", "figures", paste0(filename, ".pdf")), width = width/2.54, height = height/2.54)
  print(plot_obj)
  dev.off()
}

###############################################################
# Select samples (UNRAREFIED data, mito/chloroplast removed)
###############################################################

dataset.envfit <- subset_samples(Dt.unrare, Sample.type == "THS" & Location != "CONTROL")
dataset.envfit

###############################################################
# Metadata
###############################################################

dataset.envfit_env <- as(sample_data(dataset.envfit), "data.frame")
dataset.envfit_env$Depth.surface <- dataset.envfit_env$Depth / dataset.envfit_env$SA
head(dataset.envfit_env)

###############################################################
# OTU table (samples as rows)
###############################################################

otu.mat <- as(otu_table(dataset.envfit), "matrix")
if (taxa_are_rows(dataset.envfit)) {
  otu.mat <- t(otu.mat)
}
dim(otu.mat)

###############################################################
# AITCHISON DISTANCE
# "aitchison": Euclidean distance between CLR profiles, the CLR computed
#   within each sample after adding a pseudocount of 1 (full-rank, so the
#   % of variance of the PCoA axes is the share of the total compositional
#   variance).
# "robust.aitchison": vegan computes the rCLR of the non-zero values and, by
#   default, treats the zeros as missing and runs a rank-3 matrix completion
#   (OptSpace), as in DEICODE/RPCA. The completed rank-3 matrix replaces ALL the
#   entries, observed ones included. The distance then lives in 3 dimensions: the first three PCoA
#   axes always explain 100% and are a rotation of the leading CLR axes. vegdist
#   also drops the sample names in that case, so they are restored below.
###############################################################

dist.method <- "aitchison" # or "robust.aitchison"
if (dist.method == "aitchison") {
  ait_dist <- vegan::vegdist(otu.mat, method = "aitchison", pseudocount = 1)
} else {
  ait_dist <- vegan::vegdist(otu.mat, method = "robust.aitchison")
  attr(ait_dist, "Labels") <- rownames(otu.mat)
}
as.matrix(ait_dist)[1:10, 1:10]

###############################################################
# PCoA (single ordination used for envfit AND all plots)
###############################################################

pcoa_ait <- cmdscale(ait_dist, k = 3, eig = TRUE)
pcoa_ait

# Variance explained (relative to the sum of positive eigenvalues)
pct_exp <- 100 * pcoa_ait$eig / sum(pcoa_ait$eig[pcoa_ait$eig > 0])
pct_exp[1:3]

scores_all <- as.data.frame(pcoa_ait$points)
colnames(scores_all) <- c("Dim1", "Dim2", "Dim3")
head(scores_all)

# Make sure metadata rows are in the same order as the ordination
dataset.envfit_env <- dataset.envfit_env[rownames(scores_all), ]

# Scores re-used by immigration_speciation_rates_analysis.R
write.csv(data.frame(ID.Sample = rownames(scores_all), scores_all),
          here("4.1_Niche_Envir", "figures", "pcoa_aitchison_scores.csv"), row.names = F)
write.csv(data.frame(Axis = 1:3, Percent.variance = pct_exp[1:3]),
          here("4.1_Niche_Envir", "figures", "pcoa_aitchison_variance.csv"), row.names = F)

###############################################################
# ENVFIT VARIABLES
# (Depth.surface is created above, so it is already numeric
#  and does not need to be appended a second time)
###############################################################

numeric_cols <- names(dataset.envfit_env)[sapply(dataset.envfit_env, is.numeric)]

vars.envfit <- numeric_cols[numeric_cols %in% c("Moisture", "Max.Lum", "SA", "Volume",
                                                "Depth", "Deg.North", "Air.Temp", "SOil.Temp2",
                                                "Soil.Temp4", "Rainfall", "Dought_days",
                                                "Biomass", "Dirichlet", "Depth.surface")]
vars.envfit

###############################################################
# ENVFIT looped over variables, for each pair of axes
###############################################################

run_envfit <- function(axis_scores, env, vars) {
  set.seed(1) # reproducible permutation p-values
  out <- data.frame()
  for (col in vars) {
    keep <- !is.na(env[, col])
    fit <- envfit(axis_scores[keep, ], env[keep, col], permutations = 999)
    vec_scores <- scores(fit, "vectors")
    out <- rbind(out, data.frame(Variable = col,
                                 r2 = fit$vectors$r,
                                 P = fit$vectors$pvals,
                                 A = vec_scores[1],
                                 B = vec_scores[2]))
  }
  names(out)[4:5] <- colnames(axis_scores)
  out
}

## Axes 1-2
results_df.1.2 <- run_envfit(scores_all[, c("Dim1","Dim2")], dataset.envfit_env, vars.envfit)
results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),]

write.csv(results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes1_2_2026.csv"), row.names = F)

## Axes 1-3
results_df.1.3 <- run_envfit(scores_all[, c("Dim1","Dim3")], dataset.envfit_env, vars.envfit)
results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),]

write.csv(results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes1_3_2026.csv"), row.names = F)

## Axes 2-3
results_df.2.3 <- run_envfit(scores_all[, c("Dim2","Dim3")], dataset.envfit_env, vars.envfit)
results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),]

write.csv(results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes2_3_2026.csv"), row.names = F)


###############################################################
## Boxplots with environmental variables by metacommunity
###############################################################

env.data = as(sample_data(Dt), "data.frame")
env.data$Depth.surface = env.data$Depth / env.data$SA
env.data$Metacommunity = factor(env.data$Dirichlet, levels = c(1:10))
env.data = env.data[!is.na(env.data$Metacommunity),]
head(env.data)

custom_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                    "#6D3F5B", "#1E213D", "#1C542D", "#763C28")

make_metacomm_boxplot <- function(data, yvar, ylab_text) {
  ggplot(data, aes(x = Metacommunity, y = .data[[yvar]], fill = Metacommunity)) +
    theme_classic() +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = Metacommunity), size = 2, alpha = 0.5, width = 0.2) +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "none",
          axis.title.x = element_text(size = 14, margin = margin(t = 10)),
          axis.title.y = element_text(size = 14, margin = margin(r = 10)),
          plot.title = element_text(size = 14, face = "bold", colour = "Black"),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    ylab(ylab_text) +
    xlab("Metacommunity") +
    scale_fill_manual(values = custom_palette) +
    scale_color_manual(values = custom_palette)
}

#Volume
volume.metacomm = make_metacomm_boxplot(env.data, "Volume", expression(Volume~(cm^3)))
volume.metacomm
save_tiff_pdf(volume.metacomm, "Boxplot_volume_metacommunity", 14, 12)

#Depth
depth.metacomm = make_metacomm_boxplot(env.data, "Depth", "Depth (cm)")
depth.metacomm
save_tiff_pdf(depth.metacomm, "Boxplot_depth_metacommunity", 14, 12)

#Drought
depth.drought = make_metacomm_boxplot(env.data, "Dought_days", "Drought (days)")
depth.drought
save_tiff_pdf(depth.drought, "Boxplot_drought_metacommunity", 14, 12)

#Depth / Surface
depth.surface.metacomm = make_metacomm_boxplot(env.data, "Depth.surface", "Depth:surface ratio")
depth.surface.metacomm
save_tiff_pdf(depth.surface.metacomm, "Boxplot_depth_surface_ratio_metacommunity", 14, 12)


###############################################################
### PCoA WITH METACOMMUNITIES AND ENVFIT ARROWS
### (re-uses the PCoA and envfit results computed above)
###############################################################

library(tibble)

env.pcoa <- dataset.envfit_env
env.pcoa$Metacommunity <- factor(env.pcoa$Dirichlet, levels = 1:10)

scores_df <- cbind(scores_all, env.pcoa[rownames(scores_all), ])
colnames(scores_df)[1:3] <- c("PC1", "PC2", "PC3")

scores_df$Location <- factor(scores_df$Location, levels = 1:10)
names(scores_df)[names(scores_df) == "Location"] <- "Treehole"
scores_df <- scores_df[!is.na(scores_df$Metacommunity), ]
head(scores_df)

# Relabel/clean significant arrows consistently across axis pairs
prep_significant_arrows <- function(envfit_res) {
  sig <- envfit_res %>% filter(P < 0.05)
  sig <- sig[sig$Variable != "Dirichlet", ]
  sig$Variable <- as.factor(sig$Variable)
  levels(sig$Variable)[levels(sig$Variable) == "Dought_days"] <- "Drought days"
  levels(sig$Variable)[levels(sig$Variable) == "SA"] <- "Surface area"
  levels(sig$Variable)[levels(sig$Variable) == "Depth.surface"] <- "Depth:Surface ratio"
  sig
}

# Build a PCoA biplot for a given pair of axes
make_pcoa_biplot <- function(scores_data, arrows_data, x_dim, y_dim,
                             arrow_x, arrow_y, x_lab, y_lab) {
  # envfit arrows have length sqrt(r2) (<= 1): stretch them to the range of the scores
  mult <- 0.8 * max(abs(scores_data[, c(x_dim, y_dim)])) / max(sqrt(arrows_data$r2))
  arrows_data[, c(arrow_x, arrow_y)] <- arrows_data[, c(arrow_x, arrow_y)] * mult
  ggplot(data = scores_data, aes(x = .data[[x_dim]], y = .data[[y_dim]],
                                 color = Metacommunity, size = Volume, shape = Treehole)) +
    geom_point(alpha = 0.75, stroke = 1.5) +
    geom_segment(data = arrows_data,
                 aes(x = 0, y = 0, xend = .data[[arrow_x]], yend = .data[[arrow_y]]),
                 arrow = arrow(length = unit(0.1, "cm")),
                 inherit.aes = FALSE,
                 color = "black", linewidth = 1) +
    geom_text(data = arrows_data,
              aes(x = .data[[arrow_x]] * 1.35, y = .data[[arrow_y]] * 1.35, label = Variable),
              inherit.aes = FALSE,
              color = "black", size = 5, alpha = 1) +
    scale_shape_manual(values = c(0:14, 15:18)) +
    scale_color_manual(values = custom_palette) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2.5),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    guides(
      color = guide_legend(ncol = 2),
      shape = guide_legend(ncol = 2),
      size = guide_legend(ncol = 1)
    ) +
    xlab(x_lab) +
    ylab(y_lab)
}

axis_lab <- function(i) paste0("PCoA", i, " (", round(pct_exp[i], 1), "%)")

## Axes 1-2
significant_arrows <- prep_significant_arrows(results_df.1.2)
ss1.2 <- make_pcoa_biplot(scores_df, significant_arrows,
                          "PC1", "PC2", "Dim1", "Dim2",
                          axis_lab(1), axis_lab(2))
ss1.2
save_tiff_pdf(ss1.2, "PCOA_treehole_metacommunity_k10_axes1_2_2026", 20, 14)

## Axes 1-3
significant_arrows <- prep_significant_arrows(results_df.1.3)
ss1.3 <- make_pcoa_biplot(scores_df, significant_arrows,
                          "PC1", "PC3", "Dim1", "Dim3",
                          axis_lab(1), axis_lab(3))
ss1.3
save_tiff_pdf(ss1.3, "PCOA_treehole_metacommunity_k10_axes1_3_2026", 20, 14)

## Axes 2-3
significant_arrows <- prep_significant_arrows(results_df.2.3)
ss2.3 <- make_pcoa_biplot(scores_df, significant_arrows,
                          "PC2", "PC3", "Dim2", "Dim3",
                          axis_lab(2), axis_lab(3))
ss2.3
save_tiff_pdf(ss2.3, "PCOA_treehole_metacommunity_k10_axes2_3_2026", 20, 14)


###############################################################
## Plot PCoA axes versus environmental factors
###############################################################

mds.vectors <- scores_all
mds.vectors$Sample <- rownames(mds.vectors)
head(mds.vectors)

dataset.envfit_env$Sample <- rownames(dataset.envfit_env)
pcoa.env.dbb <- merge(mds.vectors, dataset.envfit_env, by = "Sample")
head(pcoa.env.dbb)

library(ggpmisc)

##Depth
p = ggplot(pcoa.env.dbb, aes(x = Depth, y = Dim1)) +
  theme_bw() + 
  geom_point(size = 4, alpha = 0.5) + 
  geom_smooth(method = lm, se = TRUE, color = "black", fill = "grey70") + 
  stat_poly_line(color = "black") + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")),
               size = 5,
               label.y = 0.99) +  
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 16, margin = margin(r = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  xlab("Treehole depth (cm)") +
  ylab("PCoA1")
p

cor.test(pcoa.env.dbb$Dim1, pcoa.env.dbb$Depth)

save_tiff_pdf(p, "MDS1_Depth_correlation_2026", 15.24, 15.24) # 6 in = 15.24 cm


##Depth:Surface ratio
p = ggplot(pcoa.env.dbb, aes(x = Depth.surface, y = Dim3)) +
  theme_bw() + 
  geom_point(size = 4, alpha = 0.5) + 
  geom_smooth(method = lm,
              formula = y ~ log(x),
              se = TRUE,
              color = "black",
              fill = "grey70") + 
  stat_poly_line(formula = y ~ log(x), color = "black") + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")),
               formula = y ~ log(x),
               size = 5, label.y = 0.99, label.x = 0.4) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 16, margin = margin(r = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  xlab("Depth:Surface ratio") +
  ylab("PCoA3")
p

cor.test(pcoa.env.dbb$Dim3, pcoa.env.dbb$Depth.surface)

save_tiff_pdf(p, "MDS3_Depth_surface_ratio_correlation_2026", 15.24, 15.24)


##Volume
# Fit logarithmic model
model <- lm(Dim1 ~ log(Volume), data = pcoa.env.dbb)

coef_vals <- coef(model)
intercept <- round(coef_vals[1], 3)
slope <- round(coef_vals[2], 3)

model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

p = ggplot(pcoa.env.dbb, aes(x = Volume, y = Dim1)) +
  theme_bw() + 
  geom_point(size = 4, alpha = 0.5) + 
  geom_smooth(method = lm,
              formula = y ~ log(x),
              se = TRUE,
              color = "black",
              fill = "grey70") + 
  stat_poly_line(formula = y ~ log(x), color = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("y = ", intercept,
                          ifelse(slope >= 0, " + ", " - "),
                          abs(slope), "·log(x)"),
           hjust = 1.05, vjust = 1.5, size = 5) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R² = ", r2),
           hjust = 1.1, vjust = 3.2, size = 5) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("p = ", p_val),
           hjust = 1.075, vjust = 4.7, size = 5) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 16, margin = margin(r = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  xlab(expression(Volume~(cm^3))) +
  ylab("PCoA1")
p

cor.test(pcoa.env.dbb$Dim1, pcoa.env.dbb$Volume)

save_tiff_pdf(p, "MDS1_Volume_correlation_2026", 15.24, 15.24)
