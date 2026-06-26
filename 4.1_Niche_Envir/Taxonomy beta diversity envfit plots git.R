rm(list=ls()) # clean/reset environment 

library(stringi)
library(stringr)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(here)

# READ IN DATA ------------------------------------------------------------

#sequence table from DADA2 
seqtab_treeholes <- read.table(here("2_Dirichlet", "otu_tables_K10", "otu_table_treeholes_dada2_renamed_clean_red.tsv"),
                               sep = "\t", header = T)
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
samdf_treeholes <- read.csv(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.csv"), header = T)
head(samdf_treeholes)

names(samdf_treeholes)[names(samdf_treeholes) == "ID.Sample"] = "#SampleID"
samdf_treeholes[1:10,1:10]

#write.table(samdf_treeholes, here("1_samples_properties", "Metadata_singleTree_wDirich-K10_phyloseq.tsv"), sep = "\t", quote = F, row.names = F)
#write.table(otu_table_treeholes_phyloseq, here("2_Dirichlet", "otu_tables_K10", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv"), sep = "\t", quote = F, row.names = F)

###IMPORT AS PHYLOSEQ
mapfile = here("1_samples_properties", "Metadata_singleTree_wDirich-K10_phyloseq.tsv")
otufile = here("2_Dirichlet", "otu_tables_K10", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv")

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)
physeq

set.seed(123)
Dt = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                       rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
Dt

##REMOVAL OF MITOCHONDRIA AND CHLOROPLAST
Dt <- Dt %>% subset_taxa(Family != "Mitochondria")
Dt <- Dt %>% subset_taxa(Order != "Chloroplast")
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
        axis.title.x = element_text(size =14, face = "bold", vjust = -0.1), axis.title.y = element_text(size =16, vjust = 3, face = "bold"), strip.text = element_text(size = 12),  
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
  geom_bar(stat = "identity", position = "fill") + theme_bw() + scale_fill_discrete(type = myCol) + xlab("")  + ylab("Relative abundance") + labs(fill = "Family") +
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

library(RColorBrewer)
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
###LINEAR ASSOCIATIONS BETWEEN COMMUNITY COMPOSITION AND ENVIRONMENTAL FACTORS (ENVFIT)###
##########################################################################################

##transform to relative abundance
Dt.rel <- transform_sample_counts(Dt, function(OTU) OTU/sum(OTU))
Dt.rel

##ENVIFT: LINEAR ASSOCIATIONS BETWEEN ENVIRONMENTAL VARIABLES AND COMMUNITY STRUCTURE

library(vegan)
library(microbiome)

dataset.envfit = subset_samples(Dt, Sample.type == "THS" & Location != "CONTROL")
dataset.envfit

dataset.envfit_env = meta(dataset.envfit)
head(dataset.envfit_env)

#add pseudocounts
envift_transformed = t(as.data.frame(otu_table(dataset.envfit)))
envift_transformed[envift_transformed == 0] <- 0.00000000000000000001
envift_transformed[1:10,1:10]

#calculate Aitchison distance
ait_dist = vegdist(envift_transformed, method="aitchison")
as.matrix(ait_dist)[1:10,1:10]

#create pcoa from that distance matrix
pcoa_ait <- cmdscale(ait_dist, k = 3, eig = T)
pcoa_ait

#envfit looped over all numeric factors

results_df.1.2 <- data.frame()
results_df.1.3 <- data.frame()
results_df.2.3 <- data.frame()

numeric_cols <- names(dataset.envfit_env)[sapply(dataset.envfit_env, is.numeric)]
dataset.envfit_env$Depth.surface = dataset.envfit_env$Depth / dataset.envfit_env$SA

vars.envfit = c(numeric_cols[numeric_cols %in% c("Moisture", "Max.Lum", "SA", "Volume",
                                               "Depth", "Deg.North", "Air.Temp", "SOil.Temp2",
                                               "Soil.Temp4", "Rainfall", "Dought_days",
                                               "Biomass", "Dirichlet", "Depth.surface")], "Depth.surface")
vars.envfit

## Envfit PCoA1-2
for (col in vars.envfit) {
  fit <- envfit(as.data.frame(scores(pcoa_ait)[!is.na(dataset.envfit_env[, col]), ]),
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  results_df.1.2 <- rbind(results_df.1.2, data.frame(Variable = col, r2 = r2, P = p_value, Dim1 = scores[1], Dim2 = scores[2]))
}

results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),]

#write.csv(results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes1_2_2026.csv"), row.names = F)

## Envfit PCoA1-3
scores.1.3 = as.data.frame(scores(pcoa_ait))[,c(1,3)]
scores.1.3

for (col in vars.envfit) {
  fit <- envfit(scores.1.3[!is.na(dataset.envfit_env[, col]), ],
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  results_df.1.3 <- rbind(results_df.1.3, data.frame(Variable = col, r2 = r2, P = p_value, Dim1 = scores[1], Dim3 = scores[2]))
}

results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),]

#write.csv(results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes1_3_2026.csv"), row.names = F)

## Envfit PCoA2-3
scores.2.3 = as.data.frame(scores(pcoa_ait))[,c(2,3)]
scores.2.3

for (col in vars.envfit) {
  fit <- envfit(scores.2.3[!is.na(dataset.envfit_env[, col]), ],
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  results_df.2.3 <- rbind(results_df.2.3, data.frame(Variable = col, r2 = r2, P = p_value, Dim2 = scores[1], Dim3 = scores[2]))
}

results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),]

#write.csv(results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),], here("4.1_Niche_Envir", "figures", "envfit_results_aitchison_pcoa_axes2_3_2026.csv"), row.names = F)


## Boxplots/violin plots with environmental variables by metacommunity
env.data = as.data.frame(sample_data(Dt))
env.data$Depth.surface = env.data$Depth / env.data$SA
env.data$Metacommunity = factor(env.data$Dirichlet, levels = c(1:10))
env.data = env.data[!is.na(env.data$Metacommunity),]
head(env.data)

custom_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                    "#6D3F5B", "#1E213D", "#1C542D", "#763C28")

# Helper to avoid repeating the same boxplot theme/styling for each variable
make_metacomm_boxplot <- function(data, yvar, ylab_text) {
  ggplot(data, aes(x = Metacommunity, y = .data[[yvar]], fill = Metacommunity)) +
    theme_classic() +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(aes(color = Metacommunity), size = 2, alpha = 0.5, width = 0.2) +
    theme(axis.text.x = element_text(size = 12, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "null",
          axis.title.x = element_text(size = 14, margin = margin(t = 10)),
          axis.title.y = element_text(size = 14, margin = margin(r = 10)),
          plot.title = element_text(size = 14, face = "bold", colour = "Black"),
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    ylab(ylab_text) +
    xlab("Metacommunity") +
    scale_fill_manual(values = custom_palette) +
    scale_color_manual(values = custom_palette)
}

save_tiff_pdf <- function(plot_obj, filename, width, height) {
  tiff(here("4.1_Niche_Envir", "figures", paste0(filename, ".tiff")), units = "cm", width = width, height = height, res = 300)
  print(plot_obj)
  dev.off()
  
  pdf(here("4.1_Niche_Envir", "figures", paste0(filename, ".pdf")), width = width/2.54, height = height/2.54)
  print(plot_obj)
  dev.off()
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


##Make PCoA with metacommunities and envfit results as as arrows
library(tibble) 
library(compositions) 

# Step 1: CLR transformation and Aitchison distance
otu_table_clr <- apply(otu_table(Dt), 1, function(x) clr(x + 1e-5))
dist_aitchison <- dist(otu_table_clr, method = "euclidean")

# Step 2: Perform PCoA
pcoa_res <- cmdscale(dist_aitchison, k = 3, eig = TRUE)
scores_df <- as.data.frame(pcoa_res$points)
colnames(scores_df) <- c("PC1", "PC2", "PC3")
scores_df <- cbind(scores_df, env.data) # Combine with sample metadata
head(scores_df)

# Step 3: Fit environmental factors
envfit_res1.2 <- results_df.1.2
envfit_res1.3 <- results_df.1.3
envfit_res2.3 <- results_df.2.3

# Helper to relabel/clean significant-arrows data frames consistently across axis pairs
prep_significant_arrows <- function(envfit_res, dim_cols) {
  sig <- envfit_res[, c(dim_cols, "P", "Variable", "r2")] %>%
    as.data.frame() %>%
    filter(envfit_res$P < 0.05)
  
  sig = sig[sig$Variable != "Dirichlet",]
  sig$Variable = as.factor(sig$Variable)
  levels(sig$Variable)[levels(sig$Variable) == "Dought_days"] <- "Drought days"
  levels(sig$Variable)[levels(sig$Variable) == "SA"] <- "Surface area"
  levels(sig$Variable)[levels(sig$Variable) == "Depth.surface"] <- "Depth:Surface ratio"
  sig
}

# Helper to build a PCoA biplot for a given pair of axes
make_pcoa_biplot <- function(scores_data, arrows_data, x_dim, y_dim, x_lab, y_lab) {
  ggplot(data = scores_data, aes(x = .data[[x_dim]] / 100, y = .data[[y_dim]] / 100,
                                 color = Metacommunity, size = Volume, shape = Treehole)) +
    geom_point(alpha = 0.75, stroke = 1.5) +
    geom_segment(data = arrows_data,
                aes(x = 0, y = 0, xend = .data[[x_dim]], yend = .data[[y_dim]]),
                arrow = arrow(length = unit(0.1, "cm")),
                inherit.aes = FALSE,
                color = "black", linewidth = 1) +
    geom_text(data = arrows_data,
              aes(x = .data[[x_dim]] * 1.35, y = .data[[y_dim]] * 1.35, label = Variable),
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

# Prepare scores_df for plotting (shared across axis-pair plots)
names(scores_df)[names(scores_df) %in% "Dirichlet"] = "Metacommunity"
scores_df$Metacommunity = factor(scores_df$Metacommunity, levels = c(1:10))
scores_df$Location = factor(scores_df$Location, levels = c(1:10))
names(scores_df)[names(scores_df) == "Location"] = "Treehole"
scores_df = scores_df[!is.na(scores_df$Metacommunity),]
scores_df = scores_df[,-ncol(scores_df)]

## Axes 1-2
significant_arrows = prep_significant_arrows(envfit_res1.2, c("Dim1","Dim2"))

ss1.2 = make_pcoa_biplot(scores_df, significant_arrows, "PC1", "PC2",
                         paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
                         paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)"))
ss1.2
save_tiff_pdf(ss1.2, "PCOA_treehole_metacommunity_k10_axes1_2_2026", 20, 14)

## Axes 1-3
significant_arrows = prep_significant_arrows(envfit_res1.3, c("Dim1","Dim3"))

ss1.3 = make_pcoa_biplot(scores_df, significant_arrows, "PC1", "PC3",
                         paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
                         paste0("PCoA3 (", round(pcoa_res$eig[3] / sum(pcoa_res$eig) * 100, 1), "%)"))
ss1.3
save_tiff_pdf(ss1.3, "PCOA_treehole_metacommunity_k10_axes1_3_2026", 20, 14)

## Axes 2-3
significant_arrows = prep_significant_arrows(envfit_res2.3, c("Dim2","Dim3"))

ss2.3 = make_pcoa_biplot(scores_df, significant_arrows, "PC2", "PC3",
                         paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)"),
                         paste0("PCoA3 (", round(pcoa_res$eig[3] / sum(pcoa_res$eig) * 100, 1), "%)"))
ss2.3
save_tiff_pdf(ss2.3, "PCOA_treehole_metacommunity_k10_axes2_3_2026", 20, 14)


##Plot PCoA axes versus environmental factors
mds.vectors = as.data.frame(scores(pcoa_res))
mds.vectors$Sample = rownames(mds.vectors)
head(mds.vectors)

dataset.envfit_env$Sample = rownames(dataset.envfit_env)
head(dataset.envfit_env)

pcoa.env.dbb = merge(mds.vectors, dataset.envfit_env, by = "Sample")
head(pcoa.env.dbb)

library(ggpmisc)

##Depth
p = ggplot(pcoa.env.dbb, aes(x = Depth, y = Dim1/100)) +
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

cor.test(pcoa.env.dbb$Dim1, pcoa.env.dbb$Depth) #r = 0.868, P < 0.001

save_tiff_pdf(p, "MDS1_Depth_correlation_2026", 15.24, 15.24) # 6 in ~ 15.24 cm


##Depth:Surface ratio
p = ggplot(pcoa.env.dbb, aes(x = Depth.surface, y = Dim3/100)) +
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

cor.test(pcoa.env.dbb$Dim3, pcoa.env.dbb$Depth.surface) #r = -0.676, P < 0.001

save_tiff_pdf(p, "MDS3_Depth_surface_ratio_correlation_2026", 15.24, 15.24)


##Volume
p = ggplot(pcoa.env.dbb, aes(x = Volume, y = Dim1/100)) +
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
               size = 5, label.y = 0.99, label.x = 0.99) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 16, margin = margin(r = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  xlab(expression(Volume~(cm^3))) +
  ylab("PCoA1")
p

cor.test(pcoa.env.dbb$Dim1, pcoa.env.dbb$Volume) #r = 0.693, P < 0.001

save_tiff_pdf(p, "MDS1_Volume_correlation_2026", 15.24, 15.24)
