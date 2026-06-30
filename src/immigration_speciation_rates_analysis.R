rm(list=ls()) # clean/reset environment 

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(purrr)
library(here)

# 1 Load metadata, immigration rates, speciation rates and merge in single dataframe
meta_treeholes = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv"), header = T)
meta_treeholes$Depth.surface.ratio = meta_treeholes$Depth / meta_treeholes$SA
head(meta_treeholes)

# Define environmental variables that are specific to each sample and can thus be associated with immigration rates
env.vars = c("Moisture", "Max.Lum", "Air.Temp", "Grass.Temp", "Soil.Temp2", "Soil.Temp4", "Air.Temp.max",
             "Air.Temp.min", "Grass.Temp.max", "Grass.Temp.min",  "Rainfall", "Location", "Dought_days")

# Multiple regression with treehole as factor
mult.reg <- as.formula(
  paste("Immigrate_med ~", paste(env.vars, collapse = " + "))
)
model <- lm(mult.reg, data = meta_treeholes)
summary(model)

# ---------------------------
# Step 1: Summarize data
# ---------------------------
names(meta_treeholes)[names(meta_treeholes) == "Location"] = "Treehole"
names(meta_treeholes)[names(meta_treeholes) == "Dirichlet"] = "Metacommunity"
meta_treeholes$Metacommunity = factor(meta_treeholes$Metacommunity, levels = 1:10)
meta_treeholes$Treehole = factor(meta_treeholes$Treehole, levels = 1:10)

# ---------------------------
# Fit model
# ---------------------------
model <- lm(Immigrate_med ~ log(Depth.surface.ratio), data = meta_treeholes)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R2 and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·log(x)\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes, aes(x = Depth.surface.ratio, y = Immigrate_med)) +
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
  "#6D3F5B", "#1E213D", "#1C542D", "#763C28")) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    se = TRUE,
    color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1, vjust = 1,
    size = 4
  ) +
  
  labs(
    x = "Depth:Surface ratio",
    y = "Immigration rate"
  ) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
    plot.title = element_text(face = "bold")
  ) + xlim(0, 0.3) +
  scale_shape_manual(values = c(0:14, 15:18)) +
  guides(shape = guide_legend(ncol = 2, byrow = FALSE), color = guide_legend(ncol = 2, byrow = FALSE))
p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_depth_surface_ratio_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_depth_surface_ratio_metacommunity.pdf"),
    width = 14/2.54, height = 12/2.54)
p
dev.off()


## Same with drought days
max_drought <- meta_treeholes %>%
  group_by(Treehole) %>%
  summarise(
    max_drought_days = max(Dought_days, na.rm = TRUE),
    .groups = "drop"
  )
max_drought

meta_treeholes_drought = merge(meta_treeholes, max_drought, by = "Treehole")
head(meta_treeholes_drought)

# ---------------------------
# Fit model
# ---------------------------
model <- lm(Immigrate_med ~ log(max_drought_days), data = meta_treeholes_drought)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R2 and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·log(x)\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes_drought,
            aes(x = max_drought_days, y = Immigrate_med)) +
  
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3",
                                "#497E76", "#D6AE01", "#6D3F5B", "#1E213D",
                                "#1C542D", "#763C28")) +
  
  geom_smooth(method = "lm",
              formula = y ~ log(x),
              se = TRUE,
              color = "black") +
  
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1, vjust = 1,   # aligned with reference plot
           size = 4) +
  
  labs(
    x = "Drought prevalence (total drought days)",
    y = "Immigration rate"
  ) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")) +
  
  scale_shape_manual(values = c(0:14, 15:18)) +
  
  guides(shape = guide_legend(ncol = 2, byrow = FALSE),
         color = guide_legend(ncol = 2, byrow = FALSE))
p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


## Same with volume

model <- lm(Immigrate_med ~ log(Volume), data = meta_treeholes)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R2 and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·log(x)\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes, aes(x = Volume, y = Immigrate_med)) +
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                                "#6D3F5B", "#1E213D", "#1C542D", "#763C28")) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    se = TRUE,
    color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1, vjust = 1,
    size = 4
  ) +
  
  labs(
    x = expression(Volume~(cm^3)),
    y = "Immigration rate") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")
  ) + #xlim(0, 0.3) +
  scale_shape_manual(values = c(0:14, 15:18)) +
  guides(shape = guide_legend(ncol = 2, byrow = FALSE), color = guide_legend(ncol = 2, byrow = FALSE))
p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_volume_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Immigration_volume_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


# Depth:Surface association with drought prevalence
max_drought_depth <- meta_treeholes_drought %>%
  dplyr::distinct(Treehole, Depth.surface.ratio, max_drought_days)
max_drought_depth

p = ggplot(meta_treeholes_drought,
       aes(x = Depth.surface.ratio, y = max_drought_days)) +
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3",
                                "#497E76", "#D6AE01", "#6D3F5B", "#1E213D",
                                "#1C542D", "#763C28")) +
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  annotate("text",
           x = Inf, 
           y = Inf,
           label = "Pearson's r = -0.778,\nP < 0.0001",
           hjust = 1, 
           vjust = 1,   # aligned with your other plots
           size = 4) +
  labs(
    x = "Depth:Surface ratio",
    y = "Drought prevalence (total drought days)"
  ) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")) +
  
  scale_shape_manual(values = c(0:14, 15:18)) +
  
  guides(
    shape = guide_legend(ncol = 2, byrow = FALSE),
    color = guide_legend(ncol = 2, byrow = FALSE)
  )
p

cor.test(meta_treeholes_drought$max_drought_days, meta_treeholes_drought$Depth.surface.ratio)

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Depth_surface_ratio_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Depth_surface_ratio_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


## Include immigration and speciation rates as factors in PCoA
rm(list=ls()) # clean/reset environment 

library(stringi)
library(stringr)
library(phyloseq)
library(ggplot2)
library(here)

###IMPORT AS PHYLOSEQ
mapfile = here("1_samples_properties", "Metadata_singleTree_wDirich-K10_phyloseq.tsv")
otufile = here("2_Dirichlet", "otu_tables_K10", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv")

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)
physeq

##RAREFACTION
#Rarefaction plot
#source("Rfunctions.R")
#rarefaction.plot = ggrare(physeq, step = 200, se = FALSE) + xlim(c(0,max(sample_sums(physeq)))) + ylim(0,2000) + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y=element_text(size=20), legend.position = "null",
#                                                                                                                                    axis.title.x = element_text(size = 20, margin = margin(r = 15)), axis.title.y = element_text(size =20, margin = margin(r = 10))) + xlab("Number of reads")

#tiff(here("2_Dirichlet_ImmigrationSpeciation", "Rarefaction_plot_treeholes.tiff"), units="cm", width=15, height=15, res=300)
#rarefaction.plot
#dev.off()

set.seed(123)
Dt = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                       rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
Dt

##REMOVAL OF MITOCHONDRIA AND CHLOROPLAST
Dt <- Dt %>% subset_taxa(Family!= "Mitochondria")
Dt <- Dt %>% subset_taxa(Order!= "Chloroplast")
Dt

##transform to relative abundance
Dt.rel <- transform_sample_counts(Dt, function(OTU) OTU/sum(OTU))
Dt.rel


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

# Create an empty data frame to store the results
results_df.1.2 <- data.frame()
results_df.1.3 <- data.frame()
results_df.2.3 <- data.frame()

# Get the names of columns containing numeric variables
numeric_cols <- names(dataset.envfit_env)[sapply(dataset.envfit_env, is.numeric)]
dataset.envfit_env$Depth.surface = dataset.envfit_env$Depth / dataset.envfit_env$SA

vars.envfit = c(numeric_cols[numeric_cols %in% c("Moisture", "Max.Lum", "SA", "Volume",
                                                 "Depth", "Deg.North", "Air.Temp", "SOil.Temp2",
                                                 "Soil.Temp4", "Rainfall", "Dought_days",
                                                 "Biomass", "Dirichlet", "Depth.surface")], "Depth.surface")
vars.envfit

## Envfit PCoA1-2
# Loop over the numeric columns
for (col in vars.envfit) {
  # Run the envfit function on the current column
  fit <- envfit(as.data.frame(scores(pcoa_ait)[!is.na(dataset.envfit_env[, col]), ]),
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  # Extract the r2 and Pr(>r) values from the fit object
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  # Add the results to the data frame
  results_df.1.2 <- rbind(results_df.1.2, data.frame(Variable = col, r2 = r2, P = p_value, Dim1 = scores[1], Dim2 = scores[2]))
}

results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),]

#write.csv(results_df.1.2[order(-results_df.1.2$r2, results_df.1.2$P),], here("2_Dirichlet_ImmigrationSpeciation", "envfit_results_aitchison_pcoa_axes1_2_2026.csv"), row.names = F)

## Envfit PCoA1-3
# Loop over the numeric columns
scores.1.3 = as.data.frame(scores(pcoa_ait))[,c(1,3)]
scores.1.3

for (col in vars.envfit) {
  # Run the envfit function on the current column
  fit <- envfit(scores.1.3[!is.na(dataset.envfit_env[, col]), ],
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  # Extract the r2 and Pr(>r) values from the fit object
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  # Add the results to the data frame
  results_df.1.3 <- rbind(results_df.1.3, data.frame(Variable = col, r2 = r2, P = p_value, Dim1 = scores[1], Dim3 = scores[2]))
}

results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),]

#write.csv(results_df.1.3[order(-results_df.1.3$r2, results_df.1.3$P),], here("2_Dirichlet_ImmigrationSpeciation", "envfit_results_aitchison_pcoa_axes1_3_2026.csv"), row.names = F)


## Envfit PCoA2-3
# Loop over the numeric columns
scores.2.3 = as.data.frame(scores(pcoa_ait))[,c(2,3)]
scores.2.3

for (col in vars.envfit) {
  # Run the envfit function on the current column
  fit <- envfit(scores.2.3[!is.na(dataset.envfit_env[, col]), ],
                as.vector(dataset.envfit_env[!is.na(dataset.envfit_env[, col]), col]),
                permutations = 999)
  
  # Extract the r2 and Pr(>r) values from the fit object
  r2 <- fit$vectors$r
  p_value <- fit$vectors$pvals
  scores <- scores(fit, "vectors")
  
  # Add the results to the data frame
  results_df.2.3 <- rbind(results_df.2.3, data.frame(Variable = col, r2 = r2, P = p_value, Dim2 = scores[1], Dim3 = scores[2]))
}

results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),]

#write.csv(results_df.2.3[order(-results_df.2.3$r2, results_df.2.3$P),], here("2_Dirichlet_ImmigrationSpeciation", "envfit_results_aitchison_pcoa_axes2_3_2026.csv"), row.names = F)

##Make PCoA with metacommunities and envfit results as as arrows
# Load required libraries
library(ggplot2)
library(tibble) # for data manipulation
library(dplyr) # for data manipulation
library(compositions) # for CLR transformation

# Step 1: CLR transformation and Aitchison distance
otu_table_clr <- apply(otu_table(Dt), 1, function(x) clr(x + 1e-5)) # CLR with pseudocount
dist_aitchison <- dist(otu_table_clr, method = "euclidean") # Aitchison distance

# Step 2: Perform PCoA
pcoa_res <- cmdscale(dist_aitchison, k = 3, eig = TRUE)
scores_df <- as.data.frame(pcoa_res$points)
colnames(scores_df) <- c("PC1", "PC2", "PC3")
scores_df <- cbind(scores_df, sample_data(Dt)) # Combine with sample metadata
head(scores_df)

# Step 3: Fit environmental factors (if applicable)
envfit_res1.2 <- results_df.1.2
envfit_res1.3 <- results_df.1.3
envfit_res2.3 <- results_df.2.3

## Axes 1-2
# Extract significant environmental variables
significant_arrows <- envfit_res1.2[,4:5] %>%
  as.data.frame() %>%
  filter(envfit_res1.2$P < 0.05) # Only significant factors

significant_arrows = merge(significant_arrows, envfit_res1.2, by = c("Dim1","Dim2"))
head(significant_arrows)

# Step 4: Plot PCoA
names(scores_df)[names(scores_df) %in% "Dirichlet"] = "Metacommunity"
significant_arrows = significant_arrows[significant_arrows$Variable != "Dirichlet",]
significant_arrows$Variable = as.factor(significant_arrows$Variable)
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "Dought_days"] <- "Drought days"
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "SA"] <- "Surface area"
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "Depth.surface"] <- "Depth:Surface ratio"
significant_arrows

scores_df$Metacommunity = factor(scores_df$Metacommunity, levels = c(1:10))
scores_df$Location = factor(scores_df$Location, levels = c(1:10))
names(scores_df)[names(scores_df) == "Location"] = "Treehole"
scores_df = scores_df[!is.na(scores_df$Metacommunity),]
#scores_df = scores_df[,-ncol(scores_df)]

names(scores_df)[names(scores_df) == "X.SampleID"] = "ID.Sample"

custom_palette <- c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                    "#6D3F5B", "#1E213D", "#1C542D", "#763C28")

# Read in immigration and speciation rates
meta_treeholes = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv"), header = T)
meta_treeholes.merge = meta_treeholes[,names(meta_treeholes) %in% c("ID.Sample", "Immigrate_med", "Speciation_med")]
names(meta_treeholes.merge)[names(meta_treeholes.merge) == "Immigrate_med"] = "Immigration rate"
names(meta_treeholes.merge)[names(meta_treeholes.merge) == "Speciation_med"] = "Speciation rate"
head(meta_treeholes.merge)

scores_df.plot = merge(scores_df, meta_treeholes.merge, by = "ID.Sample")
length(unique(scores_df.plot$ID.Sample)) #164, correct
head(scores_df.plot)

# PCoA axes 1-2
ss1.2 = ggplot(data = scores_df.plot, aes(x = PC1 / 100, y = PC2 / 100, color = Metacommunity, size = `Immigration rate`, shape = Treehole)) +
  geom_point(alpha = 0.75, stroke = 1.5) +
  geom_segment(data = significant_arrows, 
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               arrow = arrow(length = unit(0.1, "cm")), 
               inherit.aes = FALSE, # Prevent it from inheriting `size = Volume`
               color = "black", linewidth = 1) + # Set arrow color to black
  geom_text(data = significant_arrows, 
            aes(x = Dim1 * 1.35, y = Dim2 * 1.35, label = Variable),
            inherit.aes = FALSE, # Prevent inheritance of `size`
            color = "black", size = 5, alpha = 1) + # Set label color to black
  scale_shape_manual(values = c(0:14, 15:18)) + # Assign additional shapes as needed
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right", # Move legend to the top
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2), # Set legend for color (Metacommunity) to 2 columns
    shape = guide_legend(ncol = 2), # Set legend for shape (Treehole) to 2 columns
    size = guide_legend(ncol = 1)   # Optional: Set legend for size to 2 columns
  ) +
  xlab(paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)")) +
  ylab(paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)"))

ss1.2

tiff(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_2_immigr.tiff"), units="cm", width=20, height=14, res=300)
ss1.2
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_2_immigr.pdf"), width=20/2.54, height=14/2.54)
ss1.2
dev.off()


ss1.2 = ggplot(data = scores_df.plot, aes(x = PC1 / 100, y = PC2 / 100, color = Metacommunity, size = `Speciation rate`, shape = Treehole)) +
  geom_point(alpha = 0.75, stroke = 1.5) +
  geom_segment(data = significant_arrows, 
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               arrow = arrow(length = unit(0.1, "cm")), 
               inherit.aes = FALSE, # Prevent it from inheriting `size = Volume`
               color = "black", linewidth = 1) + # Set arrow color to black
  geom_text(data = significant_arrows, 
            aes(x = Dim1 * 1.35, y = Dim2 * 1.35, label = Variable),
            inherit.aes = FALSE, # Prevent inheritance of `size`
            color = "black", size = 5, alpha = 1) + # Set label color to black
  scale_shape_manual(values = c(0:14, 15:18)) + # Assign additional shapes as needed
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right", # Move legend to the top
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2), # Set legend for color (Metacommunity) to 2 columns
    shape = guide_legend(ncol = 2), # Set legend for shape (Treehole) to 2 columns
    size = guide_legend(ncol = 1)   # Optional: Set legend for size to 2 columns
  ) +
  xlab(paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)")) +
  ylab(paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)"))

ss1.2

tiff(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_2_speciation.tiff"), units="cm", width=20, height=14, res=300)
ss1.2
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_2_speciation.pdf"), width=20/2.54, height=14/2.54)
ss1.2
dev.off()


## Axes 1-3
# Extract significant environmental variables
significant_arrows <- envfit_res1.3[,4:5] %>%
  as.data.frame() %>%
  filter(envfit_res1.3$P < 0.05) # Only significant factors

significant_arrows = merge(significant_arrows, envfit_res1.3, by = c("Dim1","Dim3"))
head(significant_arrows)

# PCoA axes 1-3
ss1.3 = ggplot(data = scores_df.plot, aes(x = PC1 / 100, y = PC3 / 100, color = Metacommunity, size = `Immigration rate`, shape = Treehole)) +
  geom_point(alpha = 0.75, stroke = 1.5) +
  geom_segment(data = significant_arrows, 
               aes(x = 0, y = 0, xend = Dim1, yend = Dim3), 
               arrow = arrow(length = unit(0.1, "cm")), 
               inherit.aes = FALSE, # Prevent it from inheriting `size = Volume`
               color = "black", linewidth = 1) + # Set arrow color to black
  geom_text(data = significant_arrows, 
            aes(x = Dim1 * 1.35, y = Dim3 * 1.35, label = Variable),
            inherit.aes = FALSE, # Prevent inheritance of `size`
            color = "black", size = 5, alpha = 1) + # Set label color to black
  scale_shape_manual(values = c(0:14, 15:18)) + # Assign additional shapes as needed
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right", # Move legend to the top
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2), # Set legend for color (Metacommunity) to 2 columns
    shape = guide_legend(ncol = 2), # Set legend for shape (Treehole) to 2 columns
    size = guide_legend(ncol = 1)   # Optional: Set legend for size to 2 columns
  ) +
  xlab(paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)")) +
  ylab(paste0("PCoA3 (", round(pcoa_res$eig[3] / sum(pcoa_res$eig) * 100, 1), "%)"))

ss1.3

tiff(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_3_immigr.tiff"), units="cm", width=20, height=14, res=300)
ss1.3
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_3_immigr.pdf"), width=20/2.54, height=14/2.54)
ss1.3
dev.off()


ss1.3 = ggplot(data = scores_df.plot, aes(x = PC1 / 100, y = PC3 / 100, color = Metacommunity, size = `Speciation rate`, shape = Treehole)) +
  geom_point(alpha = 0.75, stroke = 1.5) +
  geom_segment(data = significant_arrows, 
               aes(x = 0, y = 0, xend = Dim1, yend = Dim3), 
               arrow = arrow(length = unit(0.1, "cm")), 
               inherit.aes = FALSE, # Prevent it from inheriting `size = Volume`
               color = "black", linewidth = 1) + # Set arrow color to black
  geom_text(data = significant_arrows, 
            aes(x = Dim1 * 1.35, y = Dim3 * 1.35, label = Variable),
            inherit.aes = FALSE, # Prevent inheritance of `size`
            color = "black", size = 5, alpha = 1) + # Set label color to black
  scale_shape_manual(values = c(0:14, 15:18)) + # Assign additional shapes as needed
  scale_color_manual(values = custom_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "right", # Move legend to the top
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", colour = "black", vjust = 2.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2), # Set legend for color (Metacommunity) to 2 columns
    shape = guide_legend(ncol = 2), # Set legend for shape (Treehole) to 2 columns
    size = guide_legend(ncol = 1)   # Optional: Set legend for size to 2 columns
  ) +
  xlab(paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)")) +
  ylab(paste0("PCoA3 (", round(pcoa_res$eig[3] / sum(pcoa_res$eig) * 100, 1), "%)"))

ss1.3

tiff(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_3_speciation.tiff"), units="cm", width=20, height=14, res=300)
ss1.3
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "PCOA_treehole_metacommunity_k10_axes1_3_speciation.pdf"), width=20/2.54, height=14/2.54)
ss1.3
dev.off()


## Speciation rate associations with environmental variables

names(meta_treeholes)[names(meta_treeholes) == "Location"] = "Treehole"
names(meta_treeholes)[names(meta_treeholes) == "Dirichlet"] = "Metacommunity"
meta_treeholes$Metacommunity = factor(meta_treeholes$Metacommunity, levels = 1:10)
meta_treeholes$Treehole = factor(meta_treeholes$Treehole, levels = 1:10)
meta_treeholes$Depth.surface.ratio = meta_treeholes$Depth / meta_treeholes$SA
head(meta_treeholes)

# Define environmental variables that are specific to each sample and can thus be associated with immigration rates
env.vars = c("Moisture", "Max.Lum", "Air.Temp", "Grass.Temp", "Soil.Temp2", "Soil.Temp4", "Air.Temp.max",
             "Air.Temp.min", "Grass.Temp.max", "Grass.Temp.min",  "Rainfall", "Dought_days", "Depth.surface.ratio")

# Multiple regression with metacommunity as factor
mult.reg <- as.formula(
  paste("Speciation_med ~", paste(env.vars, collapse = " + "))
)
model <- lm(mult.reg, data = meta_treeholes)
summary(model)

# ---------------------------
# Fit model
# ---------------------------
model <- lm(Speciation_med ~ log(Depth.surface.ratio), data = meta_treeholes)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R2 and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·log(x)\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes, aes(x = Depth.surface.ratio, y = Speciation_med)) +
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                                "#6D3F5B", "#1E213D", "#1C542D", "#763C28")) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    se = TRUE,
    color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1, vjust = 1,
    size = 4
  ) +
  
  labs(
    x = "Depth:Surface ratio",
    y = "Speciation rate"
  ) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")
  ) + xlim(0, 0.3) + ylim(450,700) +
  scale_shape_manual(values = c(0:14, 15:18)) +
  guides(shape = guide_legend(ncol = 2, byrow = FALSE), color = guide_legend(ncol = 2, byrow = FALSE))
p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_depth_surface_ratio_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_depth_surface_ratio_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


## Same with drought days
max_drought <- meta_treeholes %>%
  group_by(Treehole) %>%
  summarise(
    max_drought_days = max(Dought_days, na.rm = TRUE),
    .groups = "drop"
  )
max_drought

meta_treeholes_drought = merge(meta_treeholes, max_drought, by = "Treehole")
head(meta_treeholes_drought)

# ---------------------------
# Fit model
# ---------------------------
model <- lm(Speciation_med ~ max_drought_days,
            data = meta_treeholes_drought)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R² and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·x\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes_drought,
            aes(x = max_drought_days, y = Speciation_med)) +
  
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3",
                                "#497E76", "#D6AE01", "#6D3F5B", "#1E213D",
                                "#1C542D", "#763C28")) +
  
  geom_smooth(method = "lm",
              formula = y ~ x,   # ✅ linear instead of log
              se = TRUE,
              color = "black") +
  
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1, vjust = 1,
           size = 4) +
  
  labs(
    x = "Drought prevalence (total drought days)",
    y = "Speciation rate"
  ) +
  
  ylim(c(450, 700)) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")) +
  
  scale_shape_manual(values = c(0:14, 15:18)) +
  
  guides(shape = guide_legend(ncol = 2, byrow = FALSE),
         color = guide_legend(ncol = 2, byrow = FALSE))

p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


# Same with volume
# ---------------------------
# Fit model (log)
# ---------------------------
model <- lm(Speciation_med ~ log(Volume),
            data = meta_treeholes_drought)

# Extract coefficients
coef_vals <- coef(model)
intercept <- round(coef_vals[1], 2)
slope <- round(coef_vals[2], 2)

# Extract R² and p-value
model_summary <- summary(model)
r2 <- round(model_summary$r.squared, 3)
p_val <- signif(model_summary$coefficients[2, 4], 3)

# Create label
label_text <- paste0(
  "y = ", intercept, " + ", slope, "·log(x)\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(meta_treeholes_drought,
            aes(x = Volume, y = Speciation_med)) +
  
  geom_point(aes(color = Metacommunity, shape = Treehole), size = 3) +
  
  scale_color_manual(values = c("#C35831", "#9B111E", "#7E7B52", "#2271B3",
                                "#497E76", "#D6AE01", "#6D3F5B", "#1E213D",
                                "#1C542D", "#763C28")) +
  
  geom_smooth(method = "lm",
              formula = y ~ log(x),
              se = TRUE,
              color = "black") +
  
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1, vjust = 1,
           size = 4) +
  
  labs(
    x = expression(Volume~(cm^3)),
    y = "Speciation rate"
  ) +
  
  ylim(c(450, 700)) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        plot.title = element_text(face = "bold")) +
  
  scale_shape_manual(values = c(0:14, 15:18)) +
  
  guides(shape = guide_legend(ncol = 2, byrow = FALSE),
         color = guide_legend(ncol = 2, byrow = FALSE))

p

tiff(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_volume_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("2_Dirichlet_ImmigrationSpeciation", "Speciation_volume_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()
