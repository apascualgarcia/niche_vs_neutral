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

tiff(here("3.1_Neutral_NMGS", "figures", "Immigration_depth_surface_ratio_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Immigration_depth_surface_ratio_metacommunity.pdf"),
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

tiff(here("3.1_Neutral_NMGS", "figures", "Immigration_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Immigration_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
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

tiff(here("3.1_Neutral_NMGS", "figures", "Immigration_volume_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Immigration_volume_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()


# Depth:Surface association with drought prevalence
max_drought_depth <- meta_treeholes_drought %>%
  dplyr::distinct(Treehole, Depth.surface.ratio, max_drought_days)
max_drought_depth

model <- lm(max_drought_days ~ Depth.surface.ratio,
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
  "y = ", intercept,
  ifelse(slope >= 0, " + ", " - "),
  abs(slope), "·x\n",
  "R² = ", r2, "\n",
  "p = ", p_val
)

p <- ggplot(meta_treeholes_drought,
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
           label = label_text,
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

tiff(here("3.1_Neutral_NMGS", "figures", "Depth_surface_ratio_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Depth_surface_ratio_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
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
otufile = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv")

physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)
physeq

##RAREFACTION
#Rarefaction plot
#source("Rfunctions.R")
#rarefaction.plot = ggrare(physeq, step = 200, se = FALSE) + xlim(c(0,max(sample_sums(physeq)))) + ylim(0,2000) + theme_bw() + theme(axis.text.x = element_text(size = 20), axis.text.y=element_text(size=20), legend.position = "null",
#                                                                                                                                    axis.title.x = element_text(size = 20, margin = margin(r = 15)), axis.title.y = element_text(size =20, margin = margin(r = 10))) + xlab("Number of reads")

#tiff(here("3.1_Neutral_NMGS", "figures", "Rarefaction_plot_treeholes.tiff"), units="cm", width=15, height=15, res=300)
#rarefaction.plot
#dev.off()

library(vegan)
library(dplyr) # for data manipulation

# The PCoA scores and envfit results are those of taxonomy_betadiversity_envfit.R
# (Aitchison distance: CLR within each sample of the unrarefied counts, without
# mitochondria/chloroplasts), so both scripts show the same ordination.
dirEnvfit = here("4.1_Niche_Envir", "figures")
pct_exp = read.csv(file.path(dirEnvfit, "pcoa_aitchison_variance.csv"))$Percent.variance
envfit_res1.2 = read.csv(file.path(dirEnvfit, "envfit_results_aitchison_pcoa_axes1_2_2026.csv"))
envfit_res1.3 = read.csv(file.path(dirEnvfit, "envfit_results_aitchison_pcoa_axes1_3_2026.csv"))
envfit_res2.3 = read.csv(file.path(dirEnvfit, "envfit_results_aitchison_pcoa_axes2_3_2026.csv"))

scores_df = read.csv(file.path(dirEnvfit, "pcoa_aitchison_scores.csv"), row.names = 1)
colnames(scores_df) = c("PC1", "PC2", "PC3")
scores_df = cbind(scores_df, data.frame(sample_data(physeq))[rownames(scores_df), ]) # Combine with sample metadata

# envfit arrows have length sqrt(r2) (<= 1): stretch them to the range of the scores
arrow.mult = 0.8 * max(abs(scores_df[, c("PC1", "PC2", "PC3")])) /
  max(sqrt(c(envfit_res1.2$r2, envfit_res1.3$r2, envfit_res2.3$r2)))
envfit_res1.2[, 4:5] = envfit_res1.2[, 4:5] * arrow.mult
envfit_res1.3[, 4:5] = envfit_res1.3[, 4:5] * arrow.mult
envfit_res2.3[, 4:5] = envfit_res2.3[, 4:5] * arrow.mult


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
ss1.2 = ggplot(data = scores_df.plot, aes(x = PC1, y = PC2, color = Metacommunity, size = `Immigration rate`, shape = Treehole)) +
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
  xlab(paste0("PCoA1 (", round(pct_exp[1], 1), "%)")) +
  ylab(paste0("PCoA2 (", round(pct_exp[2], 1), "%)"))

ss1.2

tiff(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_2_immigr.tiff"), units="cm", width=20, height=14, res=300)
ss1.2
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_2_immigr.pdf"), width=20/2.54, height=14/2.54)
ss1.2
dev.off()


ss1.2 = ggplot(data = scores_df.plot, aes(x = PC1, y = PC2, color = Metacommunity, size = `Speciation rate`, shape = Treehole)) +
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
  xlab(paste0("PCoA1 (", round(pct_exp[1], 1), "%)")) +
  ylab(paste0("PCoA2 (", round(pct_exp[2], 1), "%)"))

ss1.2

tiff(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_2_speciation.tiff"), units="cm", width=20, height=14, res=300)
ss1.2
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_2_speciation.pdf"), width=20/2.54, height=14/2.54)
ss1.2
dev.off()


## Axes 1-3
# Extract significant environmental variables
significant_arrows <- envfit_res1.3[,4:5] %>%
  as.data.frame() %>%
  filter(envfit_res1.3$P < 0.05) # Only significant factors

significant_arrows = merge(significant_arrows, envfit_res1.3, by = c("Dim1","Dim3"))
head(significant_arrows)
significant_arrows = significant_arrows[significant_arrows$Variable != "Dirichlet",]
significant_arrows$Variable = as.factor(significant_arrows$Variable)
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "Dought_days"] <- "Drought days"
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "SA"] <- "Surface area"
levels(significant_arrows$Variable)[levels(significant_arrows$Variable) == "Depth.surface"] <- "Depth:Surface ratio"

# PCoA axes 1-3
ss1.3 = ggplot(data = scores_df.plot, aes(x = PC1, y = PC3, color = Metacommunity, size = `Immigration rate`, shape = Treehole)) +
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
  xlab(paste0("PCoA1 (", round(pct_exp[1], 1), "%)")) +
  ylab(paste0("PCoA3 (", round(pct_exp[3], 1), "%)"))

ss1.3

tiff(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_3_immigr.tiff"), units="cm", width=20, height=14, res=300)
ss1.3
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_3_immigr.pdf"), width=20/2.54, height=14/2.54)
ss1.3
dev.off()


ss1.3 = ggplot(data = scores_df.plot, aes(x = PC1, y = PC3, color = Metacommunity, size = `Speciation rate`, shape = Treehole)) +
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
  xlab(paste0("PCoA1 (", round(pct_exp[1], 1), "%)")) +
  ylab(paste0("PCoA3 (", round(pct_exp[3], 1), "%)"))

ss1.3

tiff(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_3_speciation.tiff"), units="cm", width=20, height=14, res=300)
ss1.3
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "PCOA_treehole_metacommunity_k10_axes1_3_speciation.pdf"), width=20/2.54, height=14/2.54)
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

tiff(here("3.1_Neutral_NMGS", "figures", "Speciation_depth_surface_ratio_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Speciation_depth_surface_ratio_metacommunity.pdf"), width=14/2.54, height=12/2.54)
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

tiff(here("3.1_Neutral_NMGS", "figures", "Speciation_drought_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Speciation_drought_metacommunity.pdf"), width=14/2.54, height=12/2.54)
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

tiff(here("3.1_Neutral_NMGS", "figures", "Speciation_volume_metacommunity.tiff"), units="cm", width=14, height=12, res=300)
p
dev.off()

pdf(here("3.1_Neutral_NMGS", "figures", "Speciation_volume_metacommunity.pdf"), width=14/2.54, height=12/2.54)
p
dev.off()
