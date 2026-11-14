###################################
#  rda_functions.R
##################################
# Helper functions shared by rda_composition_environment.R and
# rda_neutral_rates_environment.R: data loading, CLR transform, explained
# variance of a linear model on a multivariate response, and a permutation
# test that respects the tree-hole structure of the data.
#
# Tree-hole morphology (Depth, SA, Volume, orientation) is constant within each
# of the 10 tree holes, so for these variables the replication unit is the tree
# hole, not the sample. Permuting samples freely treats ~160 repeated samples
# as independent and gives p-values that are far too small; morphology_subsets()
# instead shuffles the morphology values among tree holes.

library(phyloseq)
library(vegan)
library(here)

# --- Load the K = 10 data exactly as in taxonomy_betadiversity_envfit.R ---------
# Returns a list with the (unrarefied) count matrix (samples x ASVs, mitochondria
# and chloroplasts removed, CONTROL and samples without metacommunity dropped) and
# the sample metadata. The CLR is scale invariant, so no rarefaction is needed.
load_treehole_counts = function() {
  meta = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"),
                    header = T, sep = "\t", check.names = F)
  meta = meta[, names(meta) != ""] # trailing tab in the header
  names(meta)[names(meta) == "ID.Sample"] = "#SampleID"
  mapfile = tempfile(fileext = ".tsv")
  write.table(meta, mapfile, sep = "\t", quote = F, row.names = F)
  otufile = here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv")
  physeq = import_qiime(otufilename = otufile, mapfilename = mapfile)

  Dt = subset_taxa(physeq, Family != "Mitochondria")
  Dt = subset_taxa(Dt, Order != "Chloroplast")
  Dt = subset_samples(Dt, Sample.type == "THS" & Location != "CONTROL" & !is.na(Dirichlet))
  Dt = prune_taxa(taxa_sums(Dt) > 0, Dt)

  counts = as(otu_table(Dt), "matrix")
  if (taxa_are_rows(Dt)) counts = t(counts)
  env = data.frame(sample_data(Dt), check.names = F)
  names(env)[names(env) == "X.SampleID"] = "ID.Sample"
  env$ID.Sample = rownames(env)
  list(counts = counts, env = env)
}

# --- Centred log-ratio transform within each sample (rows) ------------------------
# Euclidean distances between CLR profiles are Aitchison distances, so an RDA on
# the CLR matrix is the constrained counterpart of the Aitchison PCoA. Same
# transformation as vegdist(counts, "aitchison", pseudocount = 1) in
# taxonomy_betadiversity_envfit.R.
clr_rows = function(counts, pseudocount = 1) {
  logx = log(counts + pseudocount)
  logx - rowMeans(logx)
}

# --- Environmental predictors -------------------------------------------------
# Tree-hole morphology: constant per tree hole. Size variables are log-transformed
# (right-skewed, spanning an order of magnitude); the aspect in degrees is circular,
# so it is split into northness (cos) and eastness (sin).
# Size and Shape: orthogonal summary of log Depth and log SA. Both are
# standardised across the 10 tree holes (one value per tree hole), then
#   Size  = (z.log.Depth + z.log.SA) / sqrt(2)   ~ log Volume
#   Shape = (z.log.Depth - z.log.SA) / sqrt(2)   ~ log Depth:SA
# These are the two principal components of (log Depth, log SA) and are exactly
# uncorrelated across tree holes. They span the same space as Depth + SA, so a
# model with Size + Shape fits exactly as well, but each coefficient can be read
# on its own.
# Time-varying variables: change from sample to sample within a tree hole.
add_predictors = function(env) {
  env$Treehole = factor(env$Location, levels = 1:10)
  env$Metacommunity = factor(env$Dirichlet, levels = 1:10)
  # Deg.North changes for a couple of samples of one tree hole: use the tree-hole mode
  mode_th = tapply(env$Deg.North, env$Treehole, function(z) as.numeric(names(which.max(table(z)))))
  aspect = mode_th[as.character(env$Treehole)] * pi / 180
  env$Northness = cos(aspect)
  env$Eastness = sin(aspect)
  env$log.Depth = log(env$Depth)
  env$log.SA = log(env$SA)
  env$log.Volume = log(env$Volume)
  env$log.Depth.surface = log(env$Depth / env$SA)
  env$log.Max.Lum = log(env$Max.Lum)
  env$Drought_days = env$Dought_days # the metadata column name is misspelt
  th = env[!duplicated(env$Treehole), c("Treehole", "log.Depth", "log.SA")]
  zD = as.vector(scale(th$log.Depth)); zS = as.vector(scale(th$log.SA))
  i = match(env$Treehole, th$Treehole)
  env$Size = ((zD + zS) / sqrt(2))[i]
  env$Shape = ((zD - zS) / sqrt(2))[i]
  env
}

vars.morphology = c("Size", "Shape", "log.Depth", "log.SA", "log.Volume", "log.Depth.surface", "Northness", "Eastness")
# In the metadata Volume = SA * Depth / 3 exactly (a cone), so in log scale
#   log.Volume = log.SA + log.Depth - log(3)   and   log.Depth.surface = log.Depth - log.SA.
# The six size variables span only two dimensions: any two of them give the same
# model as Size + Shape, and three or more are exactly collinear.
vars.size = c("Size", "Shape", "log.Depth", "log.SA", "log.Volume", "log.Depth.surface")
vars.size.basis = c("Size", "Shape") # used in full models (orthogonal)
vars.temporal = c("Moisture", "log.Max.Lum", "Air.Temp", "Air.Temp.max", "Air.Temp.min",
                  "Soil.Temp4", "Rainfall", "Drought_days")

# --- Drop variables with variance inflation factor above a threshold --------------
prune_vif = function(df, vars, threshold = 10) {
  repeat {
    v = diag(solve(cor(df[, vars]), tol = 1e-12))
    names(v) = vars
    if (max(v) <= threshold) return(list(vars = vars, vif = v))
    vars = vars[-which.max(v)]
  }
}

# --- Fraction of variance of Y explained by the predictors in X -----------------
# Y: response matrix (samples x responses). Returns R2 and adjusted R2, identical
# to vegan::RsquareAdj(rda(Y ~ X)), but fast enough to be computed thousands of
# times inside a permutation test (works with the n x n Gram matrix of Y).
gram = function(Y) {
  Yc = scale(Y, center = TRUE, scale = FALSE)
  tcrossprod(Yc)
}
r2_gram = function(G, X) {
  Xc = scale(as.matrix(X), center = TRUE, scale = FALSE)
  Q = qr.Q(qr(Xc))
  r2 = sum(diag(crossprod(Q, G %*% Q))) / sum(diag(G))
  n = nrow(Xc); k = ncol(Xc)
  c(R2 = r2, adjR2 = 1 - (1 - r2) * (n - 1) / (n - k - 1))
}

# --- Morphology combinations, permutation test at the tree-hole level -------------
# Every combination of 1 to max.k morphology variables (at most two size variables,
# and then only Size + Shape, since any other pair gives the same fit).
#
# Permutations shuffle the rows of the tree-hole table (one row per tree hole) and
# map the shuffled morphology back to the samples: each sample keeps its
# composition and its tree-hole membership but gets the morphology of another tree
# hole. All morphology variables move together, so their correlations are kept.
#
# Two p-values per combination:
#  * P.treehole.perm: the combination tested on its own (R2 of the combination
#    under its own permutation distribution).
#  * P.selection.corrected: the observed adjusted R2 compared with the maximum
#    adjusted R2 over ALL combinations in each permutation (max-statistic, as in
#    Westfall & Young 1993). This accounts for having tried every combination and
#    kept the best one, so it is the valid p-value for the selected model.
morphology_subsets = function(G, env, vars = vars.morphology, max.k = 3, nperm = 9999, seed = 1) {
  subsets = list()
  for (k in 1:max.k) {
    for (cmb in combn(vars, k, simplify = FALSE)) {
      n.size = sum(cmb %in% vars.size)
      if (n.size > 2) next # exactly collinear
      if (n.size == 2 && !all(vars.size.basis %in% cmb)) next # same fit as Size + Shape
      subsets[[length(subsets) + 1]] = cmb
    }
  }
  th = as.character(env$Treehole)
  th.tab = env[!duplicated(th), vars, drop = FALSE]
  rownames(th.tab) = th[!duplicated(th)]
  fit_all = function(tab) sapply(subsets, function(cmb) r2_gram(G, tab[th, cmb, drop = FALSE]))
  obs = fit_all(th.tab) # 2 x n.subsets (R2, adjR2)
  set.seed(seed)
  null = replicate(nperm, {
    shuffled = th.tab[sample(nrow(th.tab)), , drop = FALSE]
    rownames(shuffled) = rownames(th.tab)
    fit_all(shuffled)["adjR2", ]
  }) # n.subsets x nperm
  p.own = (rowSums(null >= obs["adjR2", ] - 1e-12) + 1) / (nperm + 1)
  max.null = apply(null, 2, max)
  p.sel = sapply(obs["adjR2", ], function(o) (sum(max.null >= o - 1e-12) + 1) / (nperm + 1))
  res = data.frame(Variables = sapply(subsets, paste, collapse = " + "), k = lengths(subsets),
                   R2 = obs["R2", ], adjR2 = obs["adjR2", ],
                   P.treehole.perm = p.own, P.selection.corrected = p.sel)
  res[order(-res$adjR2), ]
}

# --- Biplot of an RDA (sites coloured by metacommunity, arrows = predictors) ---------
custom_palette = c("#C35831", "#9B111E", "#7E7B52", "#2271B3", "#497E76", "#D6AE01",
                   "#6D3F5B", "#1E213D", "#1C542D", "#763C28")

rda_biplot = function(mod, env, axes = c(1, 2), title = NULL) {
  library(ggplot2)
  sc = scores(mod, display = c("sites", "bp"), choices = axes, scaling = 2)
  sites = data.frame(sc$sites, Metacommunity = env$Metacommunity, Treehole = env$Treehole)
  bp = data.frame(sc$biplot, Variable = rownames(sc$biplot))
  mult = 0.8 * max(abs(sites[, 1:2])) / max(abs(bp[, 1:2]))
  bp[, 1:2] = bp[, 1:2] * mult
  eig = mod$CCA$eig / mod$tot.chi * 100
  ax = names(sites)[1:2]
  ggplot(sites, aes(x = .data[[ax[1]]], y = .data[[ax[2]]])) +
    geom_point(aes(color = Metacommunity, shape = Treehole), size = 2.5, alpha = 0.8, stroke = 1) +
    geom_segment(data = bp, aes(x = 0, y = 0, xend = .data[[ax[1]]], yend = .data[[ax[2]]]),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black") +
    geom_text(data = bp, aes(x = .data[[ax[1]]] * 1.12, y = .data[[ax[2]]] * 1.12, label = Variable),
              size = 3.5) +
    scale_color_manual(values = custom_palette) +
    scale_shape_manual(values = c(0:9)) +
    labs(x = sprintf("%s (%.1f%% of total)", ax[1], eig[axes[1]]),
         y = sprintf("%s (%.1f%% of total)", ax[2], eig[axes[2]]), title = title) +
    theme_bw() +
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
}

save_pdf = function(p, file, width = 20, height = 14) {
  pdf(file, width = width / 2.54, height = height / 2.54)
  print(p)
  dev.off()
}

# --- Results summary: one row per test, shared by both RDA scripts -------------------
# Each script writes <its folder>/results_summary.csv; write_combined_summary()
# stacks the summaries found in 4.1_Niche_Envir/rda/*/ into
# 4.1_Niche_Envir/rda/results_summary.csv.
summary_row = function(Step, Analysis, Response, Predictors, Permutation, adjR2 = NA,
                       Statistic = NA, P = NA, P.corrected = NA, File = NA) {
  data.frame(Step = Step, Analysis = Analysis, Response = Response, Predictors = Predictors,
             Permutation = Permutation, adjR2 = round(adjR2, 3), Statistic = Statistic,
             P = signif(P, 3), P.corrected = signif(P.corrected, 3), File = File)
}

write_combined_summary = function(dirRda = here("4.1_Niche_Envir", "rda")) {
  files = file.path(dirRda, c("composition", "neutral_rates"), "results_summary.csv")
  files = files[file.exists(files)]
  tab = do.call(rbind, lapply(files, function(f)
    data.frame(Folder = basename(dirname(f)), read.csv(f, check.names = FALSE))))
  write.csv(tab, file.path(dirRda, "results_summary.csv"), row.names = FALSE, na = "")
}

# --- Explained variance per analysis (bar chart, one panel per response) -------------
plot_explained = function(tab, title = NULL) {
  library(ggplot2)
  tab$Analysis = factor(tab$Analysis, levels = rev(unique(tab$Analysis)))
  ggplot(tab, aes(x = adjR2, y = Analysis)) +
    geom_col(width = 0.6, fill = "#2271B3") +
    geom_text(aes(label = sprintf("%.2f", adjR2)), hjust = -0.15, size = 3.2) +
    facet_wrap(~ Response, nrow = 1) +
    scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(x = "Adjusted R2", y = NULL, title = title) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
}

# Same data as plot_explained, with the responses as side-by-side bars in one panel
plot_explained_grouped = function(tab, title = NULL, colors = c("#2271B3", "#D6AE01")) {
  library(ggplot2)
  tab$Analysis = factor(tab$Analysis, levels = rev(unique(tab$Analysis)))
  tab$Response = factor(tab$Response, levels = rev(unique(tab$Response)))
  dodge = position_dodge(width = 0.75, preserve = "single")
  ggplot(tab, aes(x = adjR2, y = Analysis, fill = Response)) +
    geom_col(width = 0.7, position = dodge) +
    geom_text(aes(label = sprintf("%.2f", adjR2)), position = dodge, hjust = -0.15, size = 3) +
    scale_fill_manual(values = setNames(colors, rev(levels(tab$Response))), breaks = rev(levels(tab$Response))) +
    scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(x = "Adjusted R2", y = NULL, fill = NULL, title = title) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
}
