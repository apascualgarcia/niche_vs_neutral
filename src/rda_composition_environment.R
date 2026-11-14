###################################
#  rda_composition_environment.R
##################################
# Which combination of environmental variables best explains community
# composition? Complements the envfit analysis of taxonomy_betadiversity_envfit.R,
# where depth, surface area and volume were all correlated with PCoA1 and with
# each other. Here composition is the response and the environment the
# explanatory variables (redundancy analysis). Same data and transformation as the
# envfit: unrarefied counts without mitochondria/chloroplasts, CLR within each
# sample (pseudocount 1), i.e. Aitchison geometry. Two responses: the whole CLR
# matrix, and the first three axes of the Aitchison PCoA (the envfit axes).
#
#  1. Collinearity among predictors (correlations, VIF). Depth and SA enter the
#     models as the orthogonal pair Size + Shape (see rda_functions.R).
#  2. Forward selection (adjusted R2, Blanchet et al. 2008). Exploratory: these
#     permutation tests treat samples as independent.
#  3. Tree-hole morphology: every combination of 1-3 variables, tested by
#     permuting morphology among the 10 tree holes (the true replication unit),
#     with p-values corrected for having selected the best combination.
#     Downstream models use the a priori pair Size + Shape (the two dimensions of
#     tree-hole morphology), not the data-driven best combination.
#  4. Time-varying variables: partial RDA conditioned on tree hole, with
#     permutations restricted within tree holes.
#  5. Variation partitioning: morphology vs time-varying variables, with the
#     tree-hole identity as the upper bound of what morphology can explain.
#  6. Final model (best morphology + significant time-varying variables).
#
# Input:  1_dada2/otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv
#         1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
# Output: 4.1_Niche_Envir/rda/composition/  (results_summary.csv, tables .csv,
#         figures .pdf, summary.txt) and 4.1_Niche_Envir/rda/results_summary.csv
###################################

rm(list = ls())
library(here)
library(vegan)
library(permute)
library(ggplot2)
source(here("src", "rda_functions.R"))

dirOut = here("4.1_Niche_Envir", "rda", "composition")
dir.create(dirOut, recursive = TRUE, showWarnings = FALSE)
NPERM = 999       # permutations for the vegan tests
NPERM.TH = 9999   # permutations for the tree-hole level tests (fast)

# READ DATA -----------------------------------------------------------------
dat = load_treehole_counts()
env = add_predictors(dat$env)
Y = clr_rows(dat$counts)

# Aitchison PCoA = PCA of the CLR matrix, on all samples as in the envfit
pca = prcomp(Y, center = TRUE, scale. = FALSE)
pcoa.var = pca$sdev^2 / sum(pca$sdev^2) * 100
cat("Aitchison PCoA, % variance of axes 1-3:", round(pcoa.var[1:3], 1), "\n")
Y3 = pca$x[, 1:3]
colnames(Y3) = c("PCoA1", "PCoA2", "PCoA3")

vars.all = c(vars.morphology, vars.temporal)
keep = complete.cases(env[, vars.all])
cat("Samples:", nrow(env), "; complete cases for all predictors:", sum(keep), "\n")
env = env[keep, ]
Y = Y[keep, ]
Y = Y[, apply(Y, 2, var) > 0]
Y3 = Y3[keep, ]
G = gram(Y)
G3 = gram(Y3)
responses = c(CLR = "CLR community (all ASVs)", PCoA = "Aitchison PCoA1-3")

# 1. COLLINEARITY ---------------------------------------------------------------
cor.env = cor(env[, vars.all])
write.csv(round(cor.env, 3), file.path(dirOut, "predictor_correlations.csv"))
th.tab = env[!duplicated(env$Treehole), ]
cor.size = data.frame(Pair = c("log.Depth vs log.SA", "log.Volume vs log.Depth.surface", "Size vs Shape"),
                      r.treeholes = round(c(cor(th.tab$log.Depth, th.tab$log.SA),
                                            cor(th.tab$log.Volume, th.tab$log.Depth.surface),
                                            cor(th.tab$Size, th.tab$Shape)), 3),
                      r.samples = round(c(cor(env$log.Depth, env$log.SA),
                                          cor(env$log.Volume, env$log.Depth.surface),
                                          cor(env$Size, env$Shape)), 3))
write.csv(cor.size, file.path(dirOut, "size_shape_correlations.csv"), row.names = F)
print(cor.size)

# Figure: the 10 tree holes in the (log Depth, log SA) plane, standardised, with
# Size and Shape as the two diagonals of that plane
th.z = data.frame(Treehole = th.tab$Treehole,
                  zD = as.vector(scale(th.tab$log.Depth)), zS = as.vector(scale(th.tab$log.SA)))
ax.len = 1.1 * max(abs(c(th.z$zD, th.z$zS)))
axes.ss = data.frame(Axis = c("Size", "Shape"), x = c(ax.len, ax.len), y = c(ax.len, -ax.len))
p.ss = ggplot2::ggplot(th.z, ggplot2::aes(zD, zS)) +
  ggplot2::geom_segment(data = axes.ss, ggplot2::aes(x = -x, y = -y, xend = x, yend = y),
                        arrow = grid::arrow(length = grid::unit(0.2, "cm")), color = "#2271B3") +
  ggplot2::geom_text(data = axes.ss, ggplot2::aes(x = x, y = y, label = Axis),
                     hjust = -0.15, color = "#2271B3", size = 4) +
  ggplot2::geom_point(size = 2.5) +
  ggrepel::geom_text_repel(ggplot2::aes(label = paste0("TH", Treehole)), size = 3.2, seed = 1) +
  ggplot2::coord_equal(xlim = c(-ax.len, ax.len * 1.35), ylim = c(-ax.len, ax.len)) +
  ggplot2::labs(x = "log Depth (standardised across tree holes)",
                y = "log surface area (standardised across tree holes)",
                title = sprintf("r(log Depth, log SA) = %.2f; r(Size, Shape) = 0", cor.size$r.treeholes[1])) +
  ggplot2::theme_bw()
save_pdf(p.ss, file.path(dirOut, "size_shape_definition.pdf"), 14, 13)

# Figure: Size and Shape against the original variables. With mu and sigma the mean
# and SD of log Depth (D) and log SA (S) across tree holes, and V = SA * Depth / 3,
#   log Volume    = (mu.D + mu.S - log 3) + a * Size + b * Shape
#   log Depth:SA  = (mu.D - mu.S)         + b * Size + a * Shape
# with a = (sigma.D + sigma.S) / sqrt(2) and b = (sigma.D - sigma.S) / sqrt(2).
# b is small because sigma.D ~ sigma.S, so Size ~ log Volume and Shape ~ log Depth:SA.
sD = sd(th.tab$log.Depth); sS = sd(th.tab$log.SA)
coef.ss = data.frame(Variable = c("log.Volume", "log.Depth.surface"),
                     Intercept = c(mean(th.tab$log.Depth) + mean(th.tab$log.SA) - log(3), mean(th.tab$log.Depth) - mean(th.tab$log.SA)),
                     Size = c(sD + sS, sD - sS) / sqrt(2), Shape = c(sD - sS, sD + sS) / sqrt(2),
                     r.Size = c(cor(th.tab$Size, th.tab$log.Volume), cor(th.tab$Size, th.tab$log.Depth.surface)),
                     r.Shape = c(cor(th.tab$Shape, th.tab$log.Volume), cor(th.tab$Shape, th.tab$log.Depth.surface)))
write.csv(coef.ss, file.path(dirOut, "size_shape_vs_original.csv"), row.names = F)
old.lab = c(log.Volume = "Volume (cm3, log scale)", log.Depth.surface = "Depth:surface ratio (log scale)")
ss.long = do.call(rbind, lapply(names(old.lab), function(o) do.call(rbind, lapply(c("Size", "Shape"), function(n)
  data.frame(Original = old.lab[o], New = n, x = th.tab[[n]], y = exp(th.tab[[o]]), Treehole = th.tab$Treehole,
             r = cor(th.tab[[n]], th.tab[[o]]))))))
ss.long$Original = factor(ss.long$Original, levels = old.lab)
ss.long$New = factor(ss.long$New, levels = c("Size", "Shape"))
ss.r = ss.long[!duplicated(ss.long[, c("Original", "New")]), ]
p.ss2 = ggplot2::ggplot(ss.long, ggplot2::aes(x, y)) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey50", linewidth = 0.5) +
  ggplot2::geom_point(size = 2) +
  ggrepel::geom_text_repel(ggplot2::aes(label = paste0("TH", Treehole)), size = 2.6, seed = 1) +
  ggplot2::geom_text(data = ss.r, ggplot2::aes(label = sprintf("r = %.2f", r)), x = Inf, y = -Inf,
                     hjust = 1.1, vjust = -0.6, size = 3.4) +
  ggplot2::facet_grid(Original ~ New, scales = "free", switch = "y") +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "Tree-hole size or shape", y = NULL) +
  ggplot2::theme_bw() + ggplot2::theme(strip.placement = "outside")
save_pdf(p.ss2, file.path(dirOut, "size_shape_vs_original.pdf"), 16, 14)

vars.full = setdiff(vars.all, setdiff(vars.size, vars.size.basis)) # the other size variables are functions of Size and Shape
vif = prune_vif(env, vars.full, threshold = 10)
cat("Variables dropped for VIF > 10:", setdiff(vars.full, vif$vars), "\n")
write.csv(data.frame(Variable = names(vif$vif), VIF = round(vif$vif, 2)),
          file.path(dirOut, "predictor_vif.csv"), row.names = F)
vars.cand = vif$vars

# Standardise the predictors so that biplot arrows are comparable
envs = env
envs[, vars.all] = scale(env[, vars.all])

# 2. FORWARD SELECTION ---------------------------------------------------------
fwd_select = function(resp) {
  mod0 = rda(resp ~ 1, data = envs)
  modF = rda(as.formula(paste("resp ~", paste(vars.cand, collapse = " + "))), data = envs)
  set.seed(1)
  sel = ordiR2step(mod0, scope = formula(modF), R2scope = TRUE, direction = "forward",
                   permutations = how(nperm = NPERM), trace = FALSE)
  list(full = modF, sel = sel)
}
fs = list(CLR = fwd_select(Y), PCoA = fwd_select(Y3))

fwd.tab = do.call(rbind, lapply(names(fs), function(r)
  data.frame(Response = responses[r], Step = rownames(fs[[r]]$sel$anova), fs[[r]]$sel$anova, row.names = NULL)))
write.csv(fwd.tab, file.path(dirOut, "forward_selection.csv"), row.names = F)
vars.sel = lapply(fs, function(f) attr(terms(f$sel), "term.labels"))

marg.clr = anova(fs$CLR$sel, by = "margin", permutations = how(nperm = NPERM))
marg.pc3 = anova(fs$PCoA$sel, by = "margin", permutations = how(nperm = NPERM))
write.csv(rbind(data.frame(Response = responses["CLR"], Term = rownames(marg.clr), marg.clr),
                data.frame(Response = responses["PCoA"], Term = rownames(marg.pc3), marg.pc3)),
          file.path(dirOut, "forward_selection_margin_tests.csv"), row.names = F)

# How much of each PCoA axis does each candidate explain on its own? (envfit-like,
# with the adjusted R2)
single = do.call(rbind, lapply(vars.all, function(v) {
  data.frame(Variable = v,
             t(sapply(1:3, function(a) round(summary(lm(Y3[, a] ~ env[, v]))$adj.r.squared, 3))),
             PCoA1_3 = round(r2_gram(G3, env[, v, drop = F])["adjR2"], 3),
             CLR = round(r2_gram(G, env[, v, drop = F])["adjR2"], 3))
}))
names(single)[2:4] = c("PCoA1", "PCoA2", "PCoA3")
write.csv(single[order(-single$PCoA1_3), ], file.path(dirOut, "single_variable_adjR2.csv"), row.names = F)

# 3. TREE-HOLE MORPHOLOGY: ALL COMBINATIONS, TREE-HOLE PERMUTATIONS ----------------
# The R2 of a variable that is constant within tree holes cannot exceed the R2 of
# the tree-hole identity itself.
X.treehole = model.matrix(~ Treehole, env)[, -1]
r2.treehole = list(CLR = r2_gram(G, X.treehole), PCoA = r2_gram(G3, X.treehole))
cat("\nTree-hole identity (9 df): adj. R2 CLR =", round(r2.treehole$CLR["adjR2"], 3),
    "; PCoA1-3 =", round(r2.treehole$PCoA["adjR2"], 3), "\n")

subs.l = list(CLR = morphology_subsets(G, env, nperm = NPERM.TH),
              PCoA = morphology_subsets(G3, env, nperm = NPERM.TH))
subs = merge(subs.l$CLR, subs.l$PCoA, by = c("Variables", "k"), suffixes = c(".CLR", ".PCoA1_3"))
subs = subs[order(-subs$adjR2.CLR), ]
subs$Fraction.of.treehole.CLR = round(subs$adjR2.CLR / r2.treehole$CLR["adjR2"], 3)
subs$Fraction.of.treehole.PCoA1_3 = round(subs$adjR2.PCoA1_3 / r2.treehole$PCoA["adjR2"], 3)
write.csv(subs, file.path(dirOut, "morphology_subsets_treehole_permutation.csv"), row.names = F)
cat("\nBest morphology combinations (CLR response):\n")
print(head(subs[, c("Variables", "adjR2.CLR", "P.treehole.perm.CLR", "P.selection.corrected.CLR",
                    "adjR2.PCoA1_3", "P.treehole.perm.PCoA1_3", "P.selection.corrected.PCoA1_3")], 10))

top.morph = lapply(subs.l, function(x) strsplit(x$Variables[1], " + ", fixed = T)[[1]]) # highest adj. R2
best.morph = vars.size.basis # morphology used in the varpart and the final model
subs_row = function(r, vars) subs.l[[r]][subs.l[[r]]$Variables == paste(vars, collapse = " + "), ]

# 4. TIME-VARYING VARIABLES WITHIN TREE HOLES -----------------------------------
vars.temp.cand = intersect(vars.cand, vars.temporal)
f.temp = as.formula(paste("Y ~", paste(vars.temp.cand, collapse = " + "), "+ Condition(Treehole)"))
prda.temp = rda(f.temp, data = envs)
set.seed(1)
marg.temp = anova(prda.temp, by = "margin", permutations = how(blocks = envs$Treehole, nperm = NPERM))
cat("\nPartial RDA of time-varying variables | tree hole (permutations within tree holes):\n")
print(marg.temp)
temp.tab = data.frame(Term = rownames(marg.temp), marg.temp, check.names = F)
temp.tab$Fraction.of.total.variance = temp.tab$Variance / prda.temp$tot.chi
write.csv(temp.tab, file.path(dirOut, "temporal_partialRDA_given_treehole.csv"), row.names = F)
vars.temp.sig = rownames(marg.temp)[!is.na(marg.temp[, "Pr(>F)"]) & marg.temp[, "Pr(>F)"] < 0.05]

# 5. VARIATION PARTITIONING --------------------------------------------------------
X.temp = if (length(vars.temp.sig)) envs[, vars.temp.sig, drop = F] else envs[, vars.temp.cand, drop = F]
vp = varpart(Y, envs[, best.morph, drop = F], X.temp)
vp.th = varpart(Y, X.treehole, X.temp)
pdf(file.path(dirOut, "varpart_CLR.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(vp, Xnames = c("Morphology", "Time-varying"), bg = c("#2271B3", "#D6AE01"), digits = 2)
title("Tree-hole morphology")
plot(vp.th, Xnames = c("Tree hole", "Time-varying"), bg = c("#497E76", "#D6AE01"), digits = 2)
title("Tree-hole identity")
dev.off()
vp.tab = rbind(data.frame(Partition = "Morphology vs time-varying", Fraction = rownames(vp$part$indfract), vp$part$indfract),
               data.frame(Partition = "Tree hole vs time-varying", Fraction = rownames(vp.th$part$indfract), vp.th$part$indfract))
write.csv(vp.tab, file.path(dirOut, "varpart_CLR.csv"), row.names = F)

# 6. FINAL MODEL AND FIGURES -----------------------------------------------------------
vars.final = c(best.morph, vars.temp.sig)
mod.final = rda(as.formula(paste("Y ~", paste(vars.final, collapse = " + "))), data = envs)
set.seed(1)
axes.final = anova(mod.final, by = "axis", permutations = how(nperm = NPERM), cutoff = 0.1)
write.csv(data.frame(Axis = rownames(axes.final), axes.final), file.path(dirOut, "final_model_axis_tests.csv"), row.names = F)
# Biplot scores: correlation of each predictor with the constrained axes (the arrows)
bp.final = scores(mod.final, display = "bp", choices = 1:3)
write.csv(data.frame(Variable = rownames(bp.final), round(bp.final, 3)),
          file.path(dirOut, "final_model_biplot_scores.csv"), row.names = F)

save_pdf(rda_biplot(fs$CLR$sel, envs, c(1, 2), "Forward-selected RDA (CLR)"), file.path(dirOut, "RDA_forward_selected_axes1_2.pdf"))
save_pdf(rda_biplot(fs$CLR$sel, envs, c(1, 3), "Forward-selected RDA (CLR)"), file.path(dirOut, "RDA_forward_selected_axes1_3.pdf"))
save_pdf(rda_biplot(mod.final, envs, c(1, 2), "Morphology + time-varying RDA (CLR)"), file.path(dirOut, "RDA_final_axes1_2.pdf"))
if (length(vars.final) >= 3)
  save_pdf(rda_biplot(mod.final, envs, c(1, 3), "Morphology + time-varying RDA (CLR)"), file.path(dirOut, "RDA_final_axes1_3.pdf"))

# RESULTS SUMMARY ----------------------------------------------------------------------
fsub = "morphology_subsets_treehole_permutation.csv"
rows = list()
for (r in names(responses)) {
  rows[[length(rows) + 1]] = summary_row("1", "Tree-hole identity (ceiling for morphology)", responses[r], "Treehole (9 df)",
                                         "none", adjR2 = r2.treehole[[r]]["adjR2"], File = "summary.txt")
  rows[[length(rows) + 1]] = summary_row("2", "Forward selection (exploratory)", responses[r], paste(vars.sel[[r]], collapse = " + "),
                                         "free (samples independent)", adjR2 = RsquareAdj(fs[[r]]$sel)$adj.r.squared,
                                         P = tail(fs[[r]]$sel$anova[, "Pr(>F)"], 2)[1], File = "forward_selection.csv")
  for (v in unique(list("Size", "Shape", c("Size", "Shape"), top.morph[[r]]))) {
    sr = subs_row(r, v)
    rows[[length(rows) + 1]] = summary_row("3", if (identical(v, top.morph[[r]])) "Morphology: highest adj. R2 of all combinations" else "Morphology",
                                           responses[r], sr$Variables, "tree holes (10)", adjR2 = sr$adjR2,
                                           P = sr$P.treehole.perm, P.corrected = sr$P.selection.corrected, File = fsub)
  }
}
for (v in rownames(marg.temp)[rownames(marg.temp) != "Residual"])
  rows[[length(rows) + 1]] = summary_row("4", "Time-varying | tree hole (partial RDA, marginal)", responses["CLR"], v,
                                         "within tree holes", Statistic = sprintf("F = %.2f; fraction = %.3f", marg.temp[v, "F"],
                                                                                 temp.tab$Fraction.of.total.variance[temp.tab$Term == v]),
                                         P = marg.temp[v, "Pr(>F)"], File = "temporal_partialRDA_given_treehole.csv")
fr = vp$part$indfract$Adj.R.squared
fr.th = vp.th$part$indfract$Adj.R.squared
rows[[length(rows) + 1]] = summary_row("5", "Varpart: morphology alone", responses["CLR"], paste(best.morph, collapse = " + "), "-", adjR2 = fr[1], File = "varpart_CLR.csv")
rows[[length(rows) + 1]] = summary_row("5", "Varpart: time-varying alone", responses["CLR"], paste(colnames(X.temp), collapse = " + "), "-", adjR2 = fr[2], File = "varpart_CLR.csv")
rows[[length(rows) + 1]] = summary_row("5", "Varpart: shared morphology/time-varying", responses["CLR"], "", "-", adjR2 = fr[3], File = "varpart_CLR.csv")
rows[[length(rows) + 1]] = summary_row("5", "Varpart: tree hole alone (vs time-varying)", responses["CLR"], "Treehole", "-", adjR2 = fr.th[1], File = "varpart_CLR.csv")
rows[[length(rows) + 1]] = summary_row("6", "Final model", responses["CLR"], paste(vars.final, collapse = " + "),
                                       "free (axis tests)", adjR2 = RsquareAdj(mod.final)$adj.r.squared,
                                       Statistic = paste("significant axes:", sum(axes.final[, "Pr(>F)"] < 0.05, na.rm = T)),
                                       File = "final_model_axis_tests.csv")
res.summary = do.call(rbind, rows)
write.csv(res.summary, file.path(dirOut, "results_summary.csv"), row.names = F, na = "")
write_combined_summary()

# Figure: variance explained by each analysis
expl = do.call(rbind, lapply(names(responses), function(r) data.frame(
  Response = responses[r],
  Analysis = c("Tree-hole identity", "Size", "Shape", "Size + Shape",
               paste(top.morph[[r]], collapse = " + "), if (r == "CLR") "Time-varying (unique)"),
  adjR2 = c(r2.treehole[[r]]["adjR2"], subs_row(r, "Size")$adjR2, subs_row(r, "Shape")$adjR2,
            subs_row(r, c("Size", "Shape"))$adjR2, subs_row(r, top.morph[[r]])$adjR2, if (r == "CLR") fr[2]))))
expl$adjR2 = pmax(expl$adjR2, 0)
save_pdf(plot_explained(expl, "Composition: variance explained"), file.path(dirOut, "variance_explained.pdf"), 20, 8)
write.csv(expl, file.path(dirOut, "variance_explained.csv"), row.names = F)
save_pdf(plot_explained_grouped(expl), file.path(dirOut, "variance_explained_combined.pdf"), 12, 10)

# SUMMARY ---------------------------------------------------------------------------
sink(file.path(dirOut, "summary.txt"))
cat("RDA of community composition (CLR, Aitchison geometry) on environmental variables\n")
cat("Samples:", nrow(env), " ASVs:", ncol(Y), "\n")
cat("Aitchison PCoA % variance axes 1-3:", round(pcoa.var[1:3], 1), "\n\n")
cat("Size/Shape correlations:\n"); print(cor.size, row.names = F)
cat("\nDropped for VIF > 10:", setdiff(vars.full, vars.cand), "\n")
cat("\nForward selection, CLR response:", vars.sel$CLR, "\n")
cat("  adj. R2 =", round(RsquareAdj(fs$CLR$sel)$adj.r.squared, 3), "\n")
cat("Forward selection, PCoA1-3 response:", vars.sel$PCoA, "\n")
cat("  adj. R2 =", round(RsquareAdj(fs$PCoA$sel)$adj.r.squared, 3), "\n")
cat("\nTree-hole identity adj. R2: CLR", round(r2.treehole$CLR["adjR2"], 3), "; PCoA1-3", round(r2.treehole$PCoA["adjR2"], 3), "\n")
cat("\nMorphology combinations (tree-hole permutations, n = 10 tree holes):\n")
print(head(subs[, c("Variables", "adjR2.CLR", "P.treehole.perm.CLR", "P.selection.corrected.CLR",
                    "adjR2.PCoA1_3", "P.treehole.perm.PCoA1_3", "P.selection.corrected.PCoA1_3")], 12), row.names = F)
cat("\nHighest adj. R2: CLR", top.morph$CLR, "; PCoA1-3", top.morph$PCoA, "\n")
cat("Morphology in varpart and final model (a priori):", best.morph, "\n")
cat("\nTime-varying variables | tree hole (within-tree-hole permutations):\n")
print(marg.temp)
cat("\nVariation partitioning (adj. R2):\n")
print(vp.tab[, c("Partition", "Df", "Adj.R.squared")])
cat("\nFinal model:", vars.final, " adj. R2 =", round(RsquareAdj(mod.final)$adj.r.squared, 3), "\n")
print(axes.final)
cat("\nBiplot scores of the final model:\n")
print(round(bp.final, 3))
sink()
cat("Results written to", dirOut, "\n")
