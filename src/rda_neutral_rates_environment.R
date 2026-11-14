###################################
#  rda_neutral_rates_environment.R
##################################
# Same question as rda_composition_environment.R, but with the parameters of the
# neutral model fitted by NMGS as responses instead of community composition:
# which combination of environmental variables explains the immigration and
# speciation rates?
#
# The two rates live at different levels:
#  * Immigration rate (Immigrate_med): one value per sample.
#  * Speciation rate (Speciation_med): one value per metacommunity (Dirichlet
#    component), shared by all its samples, so there are only 10 independent
#    values. Tree-hole morphology is constant per tree hole (10 values).
#
# Analyses:
#  1. Both rates jointly (RDA on the standardised log rates), as for composition:
#     forward selection, morphology combinations with tree-hole permutations,
#     time-varying variables within tree holes, variation partitioning.
#  2. Immigration rate alone: the same steps, a variation partitioning that adds
#     the metacommunity, and a linear mixed model with tree hole as random effect
#     (correct degrees of freedom for tree-hole level predictors).
#  3. Speciation rate alone, at the metacommunity level (n = 10): Spearman
#     correlations with the metacommunity means of the predictors.
#
# Input:  1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv
# Output: 4.1_Niche_Envir/rda/neutral_rates/  (results_summary.csv, tables .csv,
#         figures .pdf, summary.txt) and 4.1_Niche_Envir/rda/results_summary.csv
###################################

rm(list = ls())
library(here)
library(vegan)
library(permute)
library(nlme)
library(ggplot2)
source(here("src", "rda_functions.R"))

dirOut = here("4.1_Niche_Envir", "rda", "neutral_rates")
dir.create(dirOut, recursive = TRUE, showWarnings = FALSE)
NPERM = 999
NPERM.TH = 9999

# READ DATA -----------------------------------------------------------------
env = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv"),
                 header = T, sep = "\t")
env = add_predictors(env)
env$log.Immigration = log(env$Immigrate_med)
env$log.Speciation = log(env$Speciation_med)

vars.all = c(vars.morphology, vars.temporal)
keep = complete.cases(env[, c(vars.all, "log.Immigration", "log.Speciation")])
cat("Samples:", nrow(env), "; complete cases:", sum(keep), "\n")
env = env[keep, ]

vars.full = setdiff(vars.all, setdiff(vars.size, vars.size.basis)) # the other size variables are functions of Size and Shape
vif = prune_vif(env, vars.full, threshold = 10)
vars.cand = vif$vars
cat("Variables dropped for VIF > 10:", setdiff(vars.full, vars.cand), "\n")

envs = env
envs[, vars.all] = scale(env[, vars.all])
X.treehole = model.matrix(~ Treehole, env)[, -1]
X.metacomm = model.matrix(~ Metacommunity, env)[, -1]

# PIPELINE (steps 2-5 of rda_composition_environment.R) --------------------------------
run_pipeline = function(Y, label) {
  Y = as.matrix(Y)
  G = gram(Y)
  out = list(label = label)

  # Forward selection (exploratory: samples treated as independent)
  mod0 = rda(Y ~ 1, data = envs)
  modF = rda(as.formula(paste("Y ~", paste(vars.cand, collapse = " + "))), data = envs)
  set.seed(1)
  sel = ordiR2step(mod0, scope = formula(modF), R2scope = TRUE, direction = "forward",
                   permutations = how(nperm = NPERM), trace = FALSE)
  out$sel = sel
  out$vars.sel = attr(terms(sel), "term.labels")
  out$fwd = data.frame(Response = label, Step = rownames(sel$anova), sel$anova, row.names = NULL)
  out$adjR2.full = RsquareAdj(modF)$adj.r.squared
  out$adjR2.sel = RsquareAdj(sel)$adj.r.squared

  # Single variables
  out$single = do.call(rbind, lapply(vars.all, function(v)
    data.frame(Response = label, Variable = v, adjR2 = round(r2_gram(G, env[, v, drop = F])["adjR2"], 3))))

  # Ceilings: everything constant within tree holes / within metacommunities
  out$adjR2.treehole = r2_gram(G, X.treehole)["adjR2"]
  out$adjR2.metacomm = r2_gram(G, X.metacomm)["adjR2"]

  # Morphology combinations, tree-hole permutations
  subs = morphology_subsets(G, env, nperm = NPERM.TH)
  subs$Fraction.of.treehole = round(subs$adjR2 / out$adjR2.treehole, 3)
  out$subs = data.frame(Response = label, subs)
  out$top.morph = strsplit(subs$Variables[1], " + ", fixed = T)[[1]] # highest adj. R2
  out$best.morph = vars.size.basis # a priori morphology (Size + Shape) for varpart and mixed model

  # Time-varying variables given tree hole
  vars.temp.cand = intersect(vars.cand, vars.temporal)
  prda = rda(as.formula(paste("Y ~", paste(vars.temp.cand, collapse = " + "), "+ Condition(Treehole)")), data = envs)
  set.seed(1)
  marg = anova(prda, by = "margin", permutations = how(blocks = envs$Treehole, nperm = NPERM))
  out$temp = data.frame(Response = label, Term = rownames(marg), marg, row.names = NULL)
  out$temp.sig = rownames(marg)[!is.na(marg[, "Pr(>F)"]) & marg[, "Pr(>F)"] < 0.05]

  # Variation partitioning
  X.temp = envs[, if (length(out$temp.sig)) out$temp.sig else vars.temp.cand, drop = F]
  X.morph = envs[, out$best.morph, drop = F]
  out$vp = varpart(Y, X.morph, X.temp)
  out$vp.th = varpart(Y, X.treehole, X.temp)
  out
}

print_pipeline = function(o) {
  cat("\n=====", o$label, "=====\n")
  cat("Forward selection:", o$vars.sel, " adj. R2 =", round(o$adjR2.sel, 3),
      "(full model:", round(o$adjR2.full, 3), ")\n")
  cat("Ceilings (adj. R2): tree hole", round(o$adjR2.treehole, 3), "; metacommunity", round(o$adjR2.metacomm, 3), "\n")
  cat("Morphology combinations (tree-hole permutations):\n")
  print(head(o$subs[, c("Variables", "adjR2", "P.treehole.perm", "P.selection.corrected", "Fraction.of.treehole")], 10), row.names = F)
  cat("Highest adj. R2:", o$top.morph, "; morphology used downstream (a priori):", o$best.morph, "\n")
  cat("Time-varying variables | tree hole (within-tree-hole permutations):\n")
  print(o$temp[, -1], row.names = F)
  cat("Variation partitioning, morphology vs time-varying (adj. R2):\n")
  print(round(o$vp$part$indfract[, c("Df", "Adj.R.squared")], 3))
  cat("Variation partitioning, tree hole vs time-varying (adj. R2):\n")
  print(round(o$vp.th$part$indfract[, c("Df", "Adj.R.squared")], 3))
}

save_pipeline = function(o, prefix) {
  write.csv(o$fwd, file.path(dirOut, paste0(prefix, "_forward_selection.csv")), row.names = F)
  write.csv(o$single[order(-o$single$adjR2), ], file.path(dirOut, paste0(prefix, "_single_variable_adjR2.csv")), row.names = F)
  write.csv(o$subs, file.path(dirOut, paste0(prefix, "_morphology_subsets_treehole_permutation.csv")), row.names = F)
  write.csv(o$temp, file.path(dirOut, paste0(prefix, "_temporal_partialRDA_given_treehole.csv")), row.names = F)
  vp = rbind(data.frame(Partition = "Morphology vs time-varying", Fraction = rownames(o$vp$part$indfract), o$vp$part$indfract),
             data.frame(Partition = "Tree hole vs time-varying", Fraction = rownames(o$vp.th$part$indfract), o$vp.th$part$indfract))
  write.csv(vp, file.path(dirOut, paste0(prefix, "_varpart.csv")), row.names = F)
}

# 1. BOTH RATES JOINTLY ----------------------------------------------------------------
Y.both = scale(env[, c("log.Immigration", "log.Speciation")])
res.both = run_pipeline(Y.both, "Immigration + speciation (standardised log rates)")
save_pipeline(res.both, "rates_joint")
if (length(res.both$vars.sel) >= 1) {
  p = rda_biplot(res.both$sel, envs, c(1, if (length(res.both$vars.sel) > 1) 2 else 1),
                 "Forward-selected RDA of the neutral rates")
  # add the response (species) scores: the two rates
  sp = data.frame(scores(res.both$sel, display = "species", choices = 1:2, scaling = 2))
  sp$Rate = c("Immigration", "Speciation")
  if (ncol(sp) == 3) save_pdf(p + geom_text(data = sp, aes(x = .data[[names(sp)[1]]], y = .data[[names(sp)[2]]], label = Rate),
                                            color = "red", fontface = "bold"),
                              file.path(dirOut, "RDA_rates_selected_axes1_2.pdf"))
}

# 2. IMMIGRATION RATE ---------------------------------------------------------------
res.imm = run_pipeline(env[, "log.Immigration", drop = F], "Immigration rate (log)")
save_pipeline(res.imm, "immigration")

# Does the metacommunity explain immigration beyond the environment?
X.final.imm = envs[, unique(c(res.imm$best.morph, res.imm$temp.sig)), drop = F]
vp.imm.meta = if (ncol(X.final.imm)) varpart(env$log.Immigration, X.final.imm, X.metacomm) else NULL
if (!is.null(vp.imm.meta)) {
  pdf(file.path(dirOut, "immigration_varpart.pdf"), width = 15, height = 5)
  par(mfrow = c(1, 3))
  plot(res.imm$vp, Xnames = c("Morphology", "Time-varying"), bg = c("#2271B3", "#D6AE01"), digits = 2)
  plot(res.imm$vp.th, Xnames = c("Tree hole", "Time-varying"), bg = c("#497E76", "#D6AE01"), digits = 2)
  plot(vp.imm.meta, Xnames = c("Environment", "Metacommunity"), bg = c("#2271B3", "#C35831"), digits = 2)
  dev.off()
}

# Mixed model: log immigration ~ selected variables + (1 | tree hole). Tree-hole level
# predictors get the between-tree-hole degrees of freedom (10 tree holes).
vars.lme = unique(c(res.imm$best.morph, res.imm$temp.sig))
if (!length(vars.lme)) vars.lme = res.imm$vars.sel
f.lme = as.formula(paste("log.Immigration ~", paste(vars.lme, collapse = " + ")))
lme.imm = lme(f.lme, random = ~ 1 | Treehole, data = envs, method = "REML")
lme.tab = data.frame(Term = rownames(summary(lme.imm)$tTable), summary(lme.imm)$tTable, row.names = NULL)
ci = intervals(lme.imm, which = "fixed")$fixed
lme.tab$CI.low = ci[, "lower"]; lme.tab$CI.high = ci[, "upper"]
write.csv(lme.tab, file.path(dirOut, "immigration_lme_coefficients.csv"), row.names = F)
vc = VarCorr(lme.imm)
icc = as.numeric(vc[1, 1]) / (as.numeric(vc[1, 1]) + as.numeric(vc[2, 1]))

coef.plot = ggplot(lme.tab[lme.tab$Term != "(Intercept)", ],
                   aes(x = Value, y = reorder(Term, Value))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = CI.low, xmax = CI.high), height = 0.2) +
  geom_point(size = 3) +
  labs(x = "Effect on log immigration rate (per SD, 95% CI)", y = NULL,
       title = "Mixed model, random intercept per tree hole") +
  theme_bw()
save_pdf(coef.plot, file.path(dirOut, "immigration_lme_coefficients.pdf"), 16, 10)

# 3. SPECIATION RATE, METACOMMUNITY LEVEL ---------------------------------------------
meta.lev = aggregate(env[, c("log.Speciation", vars.all)], by = list(Metacommunity = env$Metacommunity), FUN = mean)
meta.lev$N = as.vector(table(env$Metacommunity)[as.character(meta.lev$Metacommunity)])
spec.cor = do.call(rbind, lapply(vars.all, function(v) {
  ct = suppressWarnings(cor.test(meta.lev$log.Speciation, meta.lev[[v]], method = "spearman", exact = TRUE))
  data.frame(Variable = v, rho = round(unname(ct$estimate), 3), P = signif(ct$p.value, 3))
}))
spec.cor$P.BH = signif(p.adjust(spec.cor$P, "BH"), 3)
spec.cor = spec.cor[order(spec.cor$P), ]
write.csv(spec.cor, file.path(dirOut, "speciation_metacommunity_spearman.csv"), row.names = F)
write.csv(meta.lev, file.path(dirOut, "speciation_metacommunity_means.csv"), row.names = F)

top = spec.cor$Variable[1]
spec.plot = ggplot(meta.lev, aes(x = .data[[top]], y = exp(log.Speciation))) +
  geom_point(aes(color = Metacommunity, size = N)) +
  scale_color_manual(values = custom_palette) +
  scale_y_log10() +
  labs(x = paste("Metacommunity mean of", top), y = "Speciation rate (median, log scale)",
       title = sprintf("Spearman rho = %.2f, P = %.3g (n = 10)", spec.cor$rho[1], spec.cor$P[1])) +
  theme_bw()
save_pdf(spec.plot, file.path(dirOut, "speciation_metacommunity_top_variable.pdf"), 16, 12)

# FIGURE: IMMIGRATION RATE AGAINST TREE-HOLE SIZE AND SHAPE -----------------------------
# Samples (small points) and tree-hole means (large points); line = mixed-model slope
# when the variable is in the model, otherwise a least-squares fit of the tree-hole means.
th.means = aggregate(env[, c("log.Immigration", "Size", "Shape")], by = list(Treehole = env$Treehole), FUN = mean)
long = rbind(data.frame(Variable = "Size (~ log volume)", x = env$Size, y = env$log.Immigration, Treehole = env$Treehole, Level = "Sample"),
             data.frame(Variable = "Shape (~ log depth:surface)", x = env$Shape, y = env$log.Immigration, Treehole = env$Treehole, Level = "Sample"))
long.th = rbind(data.frame(Variable = "Size (~ log volume)", x = th.means$Size, y = th.means$log.Immigration, Treehole = th.means$Treehole),
                data.frame(Variable = "Shape (~ log depth:surface)", x = th.means$Shape, y = th.means$log.Immigration, Treehole = th.means$Treehole))
long$Variable = factor(long$Variable, levels = unique(long$Variable))
long.th$Variable = factor(long.th$Variable, levels = levels(long$Variable))
imm.plot = ggplot(long, aes(x = x, y = y)) +
  geom_point(aes(shape = Treehole), color = "grey60", size = 1.5) +
  geom_smooth(data = long.th, method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey85", linewidth = 0.7) +
  geom_point(data = long.th, aes(shape = Treehole), size = 3.5, stroke = 1.1) +
  facet_wrap(~ Variable, scales = "free_x") +
  scale_shape_manual(values = c(0:9)) +
  labs(x = "Standardised value (tree-hole level)", y = "log immigration rate (NMGS median)",
       title = "Immigration rate vs tree-hole size and shape (large points: tree-hole means, n = 10)") +
  theme_bw() + theme(plot.title = element_text(size = 10))
save_pdf(imm.plot, file.path(dirOut, "immigration_vs_size_shape.pdf"), 22, 11)

# FIGURE: IMMIGRATION RATE AGAINST THE MORPHOLOGICAL MODELS -----------------------------
# Observed log immigration against the values fitted by each model (least squares,
# samples as in Table "immigration"): the a priori model (Size + Shape) and the
# combination with the highest adjusted R2. The fitted values take one value per
# tree hole. Line: regression of observed on fitted values (slope 1 by construction).
# P: permutation of the morphology among tree holes; P.sel: corrected for the
# selection of the best combination.
fit_models = list(res.imm$best.morph, res.imm$top.morph)
fit_tab = do.call(rbind, lapply(fit_models, function(v) {
  sr = res.imm$subs[res.imm$subs$Variables == paste(v, collapse = " + "), ]
  fit = lm(as.formula(paste("log.Immigration ~", paste(v, collapse = " + "))), data = env)
  lab = paste(v, collapse = " + ")
  lab = if (identical(v, res.imm$best.morph)) paste0(lab, " (a priori)") else paste0(lab, " (highest adj. R2)")
  data.frame(Model = lab, Fitted = fitted(fit), Observed = env$log.Immigration, Treehole = env$Treehole,
             Label = if (identical(v, res.imm$best.morph)) sprintf("adj. R2 = %.2f\nP = %.3f", sr$adjR2, sr$P.treehole.perm)
                     else sprintf("adj. R2 = %.2f\nP = %.3f; P.sel = %.2f", sr$adjR2, sr$P.treehole.perm, sr$P.selection.corrected))
}))
fit_tab$Model = factor(fit_tab$Model, levels = unique(fit_tab$Model))
fit_th = aggregate(cbind(Fitted, Observed) ~ Model + Treehole, data = fit_tab, FUN = mean)
fit_lab = fit_tab[!duplicated(fit_tab$Model), c("Model", "Label")]
imm.fit.plot = ggplot(fit_tab, aes(x = Fitted, y = Observed)) +
  geom_point(aes(shape = Treehole), color = "grey60", size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey85", linewidth = 0.7) +
  geom_point(data = fit_th, aes(shape = Treehole), size = 3.5, stroke = 1.1) +
  geom_text(data = fit_lab, aes(label = Label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.2) +
  facet_wrap(~ Model) +
  scale_shape_manual(values = c(0:9)) +
  labs(x = "Fitted log immigration rate (tree-hole morphology)", y = "log immigration rate (NMGS median)") +
  theme_bw()
save_pdf(imm.fit.plot, file.path(dirOut, "immigration_vs_morphology_fit.pdf"), 22, 11)

# Standalone panel with the a priori model only (manuscript Fig. 3B)
m1 = levels(fit_tab$Model)[1]
imm.fit1.plot = ggplot(fit_tab[fit_tab$Model == m1, ], aes(x = Fitted, y = Observed)) +
  geom_point(aes(shape = Treehole), color = "grey60", size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey85", linewidth = 0.7) +
  geom_point(data = fit_th[fit_th$Model == m1, ], aes(shape = Treehole), size = 3.5, stroke = 1.1) +
  annotate("text", label = fit_lab$Label[fit_lab$Model == m1], x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3, size = 3.4) +
  scale_shape_manual(values = c(0:9)) +
  labs(x = "Log immigration rate predicted by size + shape", y = "Log immigration rate (NMGS median)") +
  theme_bw()
save_pdf(imm.fit1.plot, file.path(dirOut, "immigration_vs_size_shape_fit.pdf"), 14, 11)

# DROUGHT PREVALENCE AND TREE-HOLE MORPHOLOGY -----------------------------------------
# Drought days is tree-hole specific: days since the tree hole was last found holding
# water (0 when it held water at sampling). Drought prevalence of a tree hole =
# maximum number of drought days over the survey (as in
# immigration_speciation_rates_analysis.R); the fraction of samplings in which the
# tree hole was dry is given as an alternative. All samples are used (not only those
# with complete predictors). Which morphological variable best represents water
# retention? Pearson correlations across the 10 tree holes.
env.all = add_predictors(read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv"),
                                    header = T, sep = "\t"))
th.dr = do.call(rbind, lapply(split(env.all, env.all$Treehole), function(x) data.frame(
  Treehole = x$Treehole[1], Drought.max = max(x$Drought_days, na.rm = T),
  Fraction.dry = mean(x$Wet == "D", na.rm = T), x[1, vars.morphology])))
vars.dr = c("log.Depth", "log.SA", "log.Volume", "log.Depth.surface", "Size", "Shape", "Northness", "Eastness")
dr.cor = do.call(rbind, lapply(c("Drought.max", "Fraction.dry"), function(r) do.call(rbind, lapply(vars.dr, function(v) {
  ct = cor.test(th.dr[[r]], th.dr[[v]])
  ct9 = cor.test(th.dr[[r]][th.dr$Treehole != "1"], th.dr[[v]][th.dr$Treehole != "1"])
  data.frame(Response = r, Variable = v, r = round(unname(ct$estimate), 3), P = signif(ct$p.value, 3),
             r.without.TH1 = round(unname(ct9$estimate), 3), P.without.TH1 = signif(ct9$p.value, 3))
}))))
write.csv(dr.cor, file.path(dirOut, "drought_vs_morphology.csv"), row.names = F)
write.csv(th.dr, file.path(dirOut, "drought_treehole_table.csv"), row.names = F)

dr.long = rbind(data.frame(Variable = "Depth (cm, log scale)", x = exp(th.dr$log.Depth), y = th.dr$Drought.max, Treehole = th.dr$Treehole),
                data.frame(Variable = "Depth:surface ratio (log scale)", x = exp(th.dr$log.Depth.surface), y = th.dr$Drought.max, Treehole = th.dr$Treehole))
dr.long$Variable = factor(dr.long$Variable, levels = unique(dr.long$Variable))
dr.lab = data.frame(Variable = levels(dr.long$Variable), Label = sapply(c("log.Depth", "log.Depth.surface"), function(v) {
  d = dr.cor[dr.cor$Response == "Drought.max" & dr.cor$Variable == v, ]
  sprintf("r = %.2f, P = %.3f", d$r, d$P)
}))
dr.lab$Variable = factor(dr.lab$Variable, levels = levels(dr.long$Variable))
dr.plot = ggplot(dr.long, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey85", linewidth = 0.7) +
  geom_point(aes(shape = Treehole), size = 3.5, stroke = 1.1) +
  geom_text(data = dr.lab, aes(label = Label), x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 3.4) +
  facet_wrap(~ Variable, scales = "free_x") +
  scale_x_log10() +
  scale_shape_manual(values = c(0:9)) +
  labs(x = NULL, y = "Drought prevalence (maximum drought days)") +
  theme_bw()
save_pdf(dr.plot, file.path(dirOut, "drought_vs_depth_and_ratio.pdf"), 22, 11)

# Standalone panel against depth only (manuscript Fig. 3D)
dr1 = dr.long[dr.long$Variable == levels(dr.long$Variable)[1], ]
dr1.plot = ggplot(dr1, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "grey85", linewidth = 0.7) +
  geom_point(aes(shape = Treehole), size = 3.5, stroke = 1.1) +
  annotate("text", label = dr.lab$Label[1], x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 3.4) +
  scale_x_log10() +
  scale_shape_manual(values = c(0:9)) +
  labs(x = "Tree-hole depth (cm, log scale)", y = "Drought prevalence (maximum drought days)") +
  theme_bw()
save_pdf(dr1.plot, file.path(dirOut, "drought_vs_depth.pdf"), 14, 11)

# RESULTS SUMMARY ----------------------------------------------------------------------
pipeline_rows = function(o, prefix) {
  rows = list(
    summary_row("1", "Tree-hole identity (ceiling for morphology)", o$label, "Treehole (9 df)", "none", adjR2 = o$adjR2.treehole, File = "summary.txt"),
    summary_row("1", "Metacommunity identity", o$label, "Metacommunity (9 df)", "none", adjR2 = o$adjR2.metacomm, File = "summary.txt"),
    summary_row("2", "Forward selection (exploratory)", o$label, paste(o$vars.sel, collapse = " + "), "free (samples independent)",
                adjR2 = o$adjR2.sel, File = paste0(prefix, "_forward_selection.csv")))
  for (v in unique(list("Size", "Shape", c("Size", "Shape"), o$top.morph))) {
    sr = o$subs[o$subs$Variables == paste(v, collapse = " + "), ]
    rows[[length(rows) + 1]] = summary_row("3", if (identical(v, o$top.morph)) "Morphology: highest adj. R2 of all combinations" else "Morphology",
                                           o$label, sr$Variables, "tree holes (10)", adjR2 = sr$adjR2, P = sr$P.treehole.perm,
                                           P.corrected = sr$P.selection.corrected, File = paste0(prefix, "_morphology_subsets_treehole_permutation.csv"))
  }
  for (i in which(o$temp$Term != "Residual"))
    rows[[length(rows) + 1]] = summary_row("4", "Time-varying | tree hole (partial RDA, marginal)", o$label, o$temp$Term[i], "within tree holes",
                                           Statistic = sprintf("F = %.2f", o$temp$F[i]), P = o$temp[i, "Pr..F."],
                                           File = paste0(prefix, "_temporal_partialRDA_given_treehole.csv"))
  fr = o$vp$part$indfract$Adj.R.squared
  rows[[length(rows) + 1]] = summary_row("5", "Varpart: morphology alone", o$label, paste(o$best.morph, collapse = " + "), "-", adjR2 = fr[1], File = paste0(prefix, "_varpart.csv"))
  rows[[length(rows) + 1]] = summary_row("5", "Varpart: time-varying alone", o$label, "", "-", adjR2 = fr[2], File = paste0(prefix, "_varpart.csv"))
  do.call(rbind, rows)
}
rows = list(pipeline_rows(res.both, "rates_joint"), pipeline_rows(res.imm, "immigration"))
if (!is.null(vp.imm.meta)) {
  fm = vp.imm.meta$part$indfract$Adj.R.squared
  rows[[length(rows) + 1]] = rbind(
    summary_row("5", "Varpart: environment alone (vs metacommunity)", res.imm$label, paste(colnames(X.final.imm), collapse = " + "), "-", adjR2 = fm[1], File = "immigration_varpart.pdf"),
    summary_row("5", "Varpart: shared environment/metacommunity", res.imm$label, "", "-", adjR2 = fm[3], File = "immigration_varpart.pdf"),
    summary_row("5", "Varpart: metacommunity alone", res.imm$label, "Metacommunity", "-", adjR2 = fm[2], File = "immigration_varpart.pdf"))
}
lt = lme.tab[lme.tab$Term != "(Intercept)", ]
rows[[length(rows) + 1]] = summary_row("6", "Mixed model, random intercept per tree hole", res.imm$label, lt$Term, "none (t test, tree-hole df)",
                                       Statistic = sprintf("beta = %.3f per SD [%.3f, %.3f]; t = %.2f, df = %d; ICC = %.2f",
                                                           lt$Value, lt$CI.low, lt$CI.high, lt$t.value, lt$DF, icc),
                                       P = lt$p.value, File = "immigration_lme_coefficients.csv")
rows[[length(rows) + 1]] = summary_row("7", "Spearman, metacommunity means (n = 10)", "Speciation rate (log)", head(spec.cor$Variable, 3),
                                       "exact", Statistic = sprintf("rho = %.2f", head(spec.cor$rho, 3)), P = head(spec.cor$P, 3),
                                       P.corrected = head(spec.cor$P.BH, 3), File = "speciation_metacommunity_spearman.csv")
res.summary = do.call(rbind, rows)
write.csv(res.summary, file.path(dirOut, "results_summary.csv"), row.names = F, na = "")
write_combined_summary()

expl = do.call(rbind, lapply(list(res.both, res.imm), function(o) data.frame(
  Response = o$label,
  Analysis = c("Tree-hole identity", "Metacommunity identity", "Size", "Shape", "Size + Shape",
               paste(res.imm$top.morph, collapse = " + "), "Time-varying (unique)"),
  adjR2 = c(o$adjR2.treehole, o$adjR2.metacomm,
            sapply(c("Size", "Shape", "Size + Shape", paste(res.imm$top.morph, collapse = " + ")),
                   function(v) o$subs$adjR2[o$subs$Variables == v]),
            o$vp$part$indfract$Adj.R.squared[2]))))
expl$adjR2 = pmax(expl$adjR2, 0)
save_pdf(plot_explained(expl, "Neutral rates: variance explained"), file.path(dirOut, "variance_explained.pdf"), 20, 8)
# Same data in one panel (manuscript Fig. 3C), as variance_explained_combined.pdf for composition.
# The joint response is an RDA on both standardised log rates: its R2 is the mean of the R2
# of the two rates, and the speciation rate takes only one value per metacommunity.
expl.g = expl
expl.g$Response = c("Immigration + speciation (joint RDA)", "Immigration rate (log)")[match(expl.g$Response, c(res.both$label, res.imm$label))]
save_pdf(plot_explained_grouped(expl.g, colors = c("#2271B3", "#D6AE01")) + guides(fill = guide_legend(nrow = 2)),
         file.path(dirOut, "variance_explained_combined.pdf"), 12, 10)
save_pdf(plot_explained(expl[expl$Response == res.imm$label, ], NULL), file.path(dirOut, "variance_explained_immigration.pdf"), 10, 8)

# SUMMARY ---------------------------------------------------------------------------
sink(file.path(dirOut, "summary.txt"))
cat("Neutral-model rates (NMGS, K = 10) as responses of environmental variables\n")
cat("Samples:", nrow(env), "; tree holes: 10; metacommunities: 10\n")
cat("Correlation of log immigration and log speciation across samples:",
    round(cor(env$log.Immigration, env$log.Speciation), 3), "\n")
cat("Dropped for VIF > 10:", setdiff(vars.full, vars.cand), "\n")
print_pipeline(res.both)
print_pipeline(res.imm)
if (!is.null(vp.imm.meta)) {
  cat("\nImmigration: environment vs metacommunity (adj. R2):\n")
  print(round(vp.imm.meta$part$indfract[, c("Df", "Adj.R.squared")], 3))
}
cat("\nMixed model for log immigration (random intercept per tree hole):\n")
print(summary(lme.imm)$tTable)
cat("Intraclass correlation (tree hole):", round(icc, 3), "\n")
cat("\nSpeciation vs metacommunity means (Spearman, n = 10):\n")
print(spec.cor, row.names = F)
cat("\nDrought prevalence vs tree-hole morphology (Pearson, n = 10 tree holes):\n")
print(dr.cor, row.names = F)
sink()
cat("Results written to", dirOut, "\n")
