###################################
#  stegen_null_models.R
##################################
# Null models of the Stegen et al. (2013, ISME J) framework for every pair of
# tree-hole samples:
#  * betaNTI: standardised effect size of the beta mean nearest taxon distance
#    (betaMNTD) against a null model that shuffles the ASV labels across the
#    tips of the tree (picante::taxaShuffle). Computed both abundance-weighted
#    and unweighted (presence/absence), as picante::comdistnt with
#    exclude.conspecifics = FALSE.
#  * RC_bray: Raup-Crick index on Bray-Curtis (Stegen et al. 2013, after Chase
#    et al. 2011): null communities keep the richness and the number of reads of
#    each sample; ASVs are drawn with probability proportional to their
#    occupancy, and the remaining reads proportional to their total abundance.
#    Scaled to [-1, 1].
#
# Speed: instead of looping over pairs, each randomisation computes all pairs at
# once. betaMNTD for all pairs is (P M + t(P M)) / 2, where P holds the relative
# abundances (samples x ASVs) and M[t, s] is the distance from ASV t to the
# nearest ASV present in sample s. For RC_bray, each randomisation draws one null
# community per sample and compares all pairs; the null distribution of each pair
# is the same as in the original pairwise code.
# Randomisations run in parallel (forked processes). Each randomisation has its
# own L'Ecuyer-CMRG stream derived from `seed`, so the results do not depend on
# the number of cores.
#
# Replaces the code of J. Ramoneda (NTI_Stegen_treeholes.R, 2020), whose
# raup_crick_abundance() overrode the requested number of randomisations with
# reps = 10.
#
# Input:  1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv (unrarefied
#         ASVs x samples; mitochondria/chloroplasts and CONTROL already removed)
#         3.2_Neutral_Stegen/nj_tree_stegen.tre (phylogeny_treeholes_repset.R)
# Output (in outdir, default 3.2_Neutral_Stegen):
#         betaNTI_treeholes{suffix}.csv             weighted betaNTI
#         betaNTI_unweighted_treeholes{suffix}.csv  unweighted betaNTI
#         betaMNTD_treeholes{suffix}.csv            observed weighted betaMNTD
#         RC_bray_treeholes{suffix}.csv             Raup-Crick (Bray-Curtis)
#         stegen_null_models{suffix}.log            parameters and session info
#
# Usage:  Rscript src/stegen_null_models.R [reps=999] [cores=16] [seed=20261114]
#                                          [outdir=3.2_Neutral_Stegen] [suffix=]
###################################

suppressMessages({
  library(ape)
  library(vegan)
  library(parallel)
  library(here)
})

# --- Parameters ----------------------------------------------------------------
opt = list(reps = 999, cores = max(1, detectCores() - 4), seed = 20261114,
           outdir = here("3.2_Neutral_Stegen"), suffix = "")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv = strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opt)) stop("Unknown argument: ", a)
  opt[[kv[1]]] = if (is.numeric(opt[[kv[1]]])) as.numeric(kv[2]) else kv[2]
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out = function(name) file.path(opt$outdir, paste0(name, opt$suffix, ".csv"))

# --- Data ----------------------------------------------------------------------
tab = read.delim(here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv"),
                 row.names = 1, check.names = FALSE)
comm = t(as.matrix(tab))                       # samples x ASVs (read counts)
comm = comm[, colSums(comm) > 0, drop = FALSE]
tree = read.tree(here("3.2_Neutral_Stegen", "nj_tree_stegen.tre"))
stopifnot(all(colnames(comm) %in% tree$tip.label))
D = cophenetic.phylo(tree)[colnames(comm), colnames(comm)]

n.samp = nrow(comm); n.asv = ncol(comm)
pres = comm > 0
P.w = comm / rowSums(comm)                     # weighted: relative abundances
P.u = pres / rowSums(pres)                     # unweighted: 1/richness
present = lapply(seq_len(n.samp), function(s) which(pres[s, ]))
cat("Samples:", n.samp, " ASVs:", n.asv, " reps:", opt$reps, " cores:", opt$cores,
    " seed:", opt$seed, "\n")

# --- betaMNTD for all pairs, given a (possibly shuffled) labelling --------------
# perm = identity for the observed values. With shuffled labels, the distance
# between ASVs i and j is D[perm[i], perm[j]].
min.dist = function(perm) {
  M = matrix(0, n.asv, n.samp)
  for (s in seq_len(n.samp)) {
    cols = perm[present[[s]]]
    M[, s] = do.call(pmin, lapply(cols, function(j) D[, j]))
  }
  M[perm, , drop = FALSE]                      # row t = distance from ASV t
}
bmntd = function(P, M) { A = P %*% M; (A + t(A)) / 2 }

# --- Raup-Crick null communities --------------------------------------------------
occur = colSums(pres)                          # occupancy of each ASV
abund = colSums(comm)                          # total reads of each ASV
rich = rowSums(pres); reads = rowSums(comm)
null.comm = function() {
  t(vapply(seq_len(n.samp), function(s) {
    sp = sample.int(n.asv, rich[s], replace = FALSE, prob = occur)
    extra = sp[sample.int(length(sp), reads[s] - rich[s], replace = TRUE, prob = abund[sp])]
    as.numeric(tabulate(c(sp, extra), nbins = n.asv))
  }, numeric(n.asv)))
}
bray.obs = as.vector(vegdist(comm, method = "bray"))

# --- One randomisation -------------------------------------------------------------
one.rep = function(r) {
  assign(".Random.seed", streams[[r]], envir = globalenv())
  M = min.dist(sample.int(n.asv))
  nb = as.vector(vegdist(null.comm(), method = "bray"))
  list(w = bmntd(P.w, M), u = bmntd(P.u, M),
       less = nb < bray.obs, tie = nb == bray.obs)
}

# Independent random streams, one per randomisation
RNGkind("L'Ecuyer-CMRG")
set.seed(opt$seed)
streams = vector("list", opt$reps)
streams[[1]] = .Random.seed
for (r in seq_len(opt$reps)[-1]) streams[[r]] = nextRNGStream(streams[[r - 1]])

# --- Observed values -----------------------------------------------------------------
t0 = Sys.time()
M.obs = min.dist(seq_len(n.asv))
obs.w = bmntd(P.w, M.obs); obs.u = bmntd(P.u, M.obs)

# --- Randomisations, in chunks (one chunk per core) -------------------------------------
chunks = split(seq_len(opt$reps), (seq_len(opt$reps) - 1) %% opt$cores)
res = mclapply(chunks, function(rs) {
  acc = list(sw = 0, ssw = 0, su = 0, ssu = 0, less = 0, tie = 0)
  for (r in rs) {
    x = one.rep(r)
    acc$sw = acc$sw + x$w; acc$ssw = acc$ssw + x$w^2
    acc$su = acc$su + x$u; acc$ssu = acc$ssu + x$u^2
    acc$less = acc$less + x$less; acc$tie = acc$tie + x$tie
  }
  acc
}, mc.cores = opt$cores, mc.preschedule = FALSE)
failed = vapply(res, inherits, logical(1), "try-error")
if (any(failed)) stop("A worker failed: ", res[failed][[1]])
tot = Reduce(function(a, b) Map(`+`, a, b), res)

# --- Standardised effect sizes and Raup-Crick ----------------------------------------------
n = opt$reps
ses = function(obs, s, ss) {
  m = s / n; sdv = sqrt((ss - n * m^2) / (n - 1))
  z = (obs - m) / sdv; diag(z) = 0
  dimnames(z) = list(rownames(comm), rownames(comm)); z
}
bNTI.w = ses(obs.w, tot$sw, tot$ssw)
bNTI.u = ses(obs.u, tot$su, tot$ssu)
rc = ((tot$less + tot$tie / 2) / n - 0.5) * 2  # split ties, scaled to [-1, 1]
RC = matrix(0, n.samp, n.samp, dimnames = list(rownames(comm), rownames(comm)))
RC[lower.tri(RC)] = rc                         # vegdist order = lower triangle by columns
RC = RC + t(RC)
dimnames(obs.w) = dimnames(RC)

write.csv(bNTI.w, out("betaNTI_treeholes"), quote = FALSE)
write.csv(bNTI.u, out("betaNTI_unweighted_treeholes"), quote = FALSE)
write.csv(obs.w, out("betaMNTD_treeholes"), quote = FALSE)
write.csv(RC, out("RC_bray_treeholes"), quote = FALSE)

el = format(Sys.time() - t0)
log.file = file.path(opt$outdir, paste0("stegen_null_models", opt$suffix, ".log"))
writeLines(c(paste("Date:", Sys.time()), paste("Elapsed:", el),
             paste0(names(opt), " = ", unlist(opt)),
             paste("Samples:", n.samp, " ASVs:", n.asv),
             capture.output(sessionInfo())), log.file)
cat("Done in", el, "\n")
