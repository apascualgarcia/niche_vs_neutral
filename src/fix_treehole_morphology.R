###################################
#  fix_treehole_morphology.R
##################################
# One-off correction of the sample metadata (run once, already applied).
#
# Tree-hole morphology (SA, Volume, Depth, Distance.from.previous.TH,
# Deg.North) must be the same for all samples of a tree hole (Location).
# In the original metadata five samples (all from 17.02.15, plate P1) carry
# the morphology of another tree hole, as if those columns had been shifted by
# one row. Metadata_singleTree_wDirich-K10_rates.tsv (J. Ramoneda, 2026)
# already had the corrected values.
#
# For every sample whose morphology differs from the most common value of its
# tree hole, the five fields are replaced by the text of that most common value.
# The files are edited as text, field by field, so nothing else changes.
#
# Usage: Rscript src/fix_treehole_morphology.R file1.tsv [file2.tsv ...]
###################################

morph.cols = c("SA", "Volume", "Depth", "Distance.from.previous.TH", "Deg.North")

fix_file = function(file) {
  lines = readLines(file)
  eol = ifelse(grepl("\r$", lines), "\r", "")
  fields = strsplit(sub("\r$", "", lines), "\t", fixed = TRUE)
  header = fields[[1]]
  col.loc = match("Location", header)
  col.m = match(morph.cols, header)
  if (anyNA(c(col.loc, col.m))) stop("Missing columns in ", file)
  rows = 2:length(fields)
  loc = sapply(fields[rows], `[`, col.loc)
  key = sapply(fields[rows], function(f) paste(f[col.m], collapse = "\t"))
  ok = loc %in% as.character(1:10)
  mode.key = tapply(key[ok], loc[ok], function(k) names(which.max(table(k))))
  bad = which(ok & key != mode.key[loc])
  for (i in bad) {
    f = fields[[rows[i]]]
    f[col.m] = strsplit(mode.key[[loc[i]]], "\t", fixed = TRUE)[[1]]
    fields[[rows[i]]] = f
  }
  if (length(bad)) {
    out = paste0(sapply(fields, paste, collapse = "\t"), eol)
    writeLines(out, file, sep = "\n")
  }
  cat(sprintf("%s: %d sample(s) corrected\n", file, length(bad)))
  if (length(bad)) print(data.frame(Row = rows[bad], Location = loc[bad]), row.names = FALSE)
}

for (f in commandArgs(trailingOnly = TRUE)) fix_file(f)
