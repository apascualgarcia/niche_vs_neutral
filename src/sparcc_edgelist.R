###############################################################################
# FastSpar output -> edge list of significant correlations
#
# Keeps the pairs with a bootstrap p-value <= 0.001 (1000 permutations).
# Usage: Rscript src/sparcc_edgelist.R 4.3_Niche_Interactions/SparCC/Meta_7
#
# Input:
#   <dir>/cor.tsv, <dir>/pvalues.tsv (from src/sparcc_run_fastspar.sh)
# Output:
#   <dir>/correlations.edgelist (Sp1, Sp2, cor, pvalues; both orientations of each pair)
# Run from anywhere inside the repository: Rscript src/sparcc_edgelist.R
###############################################################################
library(here)

library(dplyr)

args = commandArgs(trailingOnly = TRUE)

dir = args[1]                                  # SparCC folder of one cluster
out = file.path(dir, "correlations.edgelist")

pvalues <- read.table(file.path(dir, "pvalues.tsv"), header = FALSE, sep = "\t", row.names = 1)
colnames(pvalues) = rownames(pvalues)
cor <- read.table(file.path(dir, "cor.tsv"), header = FALSE, sep = "\t", row.names = 1)
colnames(cor) = rownames(cor)

coords <- which(pvalues <= 0.001, arr.ind = TRUE)

row_names <- rownames(cor)[coords[, 1]]
col_names  <- colnames(cor)[coords[, 2]]

abs_values <- cor[pvalues <= 0.001]
pval <- pvalues[pvalues <= 0.001]


df <- data.frame(
  Row = row_names,
  Column = col_names,
  Correlation = abs_values,
  pvalues = pval
)

colnames(df) = c("Sp1","Sp2", "cor", "pvalues" )


write.table(df, out, quote = FALSE, row.names = FALSE, sep = "\t")


