#!/usr/bin/env bash
###############################################################################
# SparCC networks with FastSpar for every metacommunity, tree hole and the
# whole table, followed by the edge list of the significant links (p <= 0.001).
#
# Needs fastspar, fastspar_bootstrap and fastspar_pvalues
# (https://github.com/scwatts/fastspar) and GNU parallel in the PATH.
# Usage, from the repository root:
#   bash src/sparcc_run_fastspar.sh                         # all the clusters
#   bash src/sparcc_run_fastspar.sh Meta_7 Hole_1 Hole_7    # some of them
#   EDGELIST_ONLY=1 bash src/sparcc_run_fastspar.sh         # only redo the edge lists
#
# Input:  4.3_Niche_Interactions/SparCC/<cluster>/data.tsv  (src/sparcc_split_tables.R)
# Output: 4.3_Niche_Interactions/SparCC/<cluster>/{cor.tsv,pvalues.tsv,correlations.edgelist,fastspar.log}
#
# Same commands as the original 1.FastSpar.sh / 1.FastSpar_pval.sh; only the paths
# changed. Takes several hours per cluster (1000 bootstraps).
###############################################################################
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
SPARCC="$ROOT/4.3_Niche_Interactions/SparCC"
THREADS=${THREADS:-30}

if [ $# -gt 0 ]; then clusters=("$@"); else
  clusters=(); for d in "$SPARCC"/Hole_* "$SPARCC"/Meta_* "$SPARCC"/Total; do clusters+=("$(basename "$d")"); done
fi

for i in "${clusters[@]}"; do
  cd "$SPARCC/$i"
  if [ -z "${EDGELIST_ONLY:-}" ]; then
    fastspar --otu_table data.tsv --correlation cor.tsv --covariance cov.tsv
    mkdir -p bootstraps bootcor
    fastspar_bootstrap --otu_table data.tsv --number 1000 --prefix bootstraps/fake_data
    parallel -j 10 fastspar --otu_table {} --correlation bootcor/cor_{/} --covariance bootcor/cov_{/} -i 5 --threads 3 ::: bootstraps/*
    fastspar_pvalues --otu_table data.tsv --correlation cor.tsv --prefix bootcor/cor_fake_data_ \
        --permutations 1000 --outfile pvalues.tsv --threads "$THREADS"
    rm -r bootstraps bootcor cov.tsv
  fi
  Rscript "$ROOT/src/sparcc_edgelist.R" "$SPARCC/$i" > fastspar.log 2>&1
  echo "$i"
done
