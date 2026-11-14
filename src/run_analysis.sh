#!/bin/bash
# Draw the reads-per-tree-hole figures (1_samples_properties/reads_per_treehole_*.png).
# Requires: pip install pandas matplotlib   (or use `make samples_figures` from the repo root)

# Run the Python script from wherever this file lives
python3 "$(dirname "$0")/create_figure.py"
