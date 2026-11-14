# Niche vs neutral assembly in tree-hole microbial communities

Data and code used for the paper ...

How much of the assembly of tree-hole bacterial communities (16S rRNA amplicons, 10 tree
holes sampled over time) is explained by deterministic (**niche**) and by stochastic
(**neutral**) processes?

---

## Contents

1. [Pipeline at a glance](#pipeline-at-a-glance)
2. [Getting the raw data](#getting-the-raw-data)
3. [Reproducing the analysis (`make`)](#reproducing-the-analysis-make)
4. [Requirements](#requirements)
5. [Scripts (`src/`)](#scripts-src)
6. [External software](#external-software)
7. [Data folders](#data-folders)
8. [Key data files](#key-data-files)
9. [Known issues](#known-issues)

---

## Pipeline at a glance

| Step | Folder | What happens |
|---|---|---|
| 1 | `1_dada2`, `1_samples_properties` | Reads → ASV table; remove chloroplasts, mitochondria, rare ASVs, small samples and controls; rarefaction; phylogeny |
| 2 | `2_Dirichlet` | Samples clustered into *K* metacommunities with Dirichlet multinomial mixtures (K = 10 chosen; K = 7 also explored) |
| 3.1 | `3.1_Neutral_NMGS` | Hubbell's neutral model fitted to each metacommunity → immigration and speciation rates |
| 3.2 | `3.2_Neutral_Stegen` | βNTI / Raup-Crick null models (Stegen et al.) → pairwise assembly processes |
| 4.1 | `4.1_Niche_Envir` | Composition vs environment (PCoA on Aitchison distances, envfit, RDA) |
| 4.2 | `4.2_Niche_PiCRUST*` | Predicted functions (PICRUSt2) |
| 4.3 | `4.3_Niche_Interactions` | Co-occurrence network (Cytoscape session `Net.cys`) |

All code is in `src/`. Each data folder holds the outputs of its step, and most have a `figures/` subfolder.

---

## Getting the raw data

The raw reads are **not** in this repository.

| | |
|---|---|
| **Repository** | European Nucleotide Archive (ENA) |
| **Accession** | [`PRJEB35208`](https://www.ebi.ac.uk/ena/browser/view/PRJEB35208) |
| **Sample metadata** | figshare, [doi:10.6084/m9.figshare.10320713.v1](https://doi.org/10.6084/m9.figshare.10320713.v1) |
| **Where to put them** | `1_dada2/` |
| **Expected names** | `Undetermined_S0_L001_I1_001.fastq.gz`, `Undetermined_S0_L001_R1_001.fastq.gz`, `Undetermined_S0_L001_R2_001.fastq.gz` |

> ⚠️ `dada2_treeholes.R` expects the **multiplexed** Illumina run (index read `I1` plus the
> two paired reads), demultiplexed with `idemp` and a barcode table that is not in the repo.
> If ENA serves one pair of files per sample instead, skip the demultiplexing step and name
> the files `Undetermined_S0_L001_R{1,2}_001.fastq.gz_<sample>`, the names `idemp` would
> have produced.

**Everything dada2 needs goes in `1_dada2/`.** `dada2_treeholes.R` reads from and writes to
that folder (`dirInput`):

| File | Source |
|---|---|
| `Undetermined_S0_L001_{I1,R1,R2}_001.fastq.gz` | raw reads, ENA `PRJEB35208` (see above) |
| `barcodes.txt` | barcode table for `idemp` (one line per sample: barcode and sample ID, tab-separated; each 12-bp barcode is followed by `N`; not in the repo) |
| `silva_nr99_v138.1_train_set.fa.gz` | SILVA training set for `assignTaxonomy`, [Zenodo 4587955](https://doi.org/10.5281/zenodo.4587955) |

`idemp` and `cutadapt` must be in your `PATH` (or edit `idemp <-` / `cutadapt <-` at the top
of the script).

Only the DNA libraries were analysed. The RNA libraries of the same samples were ignored,
because they were almost identical.

---

## Reproducing the analysis (`make`)

A `Makefile` at the repository root chains the scripts. Run every command from the root:

```bash
make help            # list the targets and variables
make -n envfit       # preview what a target would run
make envfit          # run it
make split K=7       # targets that depend on the number of Dirichlet components take K
```

* All R scripts find the repository root with `here::here()` (the `niche_vs_neutral.Rproj`
  file), so they work both from RStudio and with `Rscript`, as long as the working
  directory is inside the repository.
* `make` does not rebuild an output that exists and is newer than its inputs. To rerun a
  step, delete its outputs, or use `make -B`, which also reruns every upstream step.
* Steps marked **[external]** need software that each user installs separately (see
  [External software](#external-software)); their targets print instructions.

| Target | Script | Main outputs |
|---|---|---|
| `dada2` **[external]** | `dada2_treeholes.R` | `1_dada2/03_tabletax/*` |
| `rename` | `rename_otu_table_after_dada2.R` | `1_dada2/otu_table_treeholes_dada2_renamed_tab.txt` |
| `asv_tables` | `filter_rare_ASVs.R` | `1_dada2/otu_table_treeholes_dada2_renamed_clean_red{.tsv,_t.tsv,.csv}` |
| `phyloseq_tables` | `rarefaction_formatting_phyloseq.R` | `1_dada2/*_phyloseq.tsv`, `1_samples_properties/Rarefaction_plot_treeholes.tiff` |
| `repset` | `format_repset_phyloseq.R` | `1_dada2/repset_treeholes_dada2_phyloseq{.txt,_clean.fasta}` |
| `tree` | `phylogeny_treeholes_repset.R` | `3.2_Neutral_Stegen/nj_tree_stegen.tre` |
| `samples_figures` | `create_figure.py` | `1_samples_properties/reads_per_treehole_*.png` |
| `dirichlet` **[external]** | DirichletMixtureGHPFit | `2_Dirichlet/otu_table_treeholes_dada2_Dirichlet_K*.{fit,z,mixture,err}` |
| `split` | `split_otu_table_Dirichlet.R` | `2_Dirichlet/otu_tables_K$K/`, `1_samples_properties/Metadata_singleTree_wDirich-K$K.tsv` |
| `nmgs` **[external]** | NMGS | `3.1_Neutral_NMGS/Results_K$K/*_NMGS{,_s}.csv` |
| `rates_metadata` | `merge_neutral_rates_metadata.R` | `1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv` |
| `neutral_stats` | `extract_neutral_stats.R` | `3.1_Neutral_NMGS/Results_K10/neutral_stats_perDirichlet.tsv` |
| `envfit` | `taxonomy_betadiversity_envfit.R` | `4.1_Niche_Envir/figures/*`, `1_samples_properties/Metadata_singleTree_wDirich-K10_phyloseq.tsv` |
| `immigration` | `immigration_speciation_rates_analysis.R` | `3.1_Neutral_NMGS/figures/*` |
| `stegen` | `stegen_assembly_process_heatmap.R` | `3.2_Neutral_Stegen/figures/AssemblyProcessHeatmap_*` |
| `rda_composition` | `rda_composition_environment.R` | `4.1_Niche_Envir/rda/composition/*`, `4.1_Niche_Envir/rda/results_summary.csv` |
| `rda_rates` | `rda_neutral_rates_environment.R` | `4.1_Niche_Envir/rda/neutral_rates/*`, `4.1_Niche_Envir/rda/results_summary.csv` |
| `picrust` **[external]** | PICRUSt2 | `4.2_Niche_PiCRUST*/` |
| `stamp` | `biom_to_stamp.py` | `4.2_Niche_PiCRUST/KEGG_pathways_out/KEGGs_vs_samples.L3.spf` |

---

## Requirements

* **R** packages
  * Bioconductor (`BiocManager::install(...)`): `dada2`, `ShortRead`, `Biostrings`,
    `phyloseq`, `DECIPHER`, `ALDEx2`
  * CRAN: `here`, `vegan`, `permute`, `nlme`, `compositions`, `ape`, `ggplot2`, `ggrepel`,
    `ggpmisc`, `ggprism`, `patchwork`, `RColorBrewer`, `colorspace`, `plotly`, `Hmisc`,
    `tidyverse` (includes `dplyr`, `tidyr`, `tibble`, `readr`, `purrr`, `stringr`),
    `reshape2`, `data.table`, `stringi`, `rprojroot`, `ggpicrust2`
  * The legacy scripts also use `microbiome`.

  All scripts use `here`. Other packages, by script:

  | Script | Packages |
  |---|---|
  | `dada2_treeholes.R` | `dada2`, `ShortRead`, `Biostrings`, `dplyr`, `tidyr`, `ggplot2`, `plotly`, `Hmisc`, `stringi` |
  | `rename_otu_table_after_dada2.R` | `stringi`, `stringr` |
  | `filter_rare_ASVs.R` | `tidyverse`, `ggplot2` |
  | `rarefaction_formatting_phyloseq.R` | `phyloseq`, `vegan`, `ggplot2`, `stringi`, `stringr` |
  | `format_repset_phyloseq.R` | `phyloseq`, `Biostrings`, `stringi`, `stringr` |
  | `phylogeny_treeholes_repset.R` | `DECIPHER`, `ape`, `Biostrings`, `phyloseq`, `ggplot2`, `stringi`, `stringr` |
  | `split_otu_table_Dirichlet.R` | `phyloseq`, `stringr` |
  | `extract_neutral_stats.R` | `dplyr` |
  | `immigration_speciation_rates_analysis.R` | `phyloseq`, `vegan`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `reshape2`, `data.table`, `stringi`, `stringr` |
  | `stegen_assembly_process_heatmap.R` | `dplyr`, `tidyr`, `reshape2`, `ggplot2` |
  | `taxonomy_betadiversity_envfit.R` | `phyloseq`, `vegan`, `compositions`, `ggpmisc`, `ggplot2`, `RColorBrewer`, `colorspace`, `dplyr`, `tibble`, `stringi`, `stringr` |
  | `rda_composition_environment.R` (+ `rda_functions.R`) | `vegan`, `permute`, `phyloseq`, `ggplot2`, `ggrepel` |
  | `rda_neutral_rates_environment.R` (+ `rda_functions.R`) | `vegan`, `permute`, `nlme`, `phyloseq`, `ggplot2` |
  | `picrust_analysis.R` | `ggpicrust2`, `ALDEx2`, `tidyverse`, `ggprism`, `patchwork`, `data.table`, `rprojroot` |
* **Python 3**: `pandas`, `matplotlib`, `numpy`, and `biom-format` for `biom_to_stamp.py`
* **External tools**: `idemp`, `cutadapt`, a SILVA training set, the Dirichlet multinomial mixtures code, NMGS and PICRUSt2; see [External software](#external-software)

Software versions were not recorded for the published results, except for the sequence
processing (DADA2 v1.32.0, cutadapt v4.6, SILVA v138.1; run on a collaborator's
computer). The RDA analyses (`rda_*.R`) were last run with R 4.6.2, vegan 2.7.6,
permute 0.9.10, nlme 3.1-171 and ggrepel 0.9.8; the previous run (R 4.5.0, vegan 2.7.1)
gave identical results except for the F statistics of the descriptive per-axis test of
the final composition RDA. Alignments and trees in particular depend on the `DECIPHER` version.

---

## Scripts (`src/`)

Scripts are listed in pipeline order. Paths are relative to the repository root.

### 1. Sequence processing

#### `dada2_treeholes.R`
Raw Illumina paired-end reads → amplicon sequence variants (ASVs) with taxonomy. Adapted
from the DADA2 big-data workflow (error rates learned on a subset of reads), but the
sample inference is pooled (`pool = TRUE`). Steps:
1. Demultiplexing with `idemp` using the index reads. The script only checks that `idemp`
   runs; the command itself (`idemp -b barcodes.txt -I1 ... -R1 ... -R2 ... -o 1_dada2/`)
   is given in a comment and must be run before the script.
2. Removal of reads with ambiguous bases (`maxN = 0`).
3. Primer trimming with `cutadapt` (`-n 2`, `--minimum-length 50`): 515F
   `GTGYCAGCMGCCGCGGTAA` (Parada et al. 2016) and 806RB `GGACTACNVGGGTWTCTAAT` (Apprill et al. 2015).
4. Quality filtering: `truncLen = c(120, 220)`, `maxEE = c(2, 2)`, `truncQ = 2`, `rm.phix = TRUE`.
5. Error rates learned from 10⁸ bases of randomly selected reads (`randomize = TRUE`, seed 100).
6. Pooled ASV inference, read merging, and chimera removal (`method = "consensus"`; 98.3% of the reads retained).
7. Taxonomy with the naive Bayesian classifier (`assignTaxonomy`, `tryRC = TRUE`) against SILVA v138.1.

Methods text for the manuscript (DADA2 v1.32.0, cutadapt v4.6; versions given by the
authors, not recorded in the script output): see the *Bioinformatics pipeline* section of
the manuscript.
* **Edit first:** the paths to `idemp`, `cutadapt` and the SILVA training set (`silva.ref`) at the top of the script. Input and output folders (`data.fp`, `project.fp`) are set to `1_dada2/` through `here::here()` and do not need editing.
* **Input:** `1_dada2/Undetermined_S0_L001_{I1,R1,R2}_001.fastq.gz` (see [Getting the raw data](#getting-the-raw-data))
* **Output:** `1_dada2/01_preprocess/`, `1_dada2/02_filter/` (quality plots such as `fwd_qual_plots_filt.png`), `1_dada2/03_tabletax/` (`seqtab_final.rds`, `tax_final.rds`, `repset.fasta`, `seqtab_wTax_mctoolsr.txt`)
* `repset.fasta` and `seqtab_wTax_mctoolsr.txt` were saved as `1_dada2/repset_treeholes_dada2.fasta` and `1_dada2/otu_table_treeholes_dada2.txt`.

#### `rename_otu_table_after_dada2.R`
Replaces the dada2 sample names with those of the metadata, using the fastq file names, so the fastq files must be in `1_dada2/`.
* **Input:** `1_dada2/otu_table_treeholes_dada2.txt`
* **Output:** `1_dada2/otu_table_treeholes_dada2_renamed_tab.txt`

#### Chloroplast and mitochondria removal (shell, in the `Makefile`)
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_tab.txt`
* **Output:** `1_dada2/otu_table_treeholes_dada2_renamed_clean.txt`, `1_dada2/taxonomy.list` (see [OTU tables](#otu-tables))

#### `filter_rare_ASVs.R`
Removes samples with ≤ 5000 reads, ASVs with < 100 reads or present in < 10 samples, and the CONTROL samples.
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean.txt`, `1_dada2/taxonomy.list`, `1_samples_properties/Metadata_singleTree.tsv`
* **Output:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv` (ASVs × samples), `..._red_t.tsv` (samples × ASVs), `..._red.csv` (input of the Dirichlet code)

#### `rarefaction_formatting_phyloseq.R`
Formats the ASV table, taxonomy and metadata for `phyloseq::import_qiime`, draws the rarefaction curves and rarefies to the smallest sample.
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv`, `1_dada2/taxonomy.list`, `1_samples_properties/Metadata_singleTree.tsv`
* **Output:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv`, `1_dada2/Metadata_singleTree_phyloseq.tsv`, `1_samples_properties/Rarefaction_plot_treeholes.tiff`

#### `Rfunctions.R`
Helper functions sourced by other scripts: `ggrare` (rarefaction curves), `top_taxa_abundance`, `multiplot` and phylogenetic-diversity helpers.

#### `format_repset_phyloseq.R`
Strips the taxonomy from the fasta headers and keeps the representative sequences of the ASVs that are not chloroplasts or mitochondria.
* **Input:** `1_dada2/repset_treeholes_dada2.fasta`, phyloseq tables in `1_dada2/`
* **Output:** `1_dada2/repset_treeholes_dada2_phyloseq.txt` (fasta), `1_dada2/repset_treeholes_dada2_phyloseq_clean.fasta`

#### `phylogeny_treeholes_repset.R`
Aligns the representative sequences (`DECIPHER`) and builds a neighbour-joining tree for the βNTI analysis.
* **Input:** `1_dada2/repset_treeholes_dada2_phyloseq.txt`, phyloseq tables in `1_dada2/`
* **Output:** `3.2_Neutral_Stegen/nj_tree_stegen.tre`

#### `create_figure.py` (wrapper: `run_analysis.sh`)
Bubble plots of the number of reads per tree hole and sampling date.
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red_t.tsv`, `1_samples_properties/Metadata_singleTree.tsv`
* **Output:** `1_samples_properties/reads_per_treehole_over_time.png`, `1_samples_properties/reads_per_treehole_evenly_spaced.png`

### 2. Metacommunities

#### `split_otu_table_Dirichlet.R`
Adds each sample's Dirichlet component to the metadata and splits the ASV table into one table per component (raw, transposed and rarefied). The rarefied tables are the NMGS input. Usage: `Rscript src/split_otu_table_Dirichlet.R 10` (number of components, default 10).
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv`, `1_samples_properties/Metadata_singleTree.tsv`, `2_Dirichlet/otu_table_treeholes_dada2_Dirichlet_K$K.z`
* **Output:** `2_Dirichlet/otu_tables_K$K/otu_table_treeholes_dada2_renamed_clean_red_DirichC{i}{.tsv,_t.tsv,_rar.csv}`, `1_dada2/otu_table_treeholes_dada2_renamed_clean_red_rar.csv` (whole table, rarefied), `1_samples_properties/Metadata_singleTree_wDirich-K$K.tsv`

### 3. Neutral processes

#### `merge_neutral_rates_metadata.R`
Adds the NMGS immigration rates (one per sample) and speciation rates (one per metacommunity) to the metadata.
* **Input:** `1_samples_properties/Metadata_singleTree_wDirich-K10.tsv`, `3.1_Neutral_NMGS/Results_K10/{Immigration,Speciation}Rates_Merged_Dirichlet_K10.tsv`
* **Output:** `1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv`

#### `extract_neutral_stats.R`
For each metacommunity: number of samples, quartiles of the immigration rate, and the speciation rate.
* **Input:** `1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv`
* **Output:** `3.1_Neutral_NMGS/Results_K10/neutral_stats_perDirichlet.tsv`

#### `immigration_speciation_rates_analysis.R`
Relates immigration and speciation rates to tree-hole properties (volume, depth:surface ratio, drought), and overlays them on the PCoA of the communities. The PCoA and envfit arrows are read from the outputs of `taxonomy_betadiversity_envfit.R`, so run `make envfit` first.
* **Input:** `1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv`, `1_samples_properties/Metadata_singleTree_wDirich-K10_phyloseq.tsv`, `4.1_Niche_Envir/figures/pcoa_aitchison_{scores,variance}.csv` and `envfit_results_aitchison_pcoa_axes*_2026.csv` (all written by `taxonomy_betadiversity_envfit.R`), `1_dada2/otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv`
* **Output:** `3.1_Neutral_NMGS/figures/{Immigration,Speciation}_*_metacommunity.{pdf,tiff}`, `Depth_surface_ratio_drought_metacommunity.*`, `PCOA_treehole_metacommunity_k10_axes*_{immigr,speciation}.*`

#### `rda_neutral_rates_environment.R`
The immigration and speciation rates as responses of the environment: RDA with forward selection, combinations of tree-hole morphology variables tested by permuting tree holes, time-varying variables tested within tree holes, variation partitioning, a mixed model for immigration (random tree hole) and metacommunity-level correlations for speciation.
It also relates the drought prevalence of each tree hole (maximum number of drought days) to its morphology, to identify the variable that best represents water retention (`drought_vs_morphology.csv`).
* **Input:** `1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv`. The metadata column `Dought_days` (misspelt) is copied to `Drought_days` by `add_predictors()` in `rda_functions.R`.
* **Output:** `4.1_Niche_Envir/rda/neutral_rates/` (`results_summary.csv`, `summary.txt`, tables `.csv`, figures `.pdf`); the manuscript uses `immigration_vs_size_shape_fit.pdf` (Fig. 3B), `variance_explained_combined.pdf` (Fig. 3C; immigration only: `variance_explained_immigration.pdf`), `drought_vs_depth.pdf` (Fig. 3D), and `immigration_vs_size_shape.pdf`, `drought_vs_depth_and_ratio.pdf` and `speciation_metacommunity_top_variable.pdf` (Suppl. Figs.). `immigration_vs_morphology_fit.pdf` adds the model with the highest adjusted R2 next to the a priori one

#### `stegen_assembly_process_heatmap.R`
Classifies each pair of samples as homogeneous or variable selection (|βNTI| > 2), homogenising dispersal or dispersal limitation (|RC| > 0.95), or drift. Plots the pairwise processes as a heatmap, with samples ordered by metacommunity.
* **Input:** `3.2_Neutral_Stegen/betaNTI_treeholes.csv`, `3.2_Neutral_Stegen/RC_bray_treeholes.csv`, `1_samples_properties/Metadata_singleTree_wDirich-K10.tsv`
* **Output:** `3.2_Neutral_Stegen/figures/AssemblyProcessHeatmap_treehole_metacommunities_k10_2026.{pdf,tiff}`

### 4. Niche processes

#### `taxonomy_betadiversity_envfit.R`
Composition barplots (phylum, family, genus) by tree hole and time (rarefied counts); boxplots of tree-hole properties per metacommunity; PCoA on Aitchison distances (unrarefied counts without mitochondria/chloroplasts, CLR within each sample, pseudocount 1) with the significant environmental vectors (envfit); correlations of PCoA axes with depth, volume and depth:surface ratio. Based on J. Ramoneda's `niche_envir_K10_git.R`. Set `dist.method <- "robust.aitchison"` for vegan's robust Aitchison distance; note that vegan then replaces the whole rCLR matrix (not only the zeros) by a rank-3 matrix completion, so the PCoA has only three axes.
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red{,_phyloseq}.tsv`, `1_dada2/taxonomy.list`, `1_samples_properties/Metadata_singleTree_wDirich-K10.tsv`
* **Output:** `4.1_Niche_Envir/figures/*` (figures, `envfit_results_aitchison_pcoa_axes{1_2,1_3,2_3}_2026.csv`, `pcoa_aitchison_{scores,variance}.csv`), `1_samples_properties/Metadata_singleTree_wDirich-K10_phyloseq.tsv`

#### `rda_composition_environment.R`
Which combination of environmental variables best explains composition? Uses RDA on the CLR-transformed community (same data as the envfit) and on the first three Aitchison PCoA axes, with the same steps as `rda_neutral_rates_environment.R`. Depth and surface area enter as the orthogonal pair Size (≈ log volume) and Shape (≈ log depth:surface).
* **Input:** `1_dada2/otu_table_treeholes_dada2_renamed_clean_red_phyloseq.tsv`, `1_samples_properties/Metadata_singleTree_wDirich-K10.tsv`
* **Output:** `4.1_Niche_Envir/rda/composition/` (`results_summary.csv`, `summary.txt`, tables `.csv`, figures `.pdf`); the manuscript uses `size_shape_definition.pdf` and `size_shape_vs_original.pdf` (Suppl. Figs.; coefficients relating Size and Shape to log volume and log depth:surface in `size_shape_vs_original.csv`) and `variance_explained_combined.pdf` (Fig. 2D)

Both RDA scripts also write `4.1_Niche_Envir/rda/results_summary.csv`: one row per test, with the response, predictors, permutation scheme, adjusted R², P, corrected P and the file holding the details.

#### `rda_functions.R`
Helpers shared by the two RDA scripts: data loading, CLR, Size/Shape, fast explained variance, tree-hole permutation test with selection-corrected p-values, results table and figures.

#### `picrust_analysis.R`
Differential abundance of predicted KEGG level-3 pathways between metacommunities with `ggpicrust2` (ALDEx2), PCA of the pathway profiles, and error-bar plots for selected pairs of metacommunities.
* **Input:** `4.2_Niche_PiCRUSTsc/KEGG_pathways_out/KEGGs_vs_samples.L3.tsv` (or `KO_metagenome_out/pred_metagenome_unstrat.tsv` with `ko2kegg = T`), `1_samples_properties/Metadata_singleTree_wDirich-K10.tsv`
* **Output:** `4.2_Niche_PiCRUSTsc/analysis/Plot_ExtErrBar_Met{i}vsMet{j}_selpaths.pdf`
> The pathways shown are chosen by row number (`id_sel`) after visual inspection; re-check them if the input table changes.

#### `biom_to_stamp.py`
Converts a BIOM table into a STAMP profile (Morgan Langille, GPL). For example:
`python3 src/biom_to_stamp.py -m KEGG_Pathways KEGGs_vs_samples.L3.biom > KEGGs_vs_samples.L3.spf`

### One-off and legacy scripts

* `fix_treehole_morphology.R`: corrects the tree-hole morphology of five samples in the metadata (already applied; see [Known issues](#known-issues)).
* `legacy/taxonomy_betadiversity_envfit_rarefiedCLR.R`: previous envfit script (rarefied counts; the CLR of its PCoA figure was computed across samples). Replaced by J. Ramoneda's version, now `taxonomy_betadiversity_envfit.R`.
* `legacy/phyloseq_analysis_private.R`: early exploratory phyloseq analysis from a previous project. It uses old metadata and absolute paths, and is not part of the pipeline.

---

## External software

Users install these tools themselves; this repository contains only their outputs and the commands used.

#### Dirichlet multinomial mixtures (step 2)
Holmes, Harris & Quince (2012) *PLoS ONE* 7:e30126. Run on the filtered table (its first column must be labelled `Taxon`, see [OTU tables](#otu-tables)) with seed 240924:
```bash
cd 2_Dirichlet
seq 1 22 | xargs -i DirichletMixtures/DirichletMixtureGHPFit \
    -in ../1_dada2/otu_table_treeholes_dada2_renamed_clean_red.csv \
    -out otu_table_treeholes_dada2_Dirichlet_K{} -k {} -l 240924
# model evidence for each K (plotted with xmgrace in figures/)
echo '#K Evidence NegativeLogDeterminant LaplaceApproximation BIC AIC' > fit_otu_table_treeholes_dada2_Dirichlet.dat
seq 1 20 | xargs -i cat otu_table_treeholes_dada2_Dirichlet_K{}.fit >> fit_otu_table_treeholes_dada2_Dirichlet.dat
# tree holes per component -> otu_table_treeholes_dada2_Dirichlet_K10.z_Location_*
DirichletMixtures/Scripts/Mixtures2Metadata.pl 10 otu_table_treeholes_dada2_Dirichlet_K10.z ../1_samples_properties/Metadata_singleTree.tsv Location
```

#### NMGS (step 3.1)
Harris et al. (2015) *Proc. IEEE* 105:516-529, <https://github.com/microbiome/NMGS>. The inputs are the rarefied tables, with their first column labelled `OTUs`:
```bash
cd 3.1_Neutral_NMGS
seq 1 10 | xargs -i NMGS/C/NMGS -in ../2_Dirichlet/otu_tables_K10/otu_table_treeholes_dada2_renamed_clean_red_DirichC{}_rar.csv \
    -out Results_K10/otu_table_treeholes_dada2_renamed_clean_red_DirichC{}_rar_NMGS -s -v
NMGS/C/NMGS -in ../1_dada2/otu_table_treeholes_dada2_renamed_clean_red_rar.csv -out Results_K10/otu_table_treeholes_dada2_renamed_clean_red_rar_NMGS -s -v
```
Post-processing uses the scripts in `NMGS/C/Scripts` (edit the paths inside each one):
* `ExtractPseudoPval.sh` / `Sig.pl` → neutrality tests (`*_DirichC-All_logLikePseudPval.txt`)
* `ExtractMedians.sh` / `Median.pl` → `*_ParMedians.txt` (first line: speciation rate; then the immigration rate of each sample)
* `ConcatenateMedians2labels.sh` → `ImmigrationRates_Merged_Dirichlet_K10.tsv`, `SpeciationRates_Merged_Dirichlet_K10.tsv`
* `ExtractMetacommunity.sh` → `*_relAbMetac.txt`

Then run `make rates_metadata`.

#### PICRUSt2 (step 4.2)
```bash
conda activate picrust2
picrust2_pipeline.py -i 1_dada2/otu_table_treeholes_dada2_renamed_clean_red.biom \
    -s 1_dada2/repset_treeholes_dada2_phyloseq_clean.fasta -o 4.2_Niche_PiCRUST -p 10 \
    --in_traits COG,EC,KO,PFAM,TIGRFAM
# PICRUSt2 >= 2.6 (bacteria + archaea references) -> 4.2_Niche_PiCRUSTsc
picrust2_pipeline.py -i 1_dada2/otu_table_treeholes_dada2_renamed_clean_red.biom \
    -s 1_dada2/repset_treeholes_dada2_phyloseq_clean.fasta -o 4.2_Niche_PiCRUSTsc -p 10 \
    --in_traits EC,KO,GO,PFAM,BIGG,CAZY,GENE_NAMES
```
* 40 of the 3537 ASVs exceed the maximum NSTI (2.0) and are dropped; the median NSTI is 0.30.
* KOs are mapped to KEGG level 3 (`KEGGs_vs_samples.L3.tsv`) with PICRUSt's `categorize_by_func.R`.
* Apostrophes were removed from pathway names, and the first column was labelled `KEGG_Pathways`.
* Converted with `biom convert --table-type="Pathway table" --to-hdf5` and `make stamp`.

#### dada2 helpers (step 1)
* **idemp**: demultiplexing of the raw reads with the index reads. Source: <https://github.com/yhwu/idemp> (compile with `make`).
* **cutadapt** (v4.6 used): primer removal. Install with `pip install cutadapt` or `conda install -c bioconda cutadapt`.
* **SILVA v138.1** DADA2-formatted training set `silva_nr99_v138.1_train_set.fa.gz`: <https://doi.org/10.5281/zenodo.4587955>.

Set the paths to these three (`idemp <-`, `cutadapt <-`, `silva.ref <-`) at the top of `dada2_treeholes.R`.

---

## Data folders

| Folder | Contents |
|---|---|
| `1_dada2` | dada2 outputs (quality plots, ASV tables, taxonomy, representative sequences) and the filtered ASV tables. **WILL EXCLUDE FILES AT THE END** |
| `1_samples_properties` | Sample metadata (`Metadata_singleTree*`; `old_metadata/` holds earlier versions) and figures of sample properties (reads per sample, rarefaction curves) |
| `2_Dirichlet` | Dirichlet mixture fits for K = 1–22: `*.fit` (model fit), `*.z` (sample assignments), `*.mixture` (component parameters), `*.err` (logs); evidence plots in `figures/`; ASV tables split by component in `otu_tables_K7/` and `otu_tables_K10/`; `*.z_Location_*` files (how tree holes distribute across components). `2_Dirichlet_noseed_sep2024.tar.bz2` archives a run without a fixed seed |
| `3.1_Neutral_NMGS` | NMGS results in `Results_K7/`, `Results_K10/`: `*_rar_NMGS.csv` (MCMC samples of θ and immigration rates), `*_rar_NMGS_s.csv` (statistics of neutral samples), `*_ParMedians.txt`, `*_relAbMetac.txt`, `*_DirichC-All_logLikePseudPval.txt`, merged rates, `neutral_stats_perDirichlet.tsv`; `figures/` |
| `3.2_Neutral_Stegen` | `betaNTI_treeholes.csv`, `RC_bray_treeholes.csv` (pairwise βNTI and Raup-Crick), `nj_tree_stegen.tre` (ASV tree), `figures/` |
| `4.1_Niche_Envir` | Composition and environment: `figures/` (figures, envfit tables, PCoA scores; `envfit_results_aitchison_pcoa_k10metadata.csv` is an older envfit table), `rda/` (`results_summary.csv`, `composition/`, `neutral_rates/`) |
| `4.2_Niche_PiCRUST` | PICRUSt2, default reference: `*_predicted.tsv.gz`, `*_metagenome_out/`, `pathways_out/`, `KEGG_pathways_out/` (KEGG level 3 in biom, tsv and STAMP formats), `intermediate/`, `figures/` |
| `4.2_Niche_PiCRUSTsc` | PICRUSt2 ≥ 2.6 with bacterial (`bac_*`), archaeal (`arc_*`) and combined (`combined_*`) references and more databases (CAZy, GO, BiGG, gene names); `analysis/` has the differential-pathway plots between metacommunities |
| `4.3_Niche_Interactions` | Co-occurrence network, Cytoscape session `Net.cys` |

---

## Key data files

### OTU tables

Generated in this order (in file names, `_t` = transposed, `_rar` = rarefied):

1. `otu_table_treeholes_dada2.txt`: output of dada2.
2. `otu_table_treeholes_dada2_renamed_tab.txt`: samples renamed as in the metadata (`rename_otu_table_after_dada2.R`).
3. `otu_table_treeholes_dada2_renamed_clean.txt`: chloroplasts and mitochondria removed (`make` targets `1_dada2/taxonomy.list` and `1_dada2/otu_table_treeholes_dada2_renamed_clean.txt`):
   ```bash
   awk -F "\t" '{print $NF}' otu_table_treeholes_dada2_renamed_tab.txt > taxonomy.list
   grep -v hloroplast otu_table_treeholes_dada2_renamed_tab.txt | grep -v itochondria | cut -f 1-214 > otu_table_treeholes_dada2_renamed_clean.txt
   ```
4. `otu_table_treeholes_dada2_renamed_clean_red*`: ASVs seen in < 10 samples or with < 100 reads, samples with ≤ 5000 reads (where the ranked sample sizes change slope) and controls removed.
5. `otu_table_treeholes_dada2_renamed_clean_red_rar*`: rarefied to the smallest sample. The same is done for each component.
6. `otu_table_treeholes_dada2_renamed_clean_red_DirichC{i}*`: the samples of Dirichlet component *i*, in `2_Dirichlet/otu_tables_K$K/`.

**Header of the first column.** The external programs need a label in the first column, so the scripts write it:

| File | First column | Why |
|---|---|---|
| `1_dada2/..._clean_red.{tsv,csv}` | `Taxon` | `DirichletMixtureGHPFit` silently drops one sample if the label is empty |
| `..._rar.csv` (whole table and per component) | `OTUs` | NMGS requires the OTU column to be labelled `OTUs` |

R scripts that read these tables rename the column themselves.

**Rarefaction seed.** All rarefactions use `set.seed(1)` followed by `rarefy_even_depth(..., rngseed = FALSE)`. This reproduces the published tables exactly: the earlier `rngseed = TRUE` silently meant `set.seed(1)` inside phyloseq.

### Metadata tables

| File | Contents |
|---|---|
| `Metadata_singleTree.tsv` | Sample metadata: tree hole (`Location`), date, moisture, light, temperatures, rainfall, drought days, tree-hole surface area (`SA`), `Depth`, `Volume`, orientation (`Deg.North`), distance to the previous tree hole, biomass |
| `Metadata_singleTree.csv` | Same content plus two index columns (`X`, `UID`); tab-separated despite the extension. Redundant, no longer used by any script |
| `Metadata_singleTree_wDirich-K{7,10}.tsv` | Plus the Dirichlet component of each sample (`Dirichlet`) |
| `Metadata_singleTree_wDirich-K10_rates.tsv` | Plus the NMGS immigration (`Immigrate_low/med/high`) and speciation (`Speciation_low/med/high`) rates |
| `Metadata_singleTree_wDirich-K10_phyloseq.tsv` | K10 metadata formatted for `phyloseq::import_qiime` |

`Volume` is not an independent measurement: `Volume = SA × Depth / 3` for every sample.

---

## Known issues

* **Corrected metadata (2026-10-27).** In the original metadata, five samples from 17.02.15
  (plate P1, tree holes 3, 4, 6, 7, 9) had the morphology (`SA`, `Volume`, `Depth`,
  `Distance.from.previous.TH`, `Deg.North`) of another tree hole. `fix_treehole_morphology.R`
  corrected them in all metadata tables.
* **Read-count filter.** `filter_rare_ASVs.R` excludes ASV_1 when counting reads per sample
  (`rowSums(x[,-1])`). Counting it would keep two more samples. The original behaviour is kept so the published tables reproduce.
  ASV_1 itself is present in every count table.
* **Single-sample component.** Dirichlet component C9 (K = 10; C4 for K = 7) has one sample. `split_otu_table_Dirichlet.R` lost the
  ASV names of its tables (rows `1, 2, ...` / `sp1, ...`); fixed with `drop = FALSE` (2026-10-28) and the tables regenerated. The counts did not change.
* **Stegen inputs.** No code in this repo computes `betaNTI_treeholes.csv` and `RC_bray_treeholes.csv`. The RC matrix has 165 samples (including a CONTROL), so it was built from an older table.
* **Tree.** The neighbour-joining tree has the same 3537 tips when regenerated, but a different topology (alignment depends on the `DECIPHER` version).
* **CLR orientation (fixed 2026-10-28).** The previous envfit and immigration scripts applied the CLR per ASV across samples for the PCoA figures. Both now use the Aitchison PCoA computed in `taxonomy_betadiversity_envfit.R`, and all figures in `4.1_Niche_Envir/figures` and the PCoA figures in `3.1_Neutral_NMGS/figures` were regenerated.
