###############################################################################
# Makefile -- reproduce the niche vs neutral analysis of the tree-hole data
#
# Run `make` (or `make help`) from the repository root to list the targets.
# All scripts live in src/ and locate the repository root with here::here()
# (the niche_vs_neutral.Rproj file), so they must be run from inside the repo.
#
# Results are already committed, so make will not rebuild an output that
# exists and is newer than its input data. Scripts are order-only
# prerequisites: editing a script does NOT trigger a rebuild. To rerun a step,
# delete its output(s) and call its target; `make -B <target>` also works but
# rebuilds every upstream step as well.
# git does not preserve modification times, and some intermediate files are not
# committed (e.g. 1_dada2/*_phyloseq.tsv), so preview with `make -n <target>`
# before running anything that would overwrite committed results.
#
# Steps marked [external] need software or data that is not in this repo
# (raw reads, Dirichlet multinomial mixtures code, NMGS, PICRUSt2).
###############################################################################

SHELL  := /bin/bash
RSCRIPT ?= Rscript
PYTHON  ?= python3
K       ?= 10# Number of Dirichlet components (metacommunities); 7 and 10 were explored
NMGS    ?= 3.1_Neutral_NMGS/NMGS/C/NMGS# NMGS binary (github.com/microbiome/NMGS)
PICRUST ?= picrust2_pipeline.py
THREADS ?= 4

SRC   := src
DADA  := 1_dada2
META  := 1_samples_properties
DIR   := 2_Dirichlet
NMGSD := 3.1_Neutral_NMGS
STEG  := 3.2_Neutral_Stegen
ENVF  := 4.1_Niche_Envir

OTU      := otu_table_treeholes_dada2_renamed_clean_red
DIRTAB   := $(DIR)/otu_tables_K$(K)
DIRICH_Z := $(DIR)/otu_table_treeholes_dada2_Dirichlet_K$(K).z
COMPS    := $(shell seq 1 $(K))

.PHONY: help all asv_tables phyloseq_tables repset tree samples_figures \
        dirichlet split nmgs rates_metadata neutral_stats immigration stegen envfit rda_composition rda_rates \
        picrust stamp dada2 rename

help:
	@echo "Targets (run in this order; outputs already in the repo are not rebuilt):"
	@echo "  dada2            [external] raw reads -> ASV table (see src/dada2_treeholes.R)"
	@echo "  rename           [external] rename ASV-table samples (needs the fastq files)"
	@echo "  asv_tables       remove chloroplast/mitochondria, rare ASVs and small samples"
	@echo "  phyloseq_tables  phyloseq-formatted tables + rarefaction curve"
	@echo "  repset           representative sequences for phyloseq"
	@echo "  tree             neighbour-joining tree of the ASVs (for bNTI)"
	@echo "  samples_figures  reads per tree hole over time (python)"
	@echo "  dirichlet        [external] Dirichlet multinomial mixtures (Holmes et al. 2012)"
	@echo "  split            split ASV table by Dirichlet component (K=$(K))"
	@echo "  nmgs             [external] fit the neutral model to each component (NMGS)"
	@echo "  rates_metadata   add the NMGS immigration/speciation rates to the metadata"
	@echo "  neutral_stats    summary of immigration/speciation rates per component"
	@echo "  immigration      immigration/speciation rate figures (3.1_Neutral_NMGS/figures)"
	@echo "  stegen           bNTI/Raup-Crick assembly heatmap (3.2_Neutral_Stegen/figures)"
	@echo "  envfit           composition barplots, PCoA, envfit (4.1_Niche_Envir/figures)"
	@echo "  rda_composition  RDA: which environmental variables explain composition (4.1_Niche_Envir/rda/composition)"
	@echo "  rda_rates        RDA/mixed models: environment vs immigration and speciation rates (4.1_Niche_Envir/rda/neutral_rates)"
	@echo "  picrust          [external] PICRUSt2 functional prediction"
	@echo "  stamp            convert the KEGG L3 biom table to STAMP format"
	@echo "  all              every non-external target"
	@echo "Variables: K=$(K) RSCRIPT=$(RSCRIPT) PYTHON=$(PYTHON) NMGS=$(NMGS)"

all: asv_tables phyloseq_tables repset tree samples_figures split rates_metadata neutral_stats \
     envfit immigration stegen rda_composition rda_rates

# --- 1. Sequence processing ----------------------------------------------------
dada2:
	@echo "dada2_treeholes.R needs the raw fastq files, a barcode table, idemp, cutadapt"
	@echo "and the SILVA v138.1 training set, none of which are in the repo (see README,"
	@echo "Getting the raw data). Run idemp first (command in the script header), edit the paths at the top of"
	@echo "$(SRC)/dada2_treeholes.R and run it. Its output $(DADA)/03_tabletax/seqtab_wTax_mctoolsr.txt"
	@echo "was saved as $(DADA)/otu_table_treeholes_dada2.txt."

# Phony on purpose: it needs the fastq file names in 1_dada2/, which are not in the repo
rename: $(DADA)/otu_table_treeholes_dada2.txt | $(SRC)/rename_otu_table_after_dada2.R
	$(RSCRIPT) $(SRC)/rename_otu_table_after_dada2.R

# Remove chloroplasts and mitochondria (taxonomy is the last column)
$(DADA)/taxonomy.list: $(DADA)/otu_table_treeholes_dada2_renamed_tab.txt
	awk -F "\t" '{print $$NF}' $< > $@
$(DADA)/otu_table_treeholes_dada2_renamed_clean.txt: $(DADA)/otu_table_treeholes_dada2_renamed_tab.txt
	ncol=$$(head -1 $< | awk -F "\t" '{print NF-1}'); \
	grep -v hloroplast $< | grep -v itochondria | cut -f 1-$$ncol > $@

# Remove samples < 5000 reads, ASVs < 100 reads or in < 10 samples, and controls
asv_tables: $(DADA)/$(OTU).tsv
$(DADA)/$(OTU).tsv $(DADA)/$(OTU)_t.tsv $(DADA)/$(OTU).csv &: \
		$(DADA)/otu_table_treeholes_dada2_renamed_clean.txt $(DADA)/taxonomy.list \
		$(META)/Metadata_singleTree.tsv | $(SRC)/filter_rare_ASVs.R
	$(RSCRIPT) $(SRC)/filter_rare_ASVs.R

phyloseq_tables: $(DADA)/$(OTU)_phyloseq.tsv
$(DADA)/$(OTU)_phyloseq.tsv $(DADA)/Metadata_singleTree_phyloseq.tsv $(META)/Rarefaction_plot_treeholes.tiff &: \
		$(DADA)/$(OTU).tsv $(DADA)/taxonomy.list $(META)/Metadata_singleTree.tsv \
		| $(SRC)/rarefaction_formatting_phyloseq.R $(SRC)/Rfunctions.R
	$(RSCRIPT) $(SRC)/rarefaction_formatting_phyloseq.R

repset: $(DADA)/repset_treeholes_dada2_phyloseq_clean.fasta
$(DADA)/repset_treeholes_dada2_phyloseq.txt $(DADA)/repset_treeholes_dada2_phyloseq_clean.fasta &: \
		$(DADA)/repset_treeholes_dada2.fasta $(DADA)/$(OTU)_phyloseq.tsv \
		$(DADA)/Metadata_singleTree_phyloseq.tsv | $(SRC)/format_repset_phyloseq.R
	$(RSCRIPT) $(SRC)/format_repset_phyloseq.R

tree: $(STEG)/nj_tree_stegen.tre
$(STEG)/nj_tree_stegen.tre: $(DADA)/repset_treeholes_dada2_phyloseq.txt $(DADA)/$(OTU)_phyloseq.tsv \
		$(DADA)/Metadata_singleTree_phyloseq.tsv | $(SRC)/phylogeny_treeholes_repset.R
	$(RSCRIPT) $(SRC)/phylogeny_treeholes_repset.R

samples_figures: $(META)/reads_per_treehole_over_time.png
$(META)/reads_per_treehole_over_time.png $(META)/reads_per_treehole_evenly_spaced.png &: \
		$(DADA)/$(OTU)_t.tsv $(META)/Metadata_singleTree.tsv | $(SRC)/create_figure.py
	$(PYTHON) $(SRC)/create_figure.py

# --- 2. Dirichlet multinomial mixtures -------------------------------------------
dirichlet: $(DIRICH_Z)
$(DIRICH_Z):
	@test -f $@ || { \
	  echo "[external] $@ is produced by the Dirichlet multinomial mixtures code of"; \
	  echo "Holmes, Harris & Quince (2012), fitted for K = 1..22 on the filtered ASV table."; \
	  echo "The fits are in $(DIR)/ and the model evidence is plotted in $(DIR)/figures/."; \
	  exit 1; }

split: $(DIRTAB)/$(OTU)_DirichC1_rar.csv
$(foreach c,$(COMPS),$(DIRTAB)/$(OTU)_DirichC$(c).tsv $(DIRTAB)/$(OTU)_DirichC$(c)_t.tsv $(DIRTAB)/$(OTU)_DirichC$(c)_rar.csv) \
$(META)/Metadata_singleTree_wDirich-K$(K).tsv $(DADA)/$(OTU)_rar.csv &: \
		$(DADA)/$(OTU).tsv $(META)/Metadata_singleTree.tsv $(DIRICH_Z) | $(SRC)/split_otu_table_Dirichlet.R
	$(RSCRIPT) $(SRC)/split_otu_table_Dirichlet.R $(K)

# --- 3.1 Neutral model fit (NMGS) -------------------------------------------------
NMGS_OUT := $(foreach c,$(COMPS),$(NMGSD)/Results_K$(K)/$(OTU)_DirichC$(c)_rar_NMGS_s.csv)
nmgs: $(NMGS_OUT)
$(NMGSD)/Results_K$(K)/$(OTU)_DirichC%_rar_NMGS_s.csv: $(DIRTAB)/$(OTU)_DirichC%_rar.csv
	@test -x $(NMGS) || { echo "[external] NMGS binary not found: set NMGS=/path/to/NMGS"; exit 1; }
	mkdir -p $(@D)
	$(NMGS) -in $< -out $(NMGSD)/Results_K$(K)/$(OTU)_DirichC$*_rar_NMGS -v -s
# The parameter medians (*_ParMedians.txt), the merged rate tables and the pseudo p-values
# were then extracted with the scripts in NMGS/C/Scripts (ExtractMedians.sh, Median.pl,
# ConcatenateMedians2labels.sh, Sig.pl). The merge of the rates into
# 1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv is done by the rates_metadata target.

rates_metadata: $(META)/Metadata_singleTree_wDirich-K10_rates.tsv
$(META)/Metadata_singleTree_wDirich-K10_rates.tsv: $(META)/Metadata_singleTree_wDirich-K10.tsv \
		$(NMGSD)/Results_K10/ImmigrationRates_Merged_Dirichlet_K10.tsv \
		$(NMGSD)/Results_K10/SpeciationRates_Merged_Dirichlet_K10.tsv | $(SRC)/merge_neutral_rates_metadata.R
	$(RSCRIPT) $(SRC)/merge_neutral_rates_metadata.R

neutral_stats: $(NMGSD)/Results_K10/neutral_stats_perDirichlet.tsv
$(NMGSD)/Results_K10/neutral_stats_perDirichlet.tsv: $(META)/Metadata_singleTree_wDirich-K10_rates.tsv \
		| $(SRC)/extract_neutral_stats.R
	$(RSCRIPT) $(SRC)/extract_neutral_stats.R

# --- 4.1 Environmental drivers (also writes the phyloseq metadata used below) ------
ENVF_FIGS := $(addprefix $(ENVF)/figures/, \
	Barplot_composition_phyla_treeholes_xtime_new.pdf \
	Barplot_composition_family_treeholes_xtime_new.pdf \
	Barplot_composition_genus_nolegend_treeholes_xtime_new.pdf \
	Boxplot_volume_metacommunity.pdf Boxplot_depth_metacommunity.pdf \
	Boxplot_drought_metacommunity.pdf Boxplot_depth_surface_ratio_metacommunity.pdf \
	PCOA_treehole_metacommunity_k10_axes1_2_2026.pdf PCOA_treehole_metacommunity_k10_axes1_3_2026.pdf \
	PCOA_treehole_metacommunity_k10_axes2_3_2026.pdf MDS1_Depth_correlation_2026.pdf \
	MDS3_Depth_surface_ratio_correlation_2026.pdf MDS1_Volume_correlation_2026.pdf)
# PCoA scores and envfit tables, re-used by immigration_speciation_rates_analysis.R
ENVF_TABS := $(addprefix $(ENVF)/figures/, \
	pcoa_aitchison_scores.csv pcoa_aitchison_variance.csv \
	envfit_results_aitchison_pcoa_axes1_2_2026.csv envfit_results_aitchison_pcoa_axes1_3_2026.csv \
	envfit_results_aitchison_pcoa_axes2_3_2026.csv)
envfit: $(ENVF_FIGS) $(ENVF_TABS)
$(ENVF_FIGS) $(ENVF_TABS) $(META)/Metadata_singleTree_wDirich-K10_phyloseq.tsv &: \
		$(DADA)/$(OTU).tsv $(DADA)/$(OTU)_phyloseq.tsv \
		$(DADA)/taxonomy.list $(META)/Metadata_singleTree_wDirich-K10.tsv \
		| $(SRC)/taxonomy_betadiversity_envfit.R
	$(RSCRIPT) $(SRC)/taxonomy_betadiversity_envfit.R

# --- 4.1 Constrained ordination: environment as explanatory variables --------------
# Both targets also update 4.1_Niche_Envir/rda/results_summary.csv (one row per test).
# Set OPENBLAS_NUM_THREADS=1 if the runs are slow (multithreaded BLAS thrashes here).
rda_composition: $(ENVF)/rda/composition/summary.txt
$(ENVF)/rda/composition/summary.txt: $(DADA)/$(OTU)_phyloseq.tsv \
		$(META)/Metadata_singleTree_wDirich-K10.tsv | $(SRC)/rda_composition_environment.R $(SRC)/rda_functions.R
	$(RSCRIPT) $(SRC)/rda_composition_environment.R

rda_rates: $(ENVF)/rda/neutral_rates/summary.txt
$(ENVF)/rda/neutral_rates/summary.txt: $(META)/Metadata_singleTree_wDirich-K10_rates.tsv \
		| $(SRC)/rda_neutral_rates_environment.R $(SRC)/rda_functions.R
	$(RSCRIPT) $(SRC)/rda_neutral_rates_environment.R

# --- 3.1 Immigration and speciation rates vs environment -----------------------------
IMM_FIGS := $(addprefix $(NMGSD)/figures/, \
	Immigration_depth_surface_ratio_metacommunity.pdf Immigration_drought_metacommunity.pdf \
	Immigration_volume_metacommunity.pdf Depth_surface_ratio_drought_metacommunity.pdf \
	PCOA_treehole_metacommunity_k10_axes1_2_immigr.pdf PCOA_treehole_metacommunity_k10_axes1_2_speciation.pdf \
	PCOA_treehole_metacommunity_k10_axes1_3_immigr.pdf PCOA_treehole_metacommunity_k10_axes1_3_speciation.pdf \
	Speciation_depth_surface_ratio_metacommunity.pdf Speciation_drought_metacommunity.pdf \
	Speciation_volume_metacommunity.pdf)
immigration: $(IMM_FIGS)
$(IMM_FIGS) &: $(META)/Metadata_singleTree_wDirich-K10_rates.tsv \
		$(META)/Metadata_singleTree_wDirich-K10_phyloseq.tsv $(DADA)/$(OTU)_phyloseq.tsv $(ENVF_TABS) \
		| $(SRC)/immigration_speciation_rates_analysis.R
	$(RSCRIPT) $(SRC)/immigration_speciation_rates_analysis.R

# --- 3.2 Stegen et al. null models ------------------------------------------------
stegen: $(STEG)/figures/AssemblyProcessHeatmap_treehole_metacommunities_k10_2026.pdf
$(STEG)/figures/AssemblyProcessHeatmap_treehole_metacommunities_k10_2026.pdf: \
		$(STEG)/betaNTI_treeholes.csv $(STEG)/RC_bray_treeholes.csv \
		$(META)/Metadata_singleTree_wDirich-K10.tsv | $(SRC)/stegen_assembly_process_heatmap.R
	$(RSCRIPT) $(SRC)/stegen_assembly_process_heatmap.R

# --- 4.2 Functional prediction ------------------------------------------------------
picrust:
	@echo "[external] PICRUSt2 was run outside this repo; outputs are in 4.2_Niche_PiCRUST"
	@echo "(default reference) and 4.2_Niche_PiCRUSTsc (bacteria + archaea references)."
	@echo "Reconstructed command (check inputs before use):"
	@echo "  $(PICRUST) -s $(DADA)/repset_treeholes_dada2_phyloseq_clean.fasta \\"
	@echo "     -i $(DADA)/$(OTU).biom -o 4.2_Niche_PiCRUST -p $(THREADS)"

stamp: 4.2_Niche_PiCRUST/KEGG_pathways_out/KEGGs_vs_samples.L3.spf
4.2_Niche_PiCRUST/KEGG_pathways_out/KEGGs_vs_samples.L3.spf: 4.2_Niche_PiCRUST/KEGG_pathways_out/KEGGs_vs_samples.L3.biom \
		| $(SRC)/biom_to_stamp.py
	$(PYTHON) $(SRC)/biom_to_stamp.py -m KEGG_Pathways $< > $@
