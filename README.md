

Data and code used for the paper ...

### Overall description

The repository contains a folder `/src` with the scripts used to analsye the data and
several data folders ordered according to the pipeline followed in the manuscript. Each
folder contains a subfolder called `figures` for plots.


### Code folder

**Scripts should contain relative paths such that any user can just clone the repo
and use the scripts. Please remove any absolute paths**

The scripts provided are:

* `dada2_treeholes.R`: **Josep add description**
* `rename_otu_table_after_dada2.R` : **Josep add description**
* `filter_rare_ASVs.R`: This script takes the ASVs table generated through the dada2 pipeline and, after removing samples with a low number of reads, filter ASVs with a low
number of reads and/or obverved in few samples.
* `split_otu_table_Dirichlet.R`: It  includes the factor identifying the Dirichlet components in the metadata table, and then split the otu table into subtables, to run NMGS.
* ``

External code published by other authors is described in the correspondent data folder.


### Data folders

The data folders contain:

* `1_dada2`: dada2 pipeline ... **Josep please describe data contained and data excluded. Since the starting data are gigantic, we will point to a repo in which data can be downloaded such that the user can download it and run the pipeline.**
* `1_samples_properties`: Different files describing the experiment's metadata and figures describing samples properties.
* `2_Dirichlet`: Output of the Dirichlet mixtures clustering. Code to generate this data is available from the reference
> Holmes, Ian, Keith Harris, and Christopher Quince. "Dirichlet multinomial mixtures: generative models for microbial metagenomics." PloS one 7.2 (2012): e30126.
* `3.1_Neutral_NMGS`: Output of the Hierarchical Dirichlet fit. Code to generate this data is available at:
> Harris, Keith, et al. "Linking statistical and ecological theory: Hubbell's unified neutral theory of biodiversity as a hierarchical Dirichlet process." Proceedings of the IEEE 105.3 (2015): 516-529.
* `3.2_Neutral_Stegen`: Analysis of selection and neutral effects following the proposal of Stegen et al.
* `4.1_Niche_Envir`: Analysis relating environmental variables with composition.
* `4.2_PiCRUST`: Output f PiCRUST 2 analysis.


### Key data files

#### OTU tables.

The following otu tables were generated in the order they appear:

* `otu_table_treeholes_dada2.csv`: OTU table, output of dada2. The label `_t` in any `otu_table` file stands for "transposed" table.
	* `otu_table_treeholes_dada2_renamed_tab.txt ` Table with the names as they appear in the metadata. Obtained with the script `Rename_otu_table_after_dada2.R`.
	* `otu_table_treeholes_dada2_clean.csv`:  The previous table after removing Chloroplasts and Mitochondria as follows:
		```bash
	         $ cat otu_table_treeholes_dada2_renamed_tab.txt | awk -F "\t" '{print $NF}' # 215 columns          
		 $ cat otu_table_treeholes_dada2_renamed_tab.txt | awk -F "\t" '{print $NF}' > taxonomy.list
		 $ cat otu_table_treeholes_dada2_renamed_tab.txt | grep -v hloroplast | grep -v itochondria | cut -f 1-214 > otu_table_treeholes_dada2_renamed_clean.txt
		```
	* `otu_table_treeholes_dada2_renamed_clean_red$label` Tables in which rare otus (observed in less than 10 samples or less than 100 reads) and samples with less than 5K reads, were removed. The latter threshold was determined after observing a change in the slope of the samples' rank at 5K.
	* `otu_table_treeholes_dada2_renamed_clean_red_rar$label` Rarefied tables. Rarefaction was performed considering the minimum number of reads in the table (same for subtables). This should not be much smaller than 5K (could be smaller since we remove rare otus after keeping 5K samples).
	* `otu_table_treeholes_dada2_renamed_clean_red_DirichC$label` Split of the otu table into the subset of samples belonging to Dirichlet Component described in label, eventually with rarefaction (label `rar`). Located in directory `2_Dirichlet`.
