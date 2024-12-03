source("src/immigration_rates.R")
source("src/process_metacommunities.R")

library(here)
library(dplyr)
library(ggplot2)

library(Biostrings)
library(zCompositions)
library(phyloseq)
library(ape)
library(vegan)



samplefile <- here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv")
otufile <- here("1_dada2", "otu_table_treeholes_dada2_renamed_clean_red.tsv")
taxfile <- here("1_dada2",  "taxonomy.list")
repsetfile <- here("1_dada2", "repset_treeholes_dada2_phyloseq_clean.fasta")
treefile <- here("3.2_Neutral_Stegen", "figures", "nj_tree_stegen.tre")

# --- Load OTUs
otu.in <- read.csv(otufile, sep = "\t", row.names = 1)
otu.pseq <- otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)

meta <- import_qiime_sample_data(samplefile)
meta <- meta[meta$Type == "DNA"]

# --- Load taxonomy and clean it 
col_names <- c("Kingdom","Phylum",
                    "Class","Order","Family","Genus","taxa_id")
tax.in <- read.table(taxfile, sep = ";", skip = 1, col.names = col_names)
rownames(tax.in) <- tax.in$taxa_id
tax.in <- subset(tax.in, select = -c(taxa_id))
tax.matrix <- as.matrix(tax.in)
tax.pseq <- tax_table(tax.matrix)

refseq_data <- readDNAStringSet(repsetfile)

phylo_tree <- read.tree(treefile)


physeq <- phyloseq(otu_table = otu.pseq, sample_data = meta,
                    tax_table = tax.pseq, phy_tree = phy_tree(phylo_tree),
                    refseq = refseq_data)

nsamples(physeq)
ntaxa(physeq)
sample_variables(physeq)
rank_names(physeq)

### First task
# Read your data
parmedian_dir <- "3.1_Neutral_NMGS/Results_K10"
otu_dir <- "2_Dirichlet/otu_tables_K10"

## Get immigration rates
imm_rates <- get_immigration_rates(parmedian_dir = parmedian_dir, dirichlet_otu_dir = otu_dir)
rownames(imm_rates)

# Extract the current sample_data
current_sample_data <- data.frame(sample_data(physeq))

# Merge the data frames by row names
merged_sample_data <- merge(current_sample_data, imm_rates, by = "row.names", all.x = TRUE)

# Set row names to SampleID and remove the extra column
rownames(merged_sample_data) <- merged_sample_data$Row.names
merged_sample_data$Row.names <- NULL

# Update the sample_data in the phyloseq object
sample_data(physeq) <- sample_data(merged_sample_data)

## Get speciation rate

sp_rates <- get_metacommunity_rates(parmedian_dir = parmedian_dir)
rownames(sp_rates) <- sp_rates$metacommunity
sp_rates$metacommunity <- NULL

current_sample_data <- data.frame(sample_data(physeq))

merged_df <- merge(data.frame(sample_data(physeq)), sp_rates, by = 0, all = TRUE,
suffixes = c(".sample", ".metacomm"))

# Combine 'rate' columns
meta_df <- merged_df %>%
  mutate(rate = coalesce(rate.sample, rate.metacomm)) %>%
  dplyr::select(-rate.sample, -rate.metacomm)

rownames(meta_df) <- meta_df$Row.names
meta_df$Row.names <- NULL

# Update the sample_data in the phyloseq object
sample_data(physeq) <- sample_data(meta_df)

write.table(meta_df, sep = '\t', row.names = TRUE, col.names = TRUE,
"1_samples_properties/Metadata_singleTree_wDirich-K10_wNMGS.tsv")


### Process metacommunities

# Manual test, check if the files match on number of ASV's
# k7 C2
# otu_table: 2370 ASV's
# Local: 5-2375 = 2370
# Full:2377-4963 = 2586

# k10 C1
# Local 3131
# otu table:2489

# k10 C2
# Local: 2370
# otu table: 2764

# k10 C9
# Local: 985
# otu table: 985

# Paths
k10_results_dir <- "3.1_Neutral_NMGS/Results_K10"
k10_otu_dir <- "2_Dirichlet/otu_tables_K10"

# Process metacommunities
df2 <- process_metacommunities(k10_results_dir, k10_otu_dir, 10)

df1 <- data.frame(otu_table(physeq))
df1 <- transform(df1, rowname = rownames(df1))

df2 <- df2 %>% dplyr::rename(rowname = ASV)

# Merge based on the rowname column
df_merged <- merge(df1, df2, by = "rowname", all = TRUE)  # all=TRUE for a full join

# Set the row names back to the merged data frame (optional)
rownames(df_merged) <- df_merged$rowname
df_merged$rowname <- NULL

outfp = "/Users/theco/Projects/niche_vs_neutral/1_dada2/otu_table_treeholes_dada2_renamed_clean_red_rar_metacom.csv"
write.csv(df_merged, outfp)

# Update the merged data in the phyloseq object
# otu_table(physeq) <- otu_table(as.matrix(df_merged), taxa_are_rows = TRUE)


### Ordination
otu.czm <- cmultRepl(t(otu_table(physeq)), method="CZM", output="p-counts", z.delete = FALSE)

zpat <- zPatterns(t(otu_table(physeq)), label = 0, type = "cgm", plot = TRUE)

nmds.res <- metaMDS(otu.czm, k = 2, plot = FALSE,
                    distance = "robust.aitchison", try = 100)

# Plotting NMDS
nmds.df <- as.data.frame(scores(nmds.res$points))

nmds.df <- merge(sample_data(physeq), nmds.df, by = 0)
rownames(nmds.df) <- nmds.df$Row.names
nmds.df$Row.names <- NULL

ggplot(nmds.df, aes(x = MDS1, y = MDS2, color = Location, size = immigration_rate)) +
    geom_point(alpha = 0.6) + 
    scale_size_continuous(range = c(0.1, 20)) +
    theme_minimal()
