###################################
#  merge_neutral_rates_metadata.R
##################################
# Adds the immigration and speciation rates estimated with NMGS for each
# Dirichlet component (metacommunity) to the sample metadata.
#
# Immigration rates are estimated per sample, speciation rates per
# metacommunity, so the latter are repeated for every sample of a component.
# The low/med/high columns are the lower confidence limit, median and upper
# confidence limit of the posterior (NMGS/C/Scripts/ExtractMedians.sh and
# ConcatenateMedians2labels.pl).
#
# Input:  1_samples_properties/Metadata_singleTree_wDirich-K10.tsv
#         3.1_Neutral_NMGS/Results_K10/ImmigrationRates_Merged_Dirichlet_K10.tsv
#         3.1_Neutral_NMGS/Results_K10/SpeciationRates_Merged_Dirichlet_K10.tsv
# Output: 1_samples_properties/Metadata_singleTree_wDirich-K10_rates.tsv
#         (or the path given as first argument)
###################################

library(here)

args = commandArgs(trailingOnly = TRUE)
file.out = if (length(args) > 0) args[1] else here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv")

# The merged rate tables mix tabs and spaces, and start with comment lines
read_rates = function(file) {
  read.table(file, header = TRUE, comment.char = "#", sep = "", stringsAsFactors = FALSE)
}

meta = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"),
                  header = TRUE, sep = "\t")
meta = meta[, !grepl("^X(\\.\\d+)?$", names(meta))] # empty column from the trailing tab
imm = read_rates(here("3.1_Neutral_NMGS", "Results_K10", "ImmigrationRates_Merged_Dirichlet_K10.tsv"))
names(imm)[names(imm) == "sample_id"] = "ID.Sample"
spe = read_rates(here("3.1_Neutral_NMGS", "Results_K10", "SpeciationRates_Merged_Dirichlet_K10.tsv"))

stopifnot(all(imm$ID.Sample %in% meta$ID.Sample))
rates = merge(meta, imm, by = c("Dirichlet", "ID.Sample"))
rates = merge(rates, spe, by = "Dirichlet")
cat("Samples with rates:", nrow(rates), "of", sum(!is.na(meta$Dirichlet)), "assigned to a component\n")

write.table(rates, file = file.out, sep = "\t", quote = FALSE, col.names = NA)
cat("Written:", file.out, "\n")
