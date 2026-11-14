rm(list=ls()) # clean/reset environment

library(dplyr)
library(here)

# 1 Load metadata with immigration and speciation rates
meta_treeholes = read.table(here("1_samples_properties", "Metadata_singleTree_wDirich-K10_rates.tsv"), header = T)

# 2 Split the samples by Dirichlet component ("sets") and summarize:
#   - number of samples in the set
#   - interquartile values of Immigrate_med (Q1, median, Q3 and IQR)
#   - speciation rates, which are constant within each set
neutral.stats = meta_treeholes %>%
  group_by(Dirichlet) %>%
  summarise(
    N = n(),
    Immigrate_med_Q1 = quantile(Immigrate_med, 0.25, na.rm = TRUE),
    Immigrate_med_Q2 = quantile(Immigrate_med, 0.50, na.rm = TRUE),
    Immigrate_med_Q3 = quantile(Immigrate_med, 0.75, na.rm = TRUE),
    Immigrate_med_IQR = Immigrate_med_Q3 - Immigrate_med_Q1,
    Speciation_low = unique(Speciation_low),
    Speciation_med = unique(Speciation_med),
    Speciation_high = unique(Speciation_high),
    .groups = "drop"
  ) %>%
  arrange(Dirichlet)

# Sanity check: speciation rates must be constant within each set
const.check = meta_treeholes %>%
  group_by(Dirichlet) %>%
  summarise(across(c(Speciation_low, Speciation_med, Speciation_high),
                   ~ length(unique(.x))), .groups = "drop")
if (any(const.check[, -1] != 1)) {
  warning("Speciation rates are not constant within all Dirichlet sets")
}

print(as.data.frame(neutral.stats))

# 3 Write the results
out.file = here("3.1_Neutral_NMGS", "Results_K10", "neutral_stats_perDirichlet.tsv")
write.table(neutral.stats, file = out.file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Results written to:", out.file, "\n")
