library(tidyverse)

remove_bidirectional_duplicates <- function(df) {
  df %>%
    mutate(
      # pmin and pmax compare the elements pairwise and return the alphabetically smaller/larger one
      t1_temp = pmin(taxa_1, taxa_2),
      t2_temp = pmax(taxa_1, taxa_2)
    ) %>%
    # Overwrite the original columns with the standardised order
    mutate(
      taxa_1 = t1_temp,
      taxa_2 = t2_temp
    ) %>%
    # Drop the temporary columns
    select(-t1_temp, -t2_temp) %>%
    # distinct removes the duplicated rows based on the taxa and the metacommunity.
    # .keep_all = TRUE keeps the 'value' column
    distinct(taxa_1, taxa_2, MetaCom, .keep_all = TRUE)
}

cors_met <- read.table("../../6.Correlation_Distribution/Results/Correlation_zscore_thres.tsv", header = T, sep = "\t")
cors_holes <- read.table("../../6.Correlation_Distribution/Results/Correlation_zscore_thres_holes.tsv", header = T, sep = "\t")


cors_met_7 <- cors_met %>% filter(MetaCom == "Meta_7" & method == "SparCC")


a <- hist(cors_met_7$zscore, breaks = 1000)

## The second one starts at around 3


view(data.frame(counts = a$counts, mids = a$mids))
