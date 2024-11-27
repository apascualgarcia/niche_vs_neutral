library(tidyverse)

# Function to extract values from LOCAL MODEL section
extract_local_model_values <- function(relabmetac_file) {
  # Read lines from file
  lines <- readLines(relabmetac_file)
  
  # Find where Local Model section starts (line 5)
  local_model_start <- which(grepl("#LocalModel", lines))
  
  if (length(local_model_start) == 0) {
    stop(paste("No Local Model section found in", relabmetac_file))
  }
  
  # Find where Full Model section starts (this marks the end of Local Model section)
  full_model_start <- which(grepl("#FullModel", lines))
  
  if (length(full_model_start) == 0) {
    stop(paste("No Full Model section found in", relabmetac_file))
  }
  
  # Extract values between Local Model and Full Model sections
  values <- lines[(local_model_start + 1):(full_model_start - 1)]
  values <- values[values != ""]  # Remove empty lines
  
  # Split each line and take second value (abundance)
  abundances <- sapply(values, function(x) {
    parts <- strsplit(trimws(x), "\\s+")[[1]]
    if (length(parts) >= 2) as.numeric(parts[2]) else NA
  })
  
  # Remove any NAs that might have come from malformed lines
  abundances <- abundances[!is.na(abundances)]
  
  return(abundances)
}

# Function to get ASV IDs from OTU table
get_asv_ids <- function(otu_file) {
  otu_table <- read.csv(otu_file)
  return(otu_table[,1])  # First column contains ASV IDs
}

# Main function to process all metacommunities
process_metacommunities <- function(k7_results_dir, k7_otu_dir) {
  # Initialize list to store metacommunity abundances
  metacom_abundances <- list()
  
  # Process each metacommunity
  for (i in 1:7) {  # K7 has 7 metacommunities
    # Construct file paths
    relabmetac_file <- file.path(k7_results_dir, 
      sprintf("otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_relAbMetac.txt", i))
    otu_file <- file.path(k7_otu_dir,
      sprintf("otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_rar.csv", i))
    
    # Get ASV IDs and abundance values
    asv_ids <- get_asv_ids(otu_file)
    abundances <- extract_local_model_values(relabmetac_file)
    
    # Check lengths match
    # if (length(asv_ids) != length(abundances)) {
    #   stop(sprintf("Mismatch in number of ASVs (%d) and abundances (%d) for metacommunity %d",
    #               length(asv_ids), length(abundances), i))
    # }
    
    # Store in list
    metacom_abundances[[sprintf("C%d", i)]] <- abundances
    
    # Store ASV IDs from first metacommunity
    if (i == 1) {
      all_asv_ids <- asv_ids
    }
  }
  
  # Convert list to data frame
  df <- data.frame(ASV = all_asv_ids)
  for (metacom in names(metacom_abundances)) {
    df[[metacom]] <- metacom_abundances[[metacom]]
  }
  
  # Save to file
  # write.csv(df, output_file, row.names = FALSE)
  # cat(sprintf("Saved combined metacommunity abundances to %s\n", output_file))
  
  return(df)
}

# Execute if running as script
if (!interactive()) {
  # Paths
  k7_results_dir <- "3.1_Neutral_NMGS/Results_K7"
  k7_otu_dir <- "2_Dirichlet/otu_tables_K7"
  
  # Process metacommunities
  df <- process_metacommunities(k7_results_dir, k7_otu_dir)
  cat("\nFirst few rows of the combined table:\n")
  print(head(df))
}
