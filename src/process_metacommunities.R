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
  
  # Remove the last element (other species) and renormalize
  if (length(abundances) > 1) {
    abundances <- abundances[-length(abundances)]  # Remove last element
    abundances <- abundances / sum(abundances)     # Renormalize to sum to 1
  }
  
  return(abundances)
}

# Function to get ASV IDs from OTU table
get_asv_ids <- function(otu_file) {
  otu_table <- read.csv(otu_file)
  return(otu_table[,1])  # First column contains ASV IDs
}

# Main function to process all metacommunities
process_metacommunities <- function(results_dir, otu_dir, n_metacom) {
  # Initialize list to store all ASV IDs and abundances
  all_asv_ids <- character()
  metacom_abundances <- list()
  
  # First pass: collect all unique ASV IDs and store abundances
  for (i in 1:n_metacom) {
    # Construct file paths
    relabmetac_file <- file.path(results_dir, 
      sprintf("otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_relAbMetac.txt", i))
    otu_file <- file.path(otu_dir,
      sprintf("otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_rar.csv", i))
    
    # Get ASV IDs and abundance values
    asv_ids <- get_asv_ids(otu_file)
    abundances <- extract_local_model_values(relabmetac_file)
    
    # Store abundances with their ASV IDs
    names(abundances) <- asv_ids[1:length(abundances)]
    metacom_abundances[[sprintf("C%d", i)]] <- abundances
    
    # Add to master list of ASV IDs
    all_asv_ids <- union(all_asv_ids, asv_ids[1:length(abundances)])
  }
  
  # Create data frame with all ASV IDs
  df <- data.frame(ASV = all_asv_ids)
  
  # Second pass: fill in abundances for each metacommunity
  for (metacom in names(metacom_abundances)) {
    abundances <- metacom_abundances[[metacom]]
    
    # Create a vector of zeros for all ASVs
    metacom_col <- rep(0, length(all_asv_ids))
    names(metacom_col) <- all_asv_ids
    
    # Fill in the actual abundances where we have them
    metacom_col[names(abundances)] <- abundances
    
    # Add to dataframe
    df[[metacom]] <- metacom_col
  }
  
  return(df)
}

# Execute if running as script
if (!interactive()) {
  # Paths
  k10_results_dir <- "3.1_Neutral_NMGS/Results_K10"
  k10_otu_dir <- "2_Dirichlet/otu_tables_K10"
  
  # Process metacommunities
  df <- process_metacommunities(k10_results_dir, k10_otu_dir, 10)
  cat("\nFirst few rows of the combined table:\n")
  print(head(df))
}
