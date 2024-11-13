`library(tidyverse)

get_immigration_rates <- function(parmedian_dir, dirichlet_otu_dir) {
  # Get all ParMedian files
  parmedian_files <- list.files(
    path = parmedian_dir,
    pattern = "ParMedians.txt$",
    full.names = TRUE
  )
  
  # Initialize empty dataframe for results
  all_rates <- data.frame()
  
  # Process each ParMedian file
  for (pfile in parmedian_files) {
    # Extract metacommunity number from filename
    metacom_num <- str_extract(basename(pfile), "(?<=DirichC)\\d+")
    metacom_name <- paste0("C", metacom_num)
    
    # Get corresponding OTU table
    base_name <- "otu_table_treeholes_dada2_renamed_clean_red"
    filename <- sprintf("%s_DirichC%s_rar.csv", base_name, metacom_num)
    otu_file <- file.path(dirichlet_otu_dir, filename)
    
    if (file.exists(otu_file)) {
      # Read files
      parmedian_data <- readLines(pfile)
      data_lines <- parmedian_data[!grepl("^#", parmedian_data)]
      otu_df <- read.csv(otu_file)
      
      # Get sample IDs from OTU table (excluding first column which contains OTUs)
      sample_ids <- colnames(otu_df)[-1]
      
      # Extract immigration rates (skip first line as it's metacommunity rate)
      rates <- sapply(data_lines[-1], function(x) {
        as.numeric(strsplit(trimws(x), "\\s+")[[1]][2])
      })
      
      # Validate we have the right number of rates
      if (length(rates) < length(sample_ids)) {
        warning(sprintf("Not enough immigration rates for all samples in %s", basename(pfile)))
        next
      }
      if (length(rates) > length(sample_ids)) {
        rates <- rates[1:length(sample_ids)]
      }
      
      # Create dataframe with sample IDs as row names
      rates_df <- data.frame(
        metacommunity = metacom_name,
        immigration_rate = rates,
        row.names = sample_ids,
        stringsAsFactors = FALSE
      )
      
      # Add to results, preserving row names
      if (nrow(all_rates) == 0) {
        all_rates <- rates_df
      } else {
        all_rates <- rbind(all_rates, rates_df)
      }
    }
  }
  
  return(all_rates)
}

get_metacommunity_rates <- function(parmedian_dir) {
  # Get all ParMedian files
  parmedian_files <- list.files(
    path = parmedian_dir,
    pattern = "ParMedians.txt$",
    full.names = TRUE
  )
  
  # Initialize empty dataframe for results
  all_rates <- data.frame()
  
  # Process each ParMedian file
  for (pfile in parmedian_files) {
    # Extract metacommunity number from filename
    metacom_num <- str_extract(basename(pfile), "(?<=DirichC)\\d+")
    metacom_name <- paste0("C", metacom_num)
    
    # Read file and get first data line (metacommunity rate)
    parmedian_data <- readLines(pfile)
    data_lines <- parmedian_data[!grepl("^#", parmedian_data)]
    
    if (length(data_lines) > 0) {
      rate <- as.numeric(strsplit(trimws(data_lines[1]), "\\s+")[[1]][2])
      
      # Add to results
      rates_df <- data.frame(
        metacommunity = metacom_name,
        speciation_rate = rate,
        stringsAsFactors = FALSE
      )
      
      all_rates <- rbind(all_rates, rates_df)
    }
  }
  
  return(all_rates)
}
