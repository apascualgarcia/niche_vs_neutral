#!/usr/bin/env Rscript

# Function to count ASVs in LocalModel section of relAbMetac file
count_local_model_asvs <- function(file_path) {
  # Read the file
  lines <- readLines(file_path)
  
  # Find the #LocalModel section
  local_model_start <- which(grepl("^#LocalModel", lines))
  if (length(local_model_start) == 0) {
    stop("No #LocalModel section found")
  }
  
  # Get lines after #LocalModel until next section or end
  next_section <- which(grepl("^#", lines[(local_model_start + 1):length(lines)]))
  if (length(next_section) > 0) {
    end_idx <- local_model_start + next_section[1] - 1
  } else {
    end_idx <- length(lines)
  }
  
  # Get the data lines (skip the #LocalModel line itself)
  local_model_lines <- lines[(local_model_start + 1):(end_idx - 1)]
  # Remove empty lines and count only non-empty lines
  local_model_lines <- local_model_lines[local_model_lines != ""]
  
  # Count lines that have two columns (index and value)
  valid_lines <- sapply(strsplit(local_model_lines, "\\s+"), length) == 2
  return(sum(valid_lines))
}

# Function to count ASVs in rar.csv file
count_rar_asvs <- function(file_path) {
  # Read the CSV file
  data <- read.csv(file_path, header = TRUE, row.names = 1)
  return(nrow(data))
}

# Print a separator line
print_separator <- function() {
  cat(sprintf("%s\n", paste(rep("-", 60), collapse="")))
}

# Print header
cat("\nComparing ASV counts between LocalModel and rar.csv files\n")
print_separator()

# Loop through metacommunities C1 to C10
for (i in 1:10) {
  # Construct file paths
  relabmetac_file <- sprintf("3.1_Neutral_NMGS/Results_K10/otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_relAbMetac.txt", i)
  rar_file <- sprintf("2_Dirichlet/otu_tables_k10/otu_table_treeholes_dada2_renamed_clean_red_DirichC%d_rar.csv", i)
  
  # Count ASVs in both files
  tryCatch({
    local_model_count <- count_local_model_asvs(relabmetac_file)
    rar_count <- count_rar_asvs(rar_file)
    
    # Print results
    cat(sprintf("Metacommunity C%d:\n", i))
    cat(sprintf("  Number of ASVs in LocalModel section: %d\n", local_model_count))
    cat(sprintf("  Number of ASVs in rar.csv file: %d\n", rar_count))
    cat(sprintf("  Difference (LocalModel - rar.csv): %d\n", local_model_count - rar_count))
    print_separator()
  }, error = function(e) {
    cat(sprintf("Error processing metacommunity C%d: %s\n", i, e$message))
    print_separator()
  })
}
