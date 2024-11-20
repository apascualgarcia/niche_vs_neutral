# Check and install required packages
required_packages <- c("phyloseq", "vegan", "ggplot2", "dplyr", "readr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Set up logging
log_message <- function(msg) {
  message(paste0(Sys.time(), ": ", msg))
}

# Import data (adjust file paths as necessary)
load_phyloseq_object <- function() {
  log_message("Starting to load phyloseq object")
  
  mapfile = "1_samples_properties/Metadata_singleTree_wDirich.tsv"
  otufile = "1_dada2/otu_table_treeholes_dada2_renamed_clean_red.tsv"
  taxfile = "1_dada2/taxonomy.list"
  
  # Import OTU table
  log_message("Importing OTU table")
  otu_table <- read_tsv(otufile, show_col_types = FALSE)
  log_message(paste("OTU table dimensions:", nrow(otu_table), "x", ncol(otu_table)))
  otu_table <- otu_table %>%
    column_to_rownames("#OTU ID") %>%
    as.matrix() %>%
    otu_table(taxa_are_rows = TRUE)
  
  # Import sample data
  log_message("Importing sample data")
  sample_data <- read_tsv(mapfile, show_col_types = FALSE)
  log_message(paste("Sample data dimensions:", nrow(sample_data), "x", ncol(sample_data)))
  sample_data <- sample_data %>%
    column_to_rownames("#SampleID") %>%
    sample_data()
  
  # Import taxonomy data
  log_message("Importing taxonomy data")
  tax_table <- read_tsv(taxfile, show_col_types = FALSE)
  log_message(paste("Taxonomy table dimensions:", nrow(tax_table), "x", ncol(tax_table)))
  tax_table <- tax_table %>%
    separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
    column_to_rownames("#OTU ID") %>%
    as.matrix() %>%
    tax_table()
  
  # Combine into a phyloseq object
  log_message("Combining data into phyloseq object")
  physeq <- phyloseq(otu_table, sample_data, tax_table)
  
  # Remove mitochondria and chloroplast
  log_message("Removing mitochondria and chloroplast")
  physeq <- physeq %>% 
    subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")
  
  log_message("Phyloseq object loaded successfully")
  return(physeq)
}

# Preprocess data
preprocess_data <- function(physeq) {
  log_message("Starting data preprocessing")
  
  # Remove low-abundance taxa and samples with few reads
  log_message("Removing low-abundance taxa and samples with few reads")
  physeq_filtered <- prune_taxa(taxa_sums(physeq) > 10, physeq)
  physeq_filtered <- prune_samples(sample_sums(physeq_filtered) > 1000, physeq_filtered)
  
  # Normalize data (using total sum scaling)
  log_message("Normalizing data using total sum scaling")
  total = median(sample_sums(physeq_filtered))
  standf = function(x, t=total) round(t * (x / sum(x)))
  physeq_normalized <- transform_sample_counts(physeq_filtered, standf)
  
  log_message("Data preprocessing completed")
  return(physeq_normalized)
}

# Perform NMDS
perform_nmds <- function(physeq_normalized) {
  log_message("Starting NMDS analysis")
  
  # Extract OTU table
  otu_table <- as(otu_table(physeq_normalized), "matrix")
  
  # Calculate distance matrix
  log_message("Calculating Bray-Curtis distance matrix")
  dist_matrix <- vegdist(t(otu_table), method = "bray")
  
  # Perform NMDS
  log_message("Performing NMDS")
  nmds_result <- metaMDS(dist_matrix, k = 2)
  
  log_message("NMDS analysis completed")
  return(nmds_result)
}

# Visualize NMDS results
visualize_nmds <- function(nmds_result, physeq_normalized) {
  log_message("Starting NMDS visualization")
  
  # Extract NMDS coordinates
  nmds_coords <- as.data.frame(nmds_result$points)
  
  # Add sample metadata
  nmds_coords$Sample <- rownames(nmds_coords)
  metadata <- as(sample_data(physeq_normalized), "data.frame")
  nmds_coords <- left_join(nmds_coords, metadata, by = c("Sample" = "Sample"))
  
  # Create plot
  log_message("Creating NMDS plot")
  p <- ggplot(nmds_coords, aes(x = MDS1, y = MDS2, color = Location)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    labs(title = "NMDS of Microbial Communities",
         x = "NMDS1",
         y = "NMDS2")
  
  # Save plot
  log_message("Saving NMDS plot as PDF")
  ggsave("nmds_plot.pdf", p, width = 10, height = 8)
  
  log_message("NMDS visualization completed")
  return(p)
}

# Main execution
main <- function() {
  tryCatch({
    physeq <- load_phyloseq_object()
    physeq_normalized <- preprocess_data(physeq)
    nmds_result <- perform_nmds(physeq_normalized)
    plot <- visualize_nmds(nmds_result, physeq_normalized)
    print(plot)
    log_message("Analysis completed successfully")
  }, error = function(e) {
    log_message("An error occurred:")
    log_message(as.character(e))
    if(inherits(e, "readr_error")) {
      log_message("Error reading file. Please check the file paths and formats.")
    }
  })
}

# Run the analysis
# log_message("Starting NMDS analysis script")
# main()
# log_message("NMDS analysis script finished")