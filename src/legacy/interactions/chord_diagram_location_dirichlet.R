# Install the packages if needed
# install.packages(c("circlize", "dplyr"))

library(circlize)
library(dplyr)
library(tidyr)

meta <- read.table(file = "../../../../6.Real_data/Single_tree/Metadata_singleTree_wDirich-K10.tsv", header = T, row.names = 1, sep = "\t")
meta <- meta[!is.na(meta$Dirichlet),]

meta$Sample_ID = rownames(meta)
meta$Dirichlet = paste0("META_", meta$Dirichlet)
meta$Location = paste0("TH_", meta$Location)

df_connections <- meta %>%
  group_by(Location, Dirichlet) %>%
  summarise(value = n(), .groups = 'drop') %>%
  rename(from = Location, to = Dirichlet)

df_clean <- df_connections %>%
  drop_na(from, to)

circos.clear()

# 1. FORCE A BASE DATA FRAME AND MAKE THE NAMES UNIQUE
df_safe <- as.data.frame(df_clean) # Drop the dplyr tibble format

# 1. Convert the data to a connection matrix
# Base R xtabs sums the 'values' for each 'from'-'to' pair
connection_matrix <- xtabs(value ~ from + to, data = df_safe)

# 2. Compute distances and do hierarchical clustering
# Order the Locations, grouping those with similar connection profiles
hc_loc <- hclust(dist(connection_matrix))
order_loc <- rownames(connection_matrix)[hc_loc$order]

# Order the Dirichlets the same way (t() transposes the matrix)
hc_dir <- hclust(dist(t(connection_matrix)))
order_dir <- colnames(connection_matrix)[hc_dir$order]

# 3. Build the new optimised order vector
optimized_order <- c(order_loc, order_dir)

# 4. Set the gaps
# (Same logic as before, but using the new lengths to be safe)
gaps <- c(rep(2, length(order_loc) - 1), 15, 
              rep(2, length(order_dir) - 1), 15)


df_safe <- df_safe %>%
  mutate(
    # Add an invisible prefix (a space or a dash) in case names overlap.
    # To be 100% safe, use explicit prefixes:
    from = paste0("Loc_", as.character(from)),
    to = paste0("Dir_", as.character(to))
  )

# 2. DEFINE THE EXACT VECTORS
locations <- unique(df_safe$from)
clusters <- unique(df_safe$to)

strict_order <- c(locations, clusters)

# 3. CREATE THE GAPS
# A normal gap (2) between elements of the same group.
# A large gap (15) after the last location and after the last cluster.
gaps <- c(rep(2, length(locations) - 1), 15, 
              rep(2, length(clusters) - 1), 15)

# 4. SET THE COLOURS (palette names updated with the prefixes)
col_locations <- setNames(hcl.colors(length(locations), "Set 2"), locations)
col_clusters <- setNames(hcl.colors(length(clusters), "Pastel 1"), clusters)
my_colors <- c(col_locations, col_clusters)

# 5. DRAW (clean reset)
circos.clear()
par(mar = c(1, 1, 1, 8))
circos.par(gap.after = gaps)

chordDiagram(
  x = df_safe, 
  order = optimized_order,      # <-- HERE THE ORDER WITHOUT CROSSINGS IS USED
  grid.col = my_colors,
  transparency = 0.5,
  link.border = NA,
  annotationTrack = c("grid"), 
  annotationTrackHeight = c(0.05)
)

# (Then run the loops for the names and legends exactly as before)

# 6. ADD THE ORIGINAL NAMES BY HAND (without the prefix)
# This loop draws the clean names around the circle
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  
  # Remove the "Loc_" or "Dir_" prefix so the text reads well
  clean_name <- gsub("Loc_|Dir_", "", si)
  
  circos.text(mean(xlim), mean(ylim) + 0.5, clean_name, 
              sector.index = si, track.index = 1, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}

# 7. CLEAN LEGENDS
# Clean the names in the legends too
legend_names_loc <- gsub("Loc_", "", names(col_locations))
legend_names_dir <- gsub("Dir_", "", names(col_clusters))

legend("topright", inset = c(-0.2, 0), legend = legend_names_loc, 
       fill = col_locations, title = "Location", bty = "n", cex = 0.8, xpd = TRUE)

legend("topright", inset = c(-0.2, 0.3), legend = legend_names_dir, 
       fill = col_clusters, title = "Dirichlet Cluster", bty = "n", cex = 0.8, xpd = TRUE)
# (The legend code used before can go here again)

# Locations legend (top-right corner)
legend(
  x = "topright", 
  inset = c(-0.2, 0),       # Move the legend into the right margin created above
  legend = names(col_locations), 
  fill = col_locations, 
  title = "Location",
  bty = "n",                # Remove the border box for a cleaner look
  cex = 0.8,                # Font size
  xpd = TRUE                # Allow drawing outside the plot area
)

# Dirichlet legend (just below)
legend(
  x = "topright", 
  inset = c(-0.2, 0.3),     # Move the legend to the right and slightly down
  legend = names(col_clusters), 
  fill = col_clusters, 
  title = "Dirichlet Cluster",
  bty = "n", 
  cex = 0.8,
  xpd = TRUE
)

# (Optional) Add an overall title to the plot
title("Sample correspondence: Location vs Dirichlet")







# 1. The original data frame (e.g. Samples and Variables)
df_original <- data.frame(
  Sample = c("M1", "M1", "M2", "M2", "M3", "M3", "M4"),
  Variable = c("V1", "V2", "V2", "V3", "V1", "V4", "V4")
)

# 2. The clustering result for the second column ("Variable")
df_clusters <- data.frame(
  Variable = c("V1", "V2", "V3", "V4"),
  Cluster = c("Group_A", "Group_A", "Group_B", "Group_B")
)


# Join the cluster information to the original data frame
df_complete <- df_original %>%
  left_join(df_clusters, by = "Variable")

# Compute the "weight" (frequency) of each Sample <-> Cluster link
chord_data <- df_complete %>%
  group_by(Sample, Cluster) %>%
  summarise(Value = n(), .groups = "drop")

# The result is an edge list, which circlize reads directly:
print(chord_data)


# Clear any previous circlize settings
circos.clear()

# Draw the chord diagram
chordDiagram(
  x = chord_data,
  transparency = 0.5,       # Chord transparency
  annotationTrack = c("name", "grid"), # Show the names and the blocks
  annotationTrackHeight = c(0.05, 0.05) # Track height
)

# Add a title (optional)
title("Relationship between samples and clusters")





























df_safe <- df_safe %>% filter(value > 0)

# 2. MATRIX AND CLUSTERING (same as before)
connection_matrix <- xtabs(value ~ from + to, data = df_safe)

hc_loc <- hclust(dist(connection_matrix))
hc_dir <- hclust(dist(t(connection_matrix)))

order_loc <- rownames(connection_matrix)[hc_loc$order]
order_dir <- colnames(connection_matrix)[hc_dir$order]

# The exact order vector
optimized_order <- c(order_loc, order_dir)

# ---------------------------------------------------------
# 3. THE KEY TRICK: A NAMED VECTOR OF GAPS
# ---------------------------------------------------------
# Create a vector of "2", one per sector, and GIVE IT THEIR NAMES
gaps <- rep(2, length(optimized_order))
names(gaps) <- optimized_order

# Now assign the large gap (15) EXACTLY
# to the last Location and the last Dirichlet, whatever their names.
gaps[order_loc[length(order_loc)]] <- 15
gaps[order_dir[length(order_dir)]] <- 15

# ---------------------------------------------------------
# 4. DRAW (now robust)
# ---------------------------------------------------------
circos.clear()
par(mar = c(1, 1, 1, 8))

# With a named vector, circlize matches the sectors internally and does not fail
circos.par(gap.after = gaps) 

chordDiagram(
  x = df_safe, 
  order = optimized_order,      
  grid.col = my_colors,
  transparency = 0.5,
  link.border = NA,
  annotationTrack = c("grid"), 
  annotationTrackHeight = c(0.05)
)

# ---------------------------------------------------------
# 5. LABELS AND LEGENDS
# ---------------------------------------------------------
# Draw the names around the circle again
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  
  clean_name <- gsub("Loc_|Dir_", "", si)
  
  circos.text(mean(xlim), mean(ylim) + 0.5, clean_name, 
              sector.index = si, track.index = 1, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}
