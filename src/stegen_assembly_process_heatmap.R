rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(here)

# ----------------------------
# Load data
# ----------------------------
# Optional arguments (defaults = committed files):
#   Rscript src/stegen_assembly_process_heatmap.R bnti=<csv> rc=<csv> out=<figure path without extension>
opt <- list(bnti = here("3.2_Neutral_Stegen", "betaNTI_treeholes.csv"),
            rc   = here("3.2_Neutral_Stegen", "RC_bray_treeholes.csv"),
            out  = here("3.2_Neutral_Stegen", "figures", "AssemblyProcessHeatmap_treehole_metacommunities_k10_2026"))
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2 || !kv[1] %in% names(opt)) stop("Unknown argument: ", a)
  opt[[kv[1]]] <- kv[2]
}
bNTI <- read.csv(opt$bnti, row.names = 1)
RC   <- read.csv(opt$rc, row.names = 1)
meta <- read.csv(here("1_samples_properties", "Metadata_singleTree_wDirich-K10.tsv"),
                 sep = "\t")

# Keep only matching samples
meta <- meta %>% filter(ID.Sample %in% colnames(bNTI))

# ----------------------------
# Melt + merge
# ----------------------------
df <- melt(as.matrix(RC), varnames = c("X1","X2"), value.name = "RC") %>%
  left_join(
    melt(as.matrix(bNTI), varnames = c("X1","X2"), value.name = "bNTI"),
    by = c("X1","X2")
  )

# ----------------------------
# Add metadata
# ----------------------------
meta_sub <- meta %>% select(ID.Sample, Dirichlet)

df <- df %>%
  left_join(meta_sub, by = c("X1" = "ID.Sample")) %>%
  left_join(meta_sub, by = c("X2" = "ID.Sample"),
            suffix = c(".x", ".y")) %>%
  drop_na(Dirichlet.x, Dirichlet.y)

# ----------------------------
# Define ecological processes
# ----------------------------
df <- df %>%
  mutate(
    Process = case_when(
      bNTI < -2      ~ "Homogeneous selection",
      bNTI >  2      ~ "Variable selection",
      RC   < -0.95   ~ "Homogenizing dispersal",
      RC   >  0.95   ~ "Dispersal limitation",
      TRUE           ~ "Ecological drift"
    ),
    Process = factor(Process, levels = c(
      "Homogeneous selection",
      "Variable selection",
      "Homogenizing dispersal",
      "Dispersal limitation",
      "Ecological drift"
    ))
  )

# ----------------------------
# Order samples
# ----------------------------
levels_X <- df %>%
  distinct(X1, Dirichlet.x) %>%
  arrange(as.numeric(Dirichlet.x)) %>%
  pull(X1)

df <- df %>%
  mutate(
    X1 = factor(X1, levels = levels_X),
    X2 = factor(X2, levels = levels_X)
  ) %>%
  complete(X1, X2) %>%
  # numeric positions (first sample at the top-left), the coordinates used by the
  # boxes and labels below; mixing them with discrete axes breaks in ggplot2 >= 4
  mutate(x = as.integer(X1), y = length(levels_X) + 1 - as.integer(X2))

# ----------------------------
# Metacommunity boxes
# ----------------------------
meta_counts <- meta %>%
  count(Dirichlet, name = "n") %>%
  arrange(as.numeric(Dirichlet))

meta_n <- meta_counts$n
total  <- sum(meta_n)

coords <- lapply(seq_along(meta_n), function(i) {
  w <- meta_n[i]
  h <- meta_n[i]
  
  x0 <- sum(meta_n[seq_len(i-1)])
  x1 <- x0 + w
  y0 <- total - sum(meta_n[seq_len(i-1)])
  y1 <- y0 - h
  
  data.frame(
    x = c(x0, x1, x1, x0),
    y = c(y0, y0, y1, y1),
    xend = c(x1, x1, x0, x0),
    yend = c(y0, y1, y1, y0)
  )
})

boxes <- bind_rows(coords)

# ----------------------------
# Labels
# ----------------------------
meta_centers <- meta_counts %>%
  mutate(center = cumsum(n) - n / 2)

label_df <- data.frame(
  x = c(meta_centers$center, rep(-10, nrow(meta_centers))),
  y = c(rep(-10, nrow(meta_centers)), total - meta_centers$center),
  label = rep(paste("Meta", meta_centers$Dirichlet), 2)
)

# ----------------------------
# Plot
# ----------------------------
metac_plot <- ggplot(df, aes(x, y, fill = Process)) +
  geom_tile() +
  
  scale_fill_manual(
    values = c(
      "Homogeneous selection"   = "gold",
      "Homogenizing dispersal"  = "#2c7bb6",
      "Variable selection"      = "#d9e5d6",
      "Ecological drift"        = "#d7191c",
      "Dispersal limitation"    = "lightblue"
    ),
    na.value = "grey80",
    name = "Assembly process",
    guide = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  
  
  geom_segment(data = boxes,
               aes(x = x+0.5, y = y+0.5,
                   xend = xend+0.5, yend = yend+0.5),
               inherit.aes = FALSE,
               color = "black", size = 0.5) +
  
  annotate("text",
           x = label_df$x,
           y = label_df$y,
           label = label_df$label,
           angle = c(rep(90, nrow(meta_centers)), rep(0, nrow(meta_centers))),
           size = 3) +
  
  coord_cartesian(xlim = c(0.5, total + 0.5), ylim = c(0.5, total + 0.5),
                  expand = FALSE, clip = "off") +
  
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "top",
    panel.background = element_blank(),
    plot.margin = margin(0.5, 1, 1.5, 1.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

metac_plot

# ----------------------------
# Save output
# ----------------------------
ggsave(
  filename = paste0(opt$out, ".tiff"),
  plot = metac_plot,
  units = "cm",
  width = 12,
  height = 12,
  dpi = 300
)

pdf(paste0(opt$out, ".pdf"),
    width = 12/2.54, height = 12/2.54)  
metac_plot
dev.off()