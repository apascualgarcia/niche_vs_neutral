###############################################################################
# Helpers to save ggplot figures and their legends as separate PDFs
#
# save_without_legend() writes the plot without its legend and, in another PDF,
# the legend alone (cropped to its size). Figures that share a legend pass the
# same `legend_file`: the legend is written only once per run; if a later
# figure has a different legend under the same name (e.g. a metacommunity is
# missing), it is written to <figure>_legend.pdf instead and a message says so.
#
# save_metacommunity_legend() writes the legend shared by all the figures that
# colour the metacommunities with meta_palette (Meta 1 ... Meta 10, without M9).
#
# Sourced by the plotting scripts: source(here("src", "legend_utils.R"))
###############################################################################
library(ggplot2)
library(grid)

.legend_registry <- new.env()

# Legend grob of a ggplot (NULL if the plot has no legend)
get_legend_grob <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "right"))
  boxes <- g$grobs[grepl("^guide-box", g$layout$name)]
  boxes <- boxes[!vapply(boxes, inherits, logical(1), "zeroGrob")]
  if (length(boxes) == 0) return(NULL)
  boxes[[1]]
}

# Content of the legend (titles, labels and mapped values of the non-position scales)
legend_signature <- function(p) {
  b <- ggplot_build(p)
  sc <- b$plot$scales$non_position_scales()$scales
  paste(vapply(sc, function(s) {
    if (identical(s$guide, "none")) return("")
    br <- s$get_breaks()
    paste(s$aesthetics[1], s$make_title(s$name), paste(s$get_labels(br), collapse = ","),
          paste(s$map(br), collapse = ","), sep = "|")
  }, character(1)), collapse = "||")
}

write_legend <- function(leg, file) {
  w <- convertWidth(sum(leg$widths), "in", valueOnly = TRUE) + 0.3
  h <- convertHeight(sum(leg$heights), "in", valueOnly = TRUE) + 0.3
  pdf(file, width = w, height = h)
  grid.newpage()
  grid.draw(leg)
  invisible(dev.off())
}

# p: ggplot; file: PDF of the figure; legend_file: PDF of the legend
# (default <figure>_legend.pdf; NULL = do not write the legend, e.g. when it is
# the shared metacommunity legend written by save_metacommunity_legend())
save_without_legend <- function(p, file, width, height,
                                legend_file = sub("\\.pdf$", "_legend.pdf", file)) {
  pdf(file, width = width, height = height)
  print(p + theme(legend.position = "none"))
  invisible(dev.off())

  if (is.null(legend_file)) return(invisible(NULL))
  leg <- get_legend_grob(p)
  if (is.null(leg)) return(invisible(NULL))

  sig <- legend_signature(p)
  key <- normalizePath(legend_file, mustWork = FALSE)
  if (exists(key, envir = .legend_registry, inherits = FALSE)) {
    if (identical(get(key, envir = .legend_registry), sig)) return(invisible(NULL))  # already written
    own <- sub("\\.pdf$", "_legend.pdf", file)
    message("Legend of ", basename(file), " differs from ", basename(legend_file),
            "; written to ", basename(own))
    legend_file <- own
    key <- normalizePath(own, mustWork = FALSE)
  }
  write_legend(leg, legend_file)
  assign(key, sig, envir = .legend_registry)
  invisible(NULL)
}

# Shared legend of the metacommunities; palette = meta_palette of the calling script
# (the same 9 colours in every script, in the order Meta 1 ... Meta 8, Meta 10)
save_metacommunity_legend <- function(file, palette) {
  key <- normalizePath(file, mustWork = FALSE)
  if (exists(key, envir = .legend_registry, inherits = FALSE)) return(invisible(NULL))
  labs_meta <- paste("Meta", c(1:8, 10))
  d <- data.frame(MetaCom = factor(labs_meta, levels = labs_meta), x = 1, y = 1)
  p <- ggplot(d, aes(x, y, fill = MetaCom)) +
    geom_tile() +
    scale_fill_manual(values = palette[seq_along(labs_meta)]) +
    labs(fill = "Metacommunity") +
    theme_bw() +
    theme(legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.key.size = unit(1.5, "lines"))
  write_legend(get_legend_grob(p), file)
  assign(key, "metacommunity", envir = .legend_registry)
  invisible(NULL)
}
