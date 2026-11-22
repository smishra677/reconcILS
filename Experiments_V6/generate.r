#!/usr/bin/env Rscript
# Generate NEW, reproducible pairs of rooted 40-tip species trees.
# Bundled published trees already exist under "all_species tree/"; do not overwrite them.

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Install the R package 'ape': install.packages('ape')")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript generate.r [output_dir]")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
output_dir <- if (length(args)) args[[1]] else file.path(script_dir, "generated")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(420)
n_trees <- 10L
tip_labels <- c(LETTERS, paste0("A", LETTERS[1:14]))
base_heights <- sample(3:50, size = n_trees, replace = TRUE)
base_trees <- lapply(base_heights, function(height) {
  tree <- ape::rtree(length(tip_labels), rooted = TRUE, tip.label = tip_labels)
  tree <- ape::compute.brlen(tree, method = "Grafen")
  tree$edge.length <- tree$edge.length * height
  stopifnot(ape::is.ultrametric(tree))
  tree
})
class(base_trees) <- "multiPhylo"

# The bundled SimPhy counterparts use exactly 10x the base branch lengths.
simphy_trees <- lapply(base_trees, function(tree) {
  tree$edge.length <- tree$edge.length * 10
  tree
})
class(simphy_trees) <- "multiPhylo"

ape::write.tree(base_trees,
  file.path(output_dir, "species_40_scaled_10.tre"), digits = 15)
ape::write.tree(simphy_trees,
  file.path(output_dir, "species_40_scaled_10_simphy.tre"), digits = 15)
writeLines(as.character(base_heights), file.path(output_dir, "scale_factors_simphy.txt"))
writeLines(as.character(base_heights * 10), file.path(output_dir, "simphy_tree_heights.txt"))
cat(sprintf("Generated %d paired 40-tip trees (base and x10 SimPhy) in %s\n",
  n_trees, normalizePath(output_dir)))
