#!/usr/bin/env Rscript
# plot_tree.R — called by simtreemaker.py to render Newick trees as PNGs
# Usage: Rscript plot_tree.R <nwk_path> <out_dir> <stem> <model_name> <model_type> <mode> [variants]
#   variants (optional): comma-separated subset of horizontal_labels,
#   horizontal_no_labels, vertical_labels, vertical_no_labels. Omit to
#   render all four (the default).
# Requires: ape

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript plot_tree.R <nwk_path> <out_dir> <stem> <model_name> <model_type> <mode> [variants]\n")
  quit(status = 1)
}

nwk_path   <- args[1]
out_dir    <- args[2]
stem       <- args[3]
model_name <- args[4]
model_type <- args[5]
mode       <- if (length(args) >= 6) args[6] else "Snapshot"
variants   <- if (length(args) >= 7 && nzchar(args[7])) strsplit(args[7], ",")[[1]] else
              c("horizontal_labels", "horizontal_no_labels", "vertical_labels", "vertical_no_labels")

if (!requireNamespace("ape", quietly = TRUE)) {
  cat("[ERROR] R package 'ape' is not installed. Run: install.packages('ape')\n")
  quit(status = 1)
}
library(ape)

tree <- read.tree(nwk_path)
if (is.null(tree)) {
  cat("[ERROR] Could not read tree from", nwk_path, "\n")
  quit(status = 1)
}

n_tips  <- length(tree$tip.label)
cex_tip <- max(0.9, min(1.8, 80 / n_tips))   # readable down to ~100 tips
w_lbl   <- max(12, n_tips * 0.22)             # wide enough for one label per tip
w_clean <- max(8,  n_tips * 0.12)             # no-label plots can be narrower
h_lbl   <- max(8,  n_tips * 0.18)             # horizontal label height
h_clean <- max(5,  n_tips * 0.12)             # horizontal no-label height
tag     <- paste0(model_name, "  [", model_type, "]  |  ", mode)

cat(sprintf("[R] Tree: %d tips | output -> %s\n", n_tips, out_dir))

# axisPhylo() mislabels "downwards"-direction trees (root/tip values come out swapped), so
# every panel below draws its depth axis manually instead: fit each node's true root-to-node
# distance against its plotted coordinate, then place ticks from that fit. Works for both
# horizontal (side=1, x-coordinate) and vertical (side=2, y-coordinate) layouts.
draw_depth_axis <- function(tree, side, ...) {
  depths_all <- node.depth.edgelength(tree)
  pp         <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
  coord      <- if (side == 1) pp$xx else pp$yy
  fit        <- lm(coord ~ depths_all)
  max_depth  <- max(depths_all)
  ticks      <- pretty(c(0, max_depth))
  ticks      <- ticks[ticks >= 0 & ticks <= max_depth]
  at_ticks   <- predict(fit, newdata = data.frame(depths_all = ticks))
  axis(side = side, at = at_ticks, labels = ticks, las = if (side == 2) 1 else 0, ...)
}

# ── 1. Horizontal | With Labels ───────────────────────────────────────────────
if ("horizontal_labels" %in% variants) {
  out1 <- file.path(out_dir, paste0(stem, "_horizontal_labels.png"))
  png(out1, width = 14, height = h_lbl, units = "in", res = 300)
  par(mar = c(4, 1, 2, 10), lend = 1, ljoin = 1)
  plot(tree,
       direction  = "rightwards",
       cex        = cex_tip,
       edge.width = 3.5,
       font       = 2,
       no.margin  = FALSE,
       main       = paste0(tag, "  —  Horizontal | With Labels"),
       cex.main   = 1.5)
  draw_depth_axis(tree, side = 1, cex.axis = 1.1, font = 2)
  mtext("Generations", side = 1, line = 2.6, font = 2, cex = 1.1)
  dev.off()
  cat("[R] Saved:", out1, "\n")
}

# ── 2. Horizontal | No Labels ─────────────────────────────────────────────────
if ("horizontal_no_labels" %in% variants) {
  out2 <- file.path(out_dir, paste0(stem, "_horizontal_no_labels.png"))
  png(out2, width = 10, height = h_clean, units = "in", res = 300)
  par(mar = c(4, 1, 2, 1), lend = 1, ljoin = 1)
  plot(tree,
       direction      = "rightwards",
       show.tip.label = FALSE,
       edge.width     = 3.5,
       no.margin      = FALSE,
       main           = paste0(tag, "  —  Horizontal | No Labels"),
       cex.main       = 1.5)
  draw_depth_axis(tree, side = 1, cex.axis = 1.1, font = 2)
  mtext("Generations", side = 1, line = 2.6, font = 2, cex = 1.1)
  dev.off()
  cat("[R] Saved:", out2, "\n")
}

# ── 3. Vertical | With Labels ─────────────────────────────────────────────────
if ("vertical_labels" %in% variants) {
  out3 <- file.path(out_dir, paste0(stem, "_vertical_labels.png"))
  png(out3, width = w_lbl, height = 14, units = "in", res = 300)
  par(mar = c(10, 5, 2, 1), lend = 1, ljoin = 1)
  plot(tree,
       direction      = "downwards",
       show.tip.label = FALSE,
       edge.width     = 3.5,
       no.margin      = FALSE,
       main           = paste0(tag, "  —  Vertical | With Labels"),
       cex.main       = 1.5)
  draw_depth_axis(tree, side = 2, cex.axis = 1.1, font = 2)
  mtext("Generations", side = 2, line = 3.2, font = 2, cex = 1.1)
  pp <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
  n  <- Ntip(tree)
  par(xpd = TRUE)
  text(x      = pp$xx[1:n],
       y      = pp$yy[1:n],
       labels = tree$tip.label,
       srt    = 90,
       adj    = c(1.05, 0.5),
       cex    = cex_tip,
       font   = 2)
  dev.off()
  cat("[R] Saved:", out3, "\n")
}

# ── 4. Vertical | No Labels ───────────────────────────────────────────────────
if ("vertical_no_labels" %in% variants) {
  out4 <- file.path(out_dir, paste0(stem, "_vertical_no_labels.png"))
  png(out4, width = w_clean, height = 10, units = "in", res = 300)
  par(mar = c(1, 5, 2, 1), lend = 1, ljoin = 1)
  plot(tree,
       direction      = "downwards",
       show.tip.label = FALSE,
       edge.width     = 3.5,
       no.margin      = FALSE,
       main           = paste0(tag, "  —  Vertical | No Labels"),
       cex.main       = 1.5)
  draw_depth_axis(tree, side = 2, cex.axis = 1.1, font = 2)
  mtext("Generations", side = 2, line = 3.2, font = 2, cex = 1.1)
  dev.off()
  cat("[R] Saved:", out4, "\n")
}

cat("[R] All plots complete.\n")
