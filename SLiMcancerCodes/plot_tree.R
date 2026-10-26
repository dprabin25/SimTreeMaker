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

# Population colours (Okabe-Ito), parsed from the "|popN|" field of each tip label.
tip_pop  <- as.integer(sub(".*\\|pop(\\d+)\\|.*", "\\1", tree$tip.label))
pops     <- sort(unique(tip_pop[!is.na(tip_pop)]))
pal      <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
pop_col  <- setNames(pal[(seq_along(pops) - 1) %% length(pal) + 1], pops)
tip_col  <- ifelse(is.na(tip_pop), "grey40", pop_col[as.character(tip_pop)])
# Colour each edge by its descendant tips' population when they all agree, else grey.
edge_col <- local({
  desc <- prop.part(tree)                         # tip sets for each internal node
  node_col <- sapply(desc, function(t) { p <- unique(tip_col[t]); if (length(p) == 1) p else "grey35" })
  all_col  <- c(tip_col, node_col)                # indexed by node number
  all_col[tree$edge[, 2]]
})
lwd_edge <- max(1, min(2.5, 300 / n_tips))       # thinner lines for bigger trees

add_pop_legend <- function() {
  if (length(pops) > 1)
    legend("topleft", legend = paste0("p", pops), col = pop_col, lwd = 3,
           bty = "n", cex = 1.1, title = "Population", inset = c(0.01, 0.01))
}

# Depth axis measured back from the most recent tip: 0 = sampling time.
draw_time_axis <- function(tree, ...) {
  depths <- node.depth.edgelength(tree)
  pp     <- get("last_plot.phylo", envir = ape:::.PlotPhyloEnv)
  fit    <- lm(pp$xx ~ depths)
  max_d  <- max(depths)
  ticks  <- pretty(c(0, max_d))
  ticks  <- ticks[ticks <= max_d]
  axis(side = 1, at = predict(fit, data.frame(depths = max_d - ticks)), labels = ticks, ...)
  mtext("Generations before sampling", side = 1, line = 2.6, font = 2, cex = 1.1)
}
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
  # One label row per tip: ~0.12 in per tip at cex 0.7 keeps labels from overlapping.
  cex_lbl <- 0.7
  png(out1, width = 11, height = max(8, n_tips * 0.12 + 2), units = "in", res = 200)
  # yaxs = "i": R's default 4% axis padding is ~2 in of blank space on a plot this tall.
  par(mar = c(4.5, 1, 3, 1), lend = 1, ljoin = 1, yaxs = "i")
  plot(tree,
       direction       = "rightwards",
       cex             = cex_lbl,
       font            = 1,
       underscore      = TRUE,      # keep gen0_ind82 as-is (ape shows "_" as a space by default)
       tip.color       = tip_col,
       edge.color      = edge_col,
       edge.width      = lwd_edge,
       label.offset    = 0.3,
       align.tip.label = FALSE,
       y.lim           = c(0.5, n_tips + 0.5),   # no blank band above/below the tips
       no.margin       = FALSE,
       main            = paste0(tag, "  |  ", n_tips, " cells"),
       cex.main        = 1.4)
  draw_time_axis(tree, cex.axis = 1.1, font = 2)
  add_pop_legend()
  dev.off()
  cat("[R] Saved:", out1, "\n")
}

# ── 2. Horizontal | No Labels ─────────────────────────────────────────────────
if ("horizontal_no_labels" %in% variants) {
  out2 <- file.path(out_dir, paste0(stem, "_horizontal_no_labels.png"))
  # Fixed, figure-sized canvas: without labels the tips don't need their own rows.
  png(out2, width = 9, height = min(12, max(6, n_tips * 0.025 + 4)), units = "in", res = 300)
  par(mar = c(4.5, 1, 3, 1), lend = 1, ljoin = 1)
  plot(tree,
       direction      = "rightwards",
       show.tip.label = FALSE,
       edge.color     = edge_col,
       edge.width     = lwd_edge,
       no.margin      = FALSE,
       main           = paste0(tag, "  |  ", n_tips, " cells"),
       cex.main       = 1.4)
  draw_time_axis(tree, cex.axis = 1.1, font = 2)
  add_pop_legend()
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
