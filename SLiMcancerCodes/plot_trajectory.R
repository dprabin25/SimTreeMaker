#!/usr/bin/env Rscript
# plot_trajectory.R — draws a treatment run's cells-per-generation trajectory.
# Called by simtreemaker.py for `nomigration treatment` runs.
#
# Usage: Rscript plot_trajectory.R <trajectory.csv> <out.png> <dose generation>
#   top:    total, sensitive and resistant cells per generation
#   bottom: resistant cells as a percentage of the tumor (the outcome after
#           the dose varies by run: resistance can fade or take over)
#   a dashed orange line marks the dose in both panels

args <- commandArgs(trailingOnly = TRUE)
d    <- read.csv(args[1])
out  <- args[2]
dose <- as.numeric(args[3])

blue <- "#0072B2"; orange <- "#D55E00"; grey <- "grey35"
d$pct <- ifelse(d$total > 0, 100 * d$resistant / d$total, NA)

png(out, width = 9, height = 7, units = "in", res = 200)
par(mfrow = c(2, 1), mar = c(4.2, 5.2, 2.4, 1), oma = c(1.2, 0, 2, 0), las = 1)

plot(d$generation, d$total, type = "n", xlab = "Generation", ylab = "",
     ylim = c(0, max(d$total) * 1.05), bty = "l", cex.axis = 0.9)
mtext("Cells", side = 2, line = 4, las = 0)
abline(v = dose, col = orange, lty = 2, lwd = 1.6)
lines(d$generation, d$total, col = grey, lwd = 2, lty = 3)
lines(d$generation, d$sensitive, col = blue, lwd = 2)
lines(d$generation, d$resistant, col = orange, lwd = 2)
legend("topleft", bty = "n", cex = 0.85,
       legend = c("Total", "Sensitive", "Resistant (m3)", "Drug dose"),
       col = c(grey, blue, orange, orange), lty = c(3, 1, 1, 2), lwd = c(2, 2, 2, 1.6))
title("Tumor size: drop at the dose and regrowth", font.main = 1, cex.main = 1)

plot(d$generation, d$pct, type = "n", xlab = "Generation", ylab = "",
     ylim = c(0, max(c(d$pct, 1), na.rm = TRUE) * 1.1), bty = "l", cex.axis = 0.9)
mtext("Resistant cells (%)", side = 2, line = 4, las = 0)
abline(v = dose, col = orange, lty = 2, lwd = 1.6)
polygon(c(d$generation, rev(d$generation)), c(ifelse(is.na(d$pct), 0, d$pct), rep(0, nrow(d))),
        col = adjustcolor(orange, 0.15), border = NA)
lines(d$generation, d$pct, col = orange, lwd = 2)
title("Resistant share of the tumor", font.main = 1, cex.main = 1)

mtext(sprintf("Single drug dose at generation %d (no migration)", dose), outer = TRUE, font = 2, cex = 1.1)
mtext("Counts after each generation's survival step; dashed orange line = dose",
      side = 1, outer = TRUE, cex = 0.75, col = "grey30")
dev.off()
