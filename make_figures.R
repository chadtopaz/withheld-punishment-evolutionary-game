## make_figures.R
## ------------------------------------------------------------------------------
## Deterministic figures for the revised manuscript (base R only, no packages).
##   fig_regimes.pdf        Figure 2: exploitation regime, (a) attracting heteroclinic
##                          cycle and (b) attracting interior equilibrium
##   fig_accommodation.pdf  Figure 3: accommodation regime, (a) regime diagram,
##                          (b) bifurcation diagram, (c) exclusion inside the band,
##                          (d) three-strategy coexistence
## Usage:  source("make_figures.R")      (runtime: well under a minute)
## Figures are written with cairo_pdf() when available so that the fonts are embedded
## (pdf() is used otherwise).
## ------------------------------------------------------------------------------
## Output location: the folder containing this script when it is source()d, otherwise the
## working directory.  Run from the "npj Revision 1" folder so that the two coincide.
out_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(out_dir) || length(out_dir) != 1 || is.na(out_dir) || !nzchar(out_dir)) out_dir <- getwd()
cat("Writing outputs to:", out_dir, "\n")
## PDF device with embedded fonts
open_pdf <- function(file, width, height) {
  if (isTRUE(capabilities("cairo"))) cairo_pdf(file, width = width, height = height, family = "Helvetica")
  else pdf(file, width = width, height = height)
}
## Everything printed below is also written to figures_log.txt in that folder.  The
## sink and any graphics device this script opens are closed on exit, whether the
## script finishes, stops at a check, or is interrupted; sinks and devices that were
## open before are left alone.
.sink_depth0 <- sink.number(); .dev0 <- dev.cur()
sink(file.path(out_dir, "figures_log.txt"), split = TRUE)
tryCatch({

## ---- replicator dynamics helpers ----------------------------------------------
rhs <- function(x, M) { p <- as.vector(M %*% x); x * (p - sum(x * p)) }
## Fixed-step RK4.  The clip at 1e-300 keeps a share from underflowing to exactly zero
## (in the continuum model a share never reaches zero; on the attracting cycle of
## Figure 2a the public share falls far below 1e-300 late in the run, so the clip is
## what keeps the cycle going numerically) and the renormalization removes rounding
## drift.  The raw RK4 proposal is inspected before both, and the attribute "rk4"
## records the smallest proposed component, the number of steps with a negative
## proposal, the largest drift of the coordinate sum from one, and the number of
## components the clip changed; rk4_check() below stops on a negative proposal or
## on drift above 1e-9 and reports the clip count.
integrate <- function(M, x0, T, dt = 0.01) {
  n_steps <- ceiling(T / dt); out <- matrix(NA_real_, n_steps + 1, 3)
  x <- x0 / sum(x0); out[1, ] <- x
  min_raw <- Inf; n_neg <- 0L; max_drift <- 0; n_clip <- 0L
  for (s in 1:n_steps) {
    k1 <- rhs(x, M); k2 <- rhs(x + 0.5 * dt * k1, M); k3 <- rhs(x + 0.5 * dt * k2, M); k4 <- rhs(x + dt * k3, M)
    x_raw <- x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    min_raw <- min(min_raw, min(x_raw)); if (any(x_raw < 0)) n_neg <- n_neg + 1L
    max_drift <- max(max_drift, abs(sum(x_raw) - 1)); n_clip <- n_clip + sum(x_raw < 1e-300)
    x <- pmax(x_raw, 1e-300); x <- x / sum(x)
    out[s + 1, ] <- x
  }
  attr(out, "rk4") <- list(min_raw = min_raw, n_neg = n_neg, max_drift = max_drift, n_clip = n_clip)
  out
}
rk4_check <- function(label, paths) {
  r <- lapply(paths, attr, "rk4")
  n_neg <- sum(vapply(r, function(z) z$n_neg, integer(1))); max_drift <- max(vapply(r, function(z) z$max_drift, numeric(1)))
  min_raw <- min(vapply(r, function(z) z$min_raw, numeric(1))); n_clip <- sum(vapply(r, function(z) z$n_clip, integer(1)))
  cat(sprintf("%s RK4 step check (raw proposals): %d negative proposals (0 required); max drift of the coordinate sum %.1e (gate 1e-9); smallest proposed component %.1e; components clipped at 1e-300: %d\n",
              label, n_neg, max_drift, min_raw, n_clip))
  if (n_neg > 0 || max_drift >= 1e-9) stop(label, ": RK4 step check failed")
}
interior_eq <- function(M) {
  A <- rbind(M[1, ] - M[2, ], M[2, ] - M[3, ], c(1, 1, 1))
  as.vector(solve(A, c(0, 0, 1)))
}

## ---- simplex drawing (vertices: D top, P bottom-left, A bottom-right) ---------
PD <- c(0.5, sqrt(3) / 2); PP <- c(0, 0); PA <- c(1, 0)
bary2cart <- function(x) x[1] * PD + x[2] * PP + x[3] * PA
## The window is kept tight around the triangle and its labels, and the vertex labels are
## centred above D and below P and A, so that the triangle fills most of each panel when
## the figure is included at the text width.
setup_ax <- function() {
  plot(NA, xlim = c(-0.26, 1.26), ylim = c(-0.34, 1.10), asp = 1, axes = FALSE, xlab = "", ylab = "")
  tri <- rbind(PD, PP, PA, PD); lines(tri[, 1], tri[, 2], col = "grey25", lwd = 1.8)
  text(PD[1], PD[2] + 0.075, expression(bolditalic(D)), cex = 1.6)
  text(PD[1], PD[2] + 0.19, "defenders", cex = 1.1, font = 3)
  text(PP[1], PP[2] - 0.075, expression(bolditalic(P)), cex = 1.6)
  text(PP[1], PP[2] - 0.19, "non-punishing\npublic", cex = 1.1, font = 3, adj = c(0.5, 1))
  text(PA[1], PA[2] - 0.075, expression(bolditalic(A)), cex = 1.6)
  text(PA[1], PA[2] - 0.19, "disruptors", cex = 1.1, font = 3, adj = c(0.5, 1))
}
plot_traj <- function(pts, col, arrows_at = c(0.25, 0.55, 0.8), lwd = 1.6, by = c("index", "arc")) {
  ## Arrowheads are placed at the given fractions either of the time points (by = "index")
  ## or of the arc length of the drawn curve (by = "arc", useful when the trajectory converges
  ## quickly).  Each arrowhead is drawn over the shortest stretch of the curve that is at least
  ## 0.02 long, so that it points along the trajectory; a fraction that falls where the
  ## trajectory has stopped moving is skipped.
  by <- match.arg(by)
  xy <- t(apply(pts, 1, bary2cart))
  lines(xy[, 1], xy[, 2], col = col, lwd = lwd)
  s <- c(0, cumsum(sqrt(rowSums(diff(xy)^2))))
  for (fr in arrows_at) {
    i <- if (by == "index") max(1, floor(fr * (nrow(xy) - 2))) else max(1, findInterval(fr * s[length(s)], s))
    j <- findInterval(s[i] + 0.02, s) + 1
    if (j > nrow(xy) || j <= i) next
    arrows(xy[i, 1], xy[i, 2], xy[j, 1], xy[j, 2], length = 0.085, angle = 25, col = col, lwd = lwd)
  }
}
edge_arrows <- function(cycle = TRUE) {
  ## boundary flow markers: D -> P, P -> A, A -> D for the cyclic case
  segs <- list(c(1, 2), c(2, 3), c(3, 1)); V <- list(PD, PP, PA)
  for (s in segs) {
    p <- V[[s[1]]]; q <- V[[s[2]]]; m <- 0.5 * (p + q); d <- (q - p) / sqrt(sum((q - p)^2))
    a <- m - 0.09 * d; b <- m + 0.09 * d
    arrows(a[1], a[2], b[1], b[2], length = 0.1, angle = 25, col = "grey45", lwd = 2.0)
  }
}
mark_eq <- function(x, col, stable = TRUE, cex = 1.3, halo = FALSE) {
  ## halo = TRUE draws a white disc behind the marker so that it stands out where
  ## converging trajectories crowd the equilibrium (Figure 2b).
  cc <- bary2cart(x)
  if (halo) points(cc[1], cc[2], pch = 21, bg = "white", col = "white", cex = 1.7 * cex)
  points(cc[1], cc[2], pch = 21, bg = if (stable) col else "white", col = col, lwd = 1.5, cex = cex)
}
panel_label <- function(lab) mtext(lab, side = 3, adj = 0, line = -0.2, font = 2, cex = 1.25)

colE <- "#1f4e8c"   # exploitation regime (blue)
colC <- "#7a1f7a"   # accommodation regime (purple)

## ==============================================================================
## Figure 2: exploitation regime
## ==============================================================================
## (a) attracting heteroclinic cycle: a_P a_A a_D = 3 > b_P b_A b_D = 0.5
aD <- 1; bD <- 0.5; aP <- 3; bP <- 1; aA <- 1; bA <- 1
M1a <- matrix(c(0, -aD, bD,  bP, 0, -aP,  -aA, bA, 0), 3, byrow = TRUE)
## (b) attracting interior equilibrium: a_P a_A a_D = 1 < b_P b_A b_D = 3.375
aD2 <- 1; bD2 <- 1.5; aP2 <- 1; bP2 <- 1.5; aA2 <- 1; bA2 <- 1.5
M1b <- matrix(c(0, -aD2, bD2,  bP2, 0, -aP2,  -aA2, bA2, 0), 3, byrow = TRUE)
stopifnot(aP * aA * aD > bP * bA * bD, aP2 * aA2 * aD2 < bP2 * bA2 * bD2)
xs_a <- interior_eq(M1a); xs_b <- interior_eq(M1b)
cat(sprintf("Fig 2a interior equilibrium (D,P,A) = (%.3f, %.3f, %.3f), repelling\n", xs_a[1], xs_a[2], xs_a[3]))
cat(sprintf("Fig 2b interior equilibrium (D,P,A) = (%.3f, %.3f, %.3f), attracting\n", xs_b[1], xs_b[2], xs_b[3]))

## 7.6 x 3.9 in, included at the text width (6.5 in): the panel height matches the
## content's aspect ratio, so the triangles fill the panels instead of floating in
## white space, and the labels print at about 9-13 pt.
## trajectories are integrated and checked before any device is opened
x0 <- xs_a + c(0.004, -0.002, -0.002)
traj_a <- integrate(M1a, x0, T = 700)
starts <- list(c(0.85, 0.10, 0.05), c(0.05, 0.85, 0.10), c(0.10, 0.05, 0.85), c(0.34, 0.33, 0.33) + c(0.3, -0.15, -0.15))
traj_b <- lapply(starts, function(s) integrate(M1b, s, T = 60))
rk4_check("Fig 2a", list(traj_a)); rk4_check("Fig 2b", traj_b)

open_pdf(file.path(out_dir, "fig_regimes.pdf"), width = 7.6, height = 3.9)
par(mfrow = c(1, 2), mar = c(0.3, 0.3, 1.4, 0.3))
## (a)
setup_ax(); edge_arrows()
## arrowheads at fractions 0.05 and 0.065 of the run fall on the interior spiral (radii about
## 0.08 and 0.22 from the equilibrium), so that the outward motion is visible before the
## trajectory reaches the boundary; the later three sit on the boundary circuit.
plot_traj(traj_a, colE, arrows_at = c(0.05, 0.065, 0.2, 0.45, 0.7))
mark_eq(xs_a, colE, stable = FALSE)
panel_label("(a)")
## (b)
setup_ax()
for (tr in traj_b) plot_traj(tr, colE, arrows_at = c(0.08, 0.3))
mark_eq(xs_b, colE, stable = TRUE, cex = 1.7, halo = TRUE)
panel_label("(b)")
dev.off()

## ==============================================================================
## Figure 3: accommodation regime
## ==============================================================================
## Drawn by make_figure3.R in the same folder (same payoff matrices, parameters,
## starting states, RK4 step and horizons as before; grid graphics, base R only).
## It writes fig_accommodation.pdf and figure3_checks.txt next to itself and prints
## the equilibrium values and the RK4 step statistics it checks.
source(file.path(out_dir, "make_figure3.R"))
cat("Wrote fig_regimes.pdf and fig_accommodation.pdf to", out_dir, "\n")

}, finally = {
  if (dev.cur() != .dev0 && dev.cur() > 1) dev.off()
  while (sink.number() > .sink_depth0) sink()
})
