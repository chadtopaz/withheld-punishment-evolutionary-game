# Figure 3: the accommodation regime.
# Uses only packages supplied with R; no add-on packages are needed.
library(grid)

# Run Rscript make_figure3.R, or source this file. Outputs (fig_accommodation.pdf and
# figure3_checks.txt) go beside this file.
script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_path)) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) script_path <- sub("^--file=", "", arg[1])
}
out_dir <- if (is.null(script_path)) getwd() else dirname(normalizePath(script_path))

make_figure3 <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Same model, parameters, starting states, RK4 step and integration horizons
  # as the supplied make_figures.R. Only graphical presentation is changed.
  aD <- bD <- bP <- aA <- cA <- 1
  cP_c <- 0.7
  cP_d <- 0.3
  cP_star <- bD * cA / (cA + aD)
  Hfun <- function(cP) bD * cA - cP * (cA + aD)
  M2 <- function(cP) matrix(c(0, -aD, bD, bP, 0, cP, -aA, cA, 0),
                            nrow = 3, byrow = TRUE)
  xdag <- function(cP) c(0, cP / (cP + cA), cA / (cP + cA))
  xint <- function(cP) {
    N <- c(Hfun(cP), bD * (aA + bP) - aA * cP,
           aD * (aA + bP) + bP * cA)
    N / sum(N)
  }
  rhs <- function(x, M) {
    payoff <- as.vector(M %*% x)
    x * (payoff - sum(x * payoff))
  }
  # Fixed-step RK4.  The clip at 1e-300 protects against underflow and the
  # renormalization removes rounding drift; the raw RK4 proposal is inspected
  # before both, and the attributes record the smallest proposed component, the
  # number of steps with a negative proposal, the largest drift of the coordinate
  # sum from one, and the number of components the clip changed.
  integrate_rk4 <- function(M, x0, horizon, dt = 0.01) {
    n <- ceiling(horizon / dt)
    ans <- matrix(NA_real_, n + 1L, 3L)
    x <- x0 / sum(x0)
    ans[1L, ] <- x
    min_raw <- Inf; n_neg <- 0L; max_drift <- 0; n_clip <- 0L
    for (s in seq_len(n)) {
      k1 <- rhs(x, M)
      k2 <- rhs(x + 0.5 * dt * k1, M)
      k3 <- rhs(x + 0.5 * dt * k2, M)
      k4 <- rhs(x + dt * k3, M)
      x_raw <- x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
      min_raw <- min(min_raw, min(x_raw))
      if (any(x_raw < 0)) n_neg <- n_neg + 1L
      max_drift <- max(max_drift, abs(sum(x_raw) - 1))
      n_clip <- n_clip + sum(x_raw < 1e-300)
      x <- pmax(x_raw, 1e-300)
      x <- x / sum(x)
      ans[s + 1L, ] <- x
    }
    attr(ans, "rk4") <- list(min_raw = min_raw, n_neg = n_neg, max_drift = max_drift, n_clip = n_clip)
    ans
  }
  starts <- list(c(0.6, 0.25, 0.15), c(0.5, 0.2, 0.3),
                 c(0.25, 0.55, 0.20), c(0.3, 0.25, 0.45),
                 c(0.45, 0.45, 0.10), c(0.2, 0.3, 0.5))
  paths_c <- lapply(starts, function(s) integrate_rk4(M2(cP_c), s, 60))
  paths_d <- lapply(starts, function(s) integrate_rk4(M2(cP_d), s, 120))

  # Check analytical equilibria and simplex preservation before drawing.
  equal_payoff <- function(M) as.vector(solve(
    rbind(M[1, ] - M[2, ], M[2, ] - M[3, ], c(1, 1, 1)), c(0, 0, 1)))
  stopifnot(Hfun(cP_c) < 0, Hfun(cP_d) > 0,
            abs(cP_star - 0.5) < 1e-12,
            max(abs(equal_payoff(M2(cP_d)) - xint(cP_d))) < 1e-12,
            max(abs(rhs(xdag(cP_c), M2(cP_c)))) < 1e-12,
            max(abs(rhs(xint(cP_d), M2(cP_d)))) < 1e-12)
  # RK4 step check on the raw proposals (before the clip and renormalization):
  # no negative proposal, and drift of the coordinate sum below 1e-9.  The clip
  # count is reported; at these horizons no component approaches 1e-300.
  rk4 <- lapply(c(paths_c, paths_d), attr, "rk4")
  n_neg <- sum(vapply(rk4, function(r) r$n_neg, integer(1)))
  max_drift <- max(vapply(rk4, function(r) r$max_drift, numeric(1)))
  min_raw <- min(vapply(rk4, function(r) r$min_raw, numeric(1)))
  n_clip <- sum(vapply(rk4, function(r) r$n_clip, integer(1)))
  if (n_neg > 0 || max_drift >= 1e-9)
    stop(sprintf("RK4 step check failed: %d negative proposals, max drift %.1e", n_neg, max_drift))

  # Canvas coordinates are explicit: panel headings and plot boxes align.
  W <- 7.2
  H <- 7.0
  ink <- "#20232A"
  muted <- "#5A6068"
  purple <- "#7B287F"
  purple_fill <- "#EFE0EF"
  blue_fill <- "#E1EBF5"
  rule <- "#737780"
  FONT <- "sans"
  txt <- function(label, x, y, size = 12.5, col = ink, face = "plain",
                  just = "centre", rot = 0) {
    grid.text(label, x = x, y = y, default.units = "npc", just = just,
              rot = rot, gp = gpar(fontfamily = FONT, fontsize = size,
                                  col = col, fontface = face, lineheight = 1.08))
  }
  ln <- function(x, y, col = ink, width = 0.9, lty = 1, arrow = NULL) {
    grid.lines(x, y, default.units = "npc", arrow = arrow,
               gp = gpar(col = col, lwd = width, lty = lty,
                         linejoin = "round", lineend = "round"))
  }
  dot <- function(x, y, fill = purple, size = 5.8, shape = 21) {
    grid.points(x, y, default.units = "npc", pch = shape, size = unit(size, "pt"),
                gp = gpar(col = purple, fill = fill, lwd = 1.35))
  }
  # Panel letter only: the journal requires panel titles to appear in the figure
  # legend rather than in the image, so the title argument is not drawn.
  heading <- function(letter, title, col, y) {
    left <- if (col == 1) 0.055 else 0.555
    txt(paste0("(", letter, ")"), left, y, size = 14, face = "bold", just = "left")
  }
  # Both upper data rectangles have identical physical dimensions.
  A <- list(x = c(0.115, 0.475), y = c(0.615, 0.895),
            xr = c(0, 3), yr = c(0, 1.6))
  B <- list(x = c(0.615, 0.975), y = A$y,
            xr = c(0, 1.2), yr = c(-0.022, 0.19))
  mapx <- function(p, z) p$x[1] + diff(p$x) * (z - p$xr[1]) / diff(p$xr)
  mapy <- function(p, z) p$y[1] + diff(p$y) * (z - p$yr[1]) / diff(p$yr)
  ptxt <- function(p, label, x, y, ...) txt(label, mapx(p, x), mapy(p, y), ...)
  pline <- function(p, x, y, ...) ln(mapx(p, x), mapy(p, y), ...)
  axes <- function(p, xticks, xlabels, yticks, ylabels, xlab, ylab) {
    grid.rect(x = mean(p$x), y = mean(p$y), width = diff(p$x), height = diff(p$y),
              gp = gpar(fill = NA, col = ink, lwd = 0.85))
    for (k in seq_along(xticks)) {
      xx <- mapx(p, xticks[k])
      ln(c(xx, xx), c(p$y[1], p$y[1] - 0.007), width = 0.85)
      txt(xlabels[k], xx, p$y[1] - 0.027, size = 12)
    }
    for (k in seq_along(yticks)) {
      yy <- mapy(p, yticks[k])
      ln(c(p$x[1] - 0.007, p$x[1]), c(yy, yy), width = 0.85)
      txt(ylabels[k], p$x[1] - 0.015, yy, size = 12, just = "right")
    }
    txt(xlab, mean(p$x), p$y[1] - 0.077, size = 14)
    txt(ylab, p$x[1] - 0.083, mean(p$y), size = 14, rot = 90)
  }
  triangle_xy <- function(z, centre) {
    # Physical aspect ratio is one; the triangle is exactly equilateral.
    tri_width <- 0.345
    bottom <- 0.083
    cbind(centre - tri_width / 2 + tri_width * (0.5 * z[, 1] + z[, 3]),
          bottom + tri_width * W / H * sqrt(3) / 2 * z[, 1])
  }
  phase_panel <- function(paths, cP, centre, coexist = FALSE) {
    v <- triangle_xy(diag(3), centre)
    ln(v[c(1, 2, 3, 1), 1], v[c(1, 2, 3, 1), 2], width = 1.15)
    txt("D", v[1, 1], v[1, 2] + 0.024, size = 14.5, face = "bold")
    txt("P", v[2, 1] - 0.013, v[2, 2] - 0.025, size = 14.5, face = "bold")
    txt("A", v[3, 1] + 0.013, v[3, 2] - 0.025, size = 14.5, face = "bold")
    for (j in seq_along(paths)) {
      xy <- triangle_xy(paths[[j]], centre)
      # Draw only geometrically distinct points, preserving the numerical path.
      arc <- c(0, cumsum(sqrt(rowSums(diff(sweep(xy, 2, c(W, H), "*"))^2))))
      keep <- c(TRUE, diff(floor(arc / 0.0012)) > 0)
      keep[length(keep)] <- TRUE
      ln(xy[keep, 1], xy[keep, 2], col = purple, width = 1.0)
      # One clear arrow per trajectory, before the curves bunch at equilibrium.
      fr <- c(0.40, 0.48, 0.52, 0.43, 0.55, 0.61)[j]
      i <- max(1L, findInterval(fr * tail(arc, 1), arc))
      end <- min(nrow(xy), findInterval(arc[i] + 0.07, arc) + 1L)
      if (end > i) ln(xy[c(i, end), 1], xy[c(i, end), 2], col = purple,
                     width = 1.0,
                     arrow = arrow(length = unit(3.7, "pt"), angle = 24,
                                   type = "open"))
    }
    eq <- triangle_xy(matrix(if (coexist) xint(cP) else xdag(cP), nrow = 1), centre)
    # White halo keeps the attracting point visible among converging curves.
    grid.points(eq[1], eq[2], default.units = "npc", pch = 21, size = unit(9.5, "pt"),
                gp = gpar(col = "white", fill = "white"))
    dot(eq[1], eq[2], size = 7)
    if (coexist) {
      saddle <- triangle_xy(matrix(xdag(cP), nrow = 1), centre)
      dot(saddle[1], saddle[2], fill = "white", size = 6.7)
    }
  }
  draw_figure <- function() {
    grid.newpage()
    grid.rect(gp = gpar(fill = "white", col = NA))
    heading("a", "Exclusion threshold", 1, 0.958)
    heading("b", "Equilibrium defender share", 2, 0.958)
    heading("c", "Defender exclusion", 1, 0.492)
    heading("d", "Three-strategy coexistence", 2, 0.492)

    # (a) Exact threshold and the stronger sufficient condition.
    u <- seq(0, 3, length.out = 600)
    threshold <- 1 / (1 + u)
    grid.polygon(mapx(A, c(u, rev(u))), mapy(A, c(threshold, rep(1.6, length(u)))),
                 gp = gpar(fill = purple_fill, col = NA))
    grid.polygon(mapx(A, c(u, rev(u))), mapy(A, c(threshold, rep(0, length(u)))),
                 gp = gpar(fill = blue_fill, col = NA))
    pline(A, u, threshold, col = purple, width = 1.8)
    pline(A, c(0, 3), c(1, 1), col = rule, width = 0.9, lty = "dashed")
    ptxt(A, "Defenders excluded", 1.63, 1.34, size = 12.8, face = "bold")
    ptxt(A, "All three coexist", 1.83, 0.115, size = 12.8, face = "bold")
    ptxt(A, expression(c[P] == b[D]), 0.16, 1.10, size = 12, col = muted, just = "left")
    ptxt(A, expression(H == 0), 2.23, 0.43, size = 12, col = purple)
    dot(mapx(A, 1), mapy(A, cP_c), shape = 17, size = 7)
    dot(mapx(A, 1), mapy(A, cP_d), shape = 16, size = 6)
    ptxt(A, "(c)", 0.87, cP_c, size = 12, just = "right")
    ptxt(A, "(d)", 1.13, cP_d, size = 12, just = "left")
    axes(A, seq(0, 3, .5), sprintf("%.1f", seq(0, 3, .5)),
         c(0, .5, 1, 1.5), c("0.0", "0.5", "1.0", "1.5"),
         expression(a[D] / c[A]), expression(c[P] / b[D]))

    # (b) Defender share. Stable branches solid; saddle branch dotted.
    cp <- seq(0.001, cP_star, length.out = 600)
    xd <- vapply(cp, function(z) xint(z)[1], numeric(1))
    pline(B, c(1, 1), B$yr, col = rule, width = 0.9, lty = "dashed")
    pline(B, cp, xd, col = purple, width = 1.9)
    pline(B, c(cP_star, 1.2), c(0, 0), col = purple, width = 1.9)
    pline(B, c(0, cP_star), c(0, 0), col = purple, width = 1.8, lty = "dotted")
    dot(mapx(B, cP_star), mapy(B, 0), fill = "white", size = 6.8)
    ptxt(B, "Three-strategy\ncoexistence", 0.67, 0.127, size = 12.4, col = purple)
    pline(B, c(.44, .31), c(.110, xint(.31)[1] + .003),
          col = purple, width = 0.85,
          arrow = arrow(length = unit(3.2, "pt"), angle = 24))
    ptxt(B, "Exclusion", .76, .024, size = 12.4, col = purple)
    ptxt(B, "Saddle", .23, -.013, size = 11.8, col = purple)
    ptxt(B, expression(c[P]^"*" == 0.5), .54, .061, size = 12, just = "left")
    pline(B, c(.5, .5), c(.049, .006), col = muted, width = .7)
    ptxt(B, expression(c[P] == b[D]), .97, .176, size = 11.8,
          col = muted, just = "right")
    axes(B, seq(0, 1.2, .3), sprintf("%.1f", seq(0, 1.2, .3)),
         c(0, .05, .10, .15), c("0.00", "0.05", "0.10", "0.15"),
         expression(c[P]), expression(x[D]))

    # (c,d) Much larger, aligned simplex portraits.
    phase_panel(paths_c, cP_c, centre = 0.285)
    phase_panel(paths_d, cP_d, centre = 0.785, coexist = TRUE)
    # A single shared key replaces six repeated labels around the triangles.
    txt("D  defenders", .16, .023, size = 12.4)
    txt("P  non-punishing public", .50, .023, size = 12.4)
    txt("A  disruptors", .85, .023, size = 12.4)
  }
  # The PDF device is closed on exit even if drawing fails; a device that was
  # open before this function ran is left alone.
  pdf_path <- file.path(out_dir, "fig_accommodation.pdf")
  dev_before <- dev.cur()
  if (isTRUE(capabilities("cairo"))) {
    cairo_pdf(pdf_path, width = W, height = H, family = FONT)
  } else {
    pdf(pdf_path, width = W, height = H, family = "Helvetica", useDingbats = FALSE,
        title = "Figure 3. The accommodation regime")
  }
  on.exit(if (dev.cur() != dev_before && dev.cur() > 1) dev.off(), add = TRUE)
  draw_figure()
  dev.off()
  report <- c(
    "Figure 3 regenerated with the original deterministic model and parameters.",
    sprintf("Threshold cP* = %.12f", cP_star),
    sprintf("Panel c: H = %.12f; coalition (D,P,A) = %s", Hfun(cP_c),
            paste(sprintf("%.12f", xdag(cP_c)), collapse = ", ")),
    sprintf("Panel d: H = %.12f; interior (D,P,A) = %s", Hfun(cP_d),
            paste(sprintf("%.12f", xint(cP_d)), collapse = ", ")),
    sprintf("Panel d: saddle coalition (D,P,A) = %s",
            paste(sprintf("%.12f", xdag(cP_d)), collapse = ", ")),
    "Checks passed: signs, equilibrium residuals, direct linear solve.",
    sprintf("RK4 step check (raw proposals): %d negative proposals (0 required); max drift of the coordinate sum %.1e (gate 1e-9); smallest proposed component %.1e; components clipped at 1e-300: %d.",
            n_neg, max_drift, min_raw, n_clip),
    "Integration: fixed-step RK4, dt=0.01; horizons 60 (c), 120 (d).",
    "Plot geometry is subsampled by arc length only; computed trajectories are unchanged."
  )
  writeLines(report, file.path(out_dir, "figure3_checks.txt"))
  cat(paste(report, collapse = "\n"), "\n")
  invisible(list(pdf = pdf_path, paths_c = paths_c, paths_d = paths_d))
}

make_figure3(out_dir)
