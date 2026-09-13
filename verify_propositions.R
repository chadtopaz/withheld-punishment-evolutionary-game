## verify_propositions.R
## ------------------------------------------------------------------------------
## Numerical checks of the analytical claims in the revised manuscript (base R only).
##   1. Referee 1's check of the exploitation-regime interior equilibrium (Fig. 2a).
##   2. Accommodation regime (Proposition 2): over random parameter draws,
##      (i)   an interior equilibrium exists iff H > 0, and equals (N_D, N_P, N_A)/S;
##      (ii)  at it, tr J = -pibar(x*) and det J = S x_D x_P x_A (so it is a sink);
##      (iii) the transverse exponent of defenders at the coalition equals H/(c_P + c_A);
##      (iv)  the disruptor share at coexistence is below its coalition value.
##   3. Global convergence by direct integration in each case (H > 0, the band, H = 0,
##      and strong accommodation), and at H = 0 the algebraic rate x_D ~ L0/(S0 t) of
##      the center-manifold calculation in the text (L0 = c_P + c_A, S0 = S at H = 0).
## Usage:  source("verify_propositions.R")      (runtime: about a minute)
##
## Gates.  The exact identities in 1 and 2 (existence iff H > 0, the closed form, the
## trace and determinant identities, the sink classification, the disruptor-share
## inequality, the transverse exponent) stop the script if violated, with the
## tolerances stated at each check.  The integrations in 3 are finite-time
## diagnostics: their deviations and the algebraic-rate ratio are gated only loosely,
## and the RK4 steps are checked for negative proposals and drift of the coordinate
## sum before the clip-and-renormalize step, so that the step itself is verified
## rather than the property it enforces.
## ------------------------------------------------------------------------------
set.seed(1)
## Output location: the folder containing this script when it is source()d, otherwise the
## working directory.  Run from the "npj Revision 1" folder so that the two coincide.
out_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(out_dir) || length(out_dir) != 1 || is.na(out_dir) || !nzchar(out_dir)) out_dir <- getwd()
cat("Writing outputs to:", out_dir, "\n")
## Everything printed below is also written to verify_report.txt in that folder.  The
## sink this script opens is closed on exit, whether the script finishes, stops at a
## gate, or is interrupted; sinks that were open before are left alone.
.sink_depth0 <- sink.number()
sink(file.path(out_dir, "verify_report.txt"), split = TRUE)
tryCatch({

rhs <- function(x, M) { p <- as.vector(M %*% x); x * (p - sum(x * p)) }
## RK4 integration to time T.  Besides the end state, the result records what the raw
## RK4 update proposed before the clip at 1e-300 (underflow protection) and the
## renormalization: the smallest proposed component, the number of steps with a
## negative proposal, the largest drift of the coordinate sum from one, and the number
## of components the clip changed.
integrate_end <- function(M, x0, T, dt = 0.01) {
  x <- x0 / sum(x0); min_raw <- Inf; n_neg <- 0L; max_drift <- 0; n_clip <- 0L
  for (s in 1:ceiling(T / dt)) {
    k1 <- rhs(x, M); k2 <- rhs(x + 0.5 * dt * k1, M); k3 <- rhs(x + 0.5 * dt * k2, M); k4 <- rhs(x + dt * k3, M)
    x_raw <- x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    min_raw <- min(min_raw, min(x_raw)); if (any(x_raw < 0)) n_neg <- n_neg + 1L
    max_drift <- max(max_drift, abs(sum(x_raw) - 1)); n_clip <- n_clip + sum(x_raw < 1e-300)
    x <- pmax(x_raw, 1e-300); x <- x / sum(x)
  }
  list(x = x, min_raw = min_raw, n_neg = n_neg, max_drift = max_drift, n_clip = n_clip)
}
interior_eq <- function(M) {
  A <- rbind(M[1, ] - M[2, ], M[2, ] - M[3, ], c(1, 1, 1))
  as.vector(solve(A, c(0, 0, 1)))
}
## Jacobian of the planar system in the (x_P, x_A) chart, by central differences
jac_chart <- function(M, x, h = 1e-6) {
  f <- function(y) { xx <- c(1 - y[1] - y[2], y[1], y[2]); p <- as.vector(M %*% xx); (xx * (p - sum(xx * p)))[2:3] }
  y0 <- x[2:3]; J <- matrix(0, 2, 2)
  for (j in 1:2) { e <- c(0, 0); e[j] <- h; J[, j] <- (f(y0 + e) - f(y0 - e)) / (2 * h) }
  J
}

## ---- 1. exploitation regime, Figure 2a parameters ------------------------------
M1 <- matrix(c(0, -1, 0.5,  1, 0, -3,  -1, 1, 0), 3, byrow = TRUE)
xs <- interior_eq(M1); pis <- as.vector(M1 %*% xs); J <- jac_chart(M1, xs)
cat(sprintf("Exploitation regime (Fig. 2a): x* = (%.3f, %.3f, %.3f); payoffs at x* = (%.3f, %.3f, %.3f)\n",
            xs[1], xs[2], xs[3], pis[1], pis[2], pis[3]))
## S is the common denominator of the interior equilibrium, eq. (9) of the manuscript
aD <- 1; bD <- 0.5; aP <- 3; bP <- 1; aA <- 1; bA <- 1
S1 <- (aP * aD + aP * bA + bA * bD) + (aA * aP + aA * bD + bP * bD) + (aA * aD + aD * bP + bA * bP)
cat(sprintf("  tr J = %.4f (formula (a_P a_A a_D - b_P b_A b_D)/S gives %.4f); det J = %.4f (formula S x_D x_P x_A gives %.4f); focus: %s\n",
            sum(diag(J)), (aP * aA * aD - bP * bA * bD) / S1, det(J), S1 * prod(xs),
            if (sum(diag(J))^2 < 4 * det(J)) "yes" else "no"))

## ---- 2. accommodation regime: random parameter draws ---------------------------
n_trials <- 20000; n_int <- 0; n_sink <- 0; n_focus <- 0; n_band <- 0
worst_tr <- 0; worst_det <- 0; worst_x <- 0; worst_lam <- 0; n_share_viol <- 0
for (t in 1:n_trials) {
  p <- runif(6, 0.05, 5); aD <- p[1]; bD <- p[2]; bP <- p[3]; aA <- p[4]; cP <- p[5]; cA <- p[6]
  M <- matrix(c(0, -aD, bD,  bP, 0, cP,  -aA, cA, 0), 3, byrow = TRUE)
  H <- bD * cA - cP * (cA + aD)
  ND <- H; NP <- bD * (aA + bP) - aA * cP; NA_ <- aD * (aA + bP) + bP * cA; S <- ND + NP + NA_
  ## direct solve of the equal-payoff system (may be singular when S = 0; treat as no interior eq)
  xs <- tryCatch(interior_eq(M), error = function(e) c(NA, NA, NA))
  interior <- !any(is.na(xs)) && all(xs > 0)
  if (interior != (H > 0)) stop("existence mismatch at trial ", t)
  if (H > 0) {
    n_int <- n_int + 1
    if (NP <= 0) stop("N_P not positive although H > 0")
    xf <- c(ND, NP, NA_) / S
    worst_x <- max(worst_x, max(abs(xf - xs)))
    J <- jac_chart(M, xs); pibar <- as.numeric(t(xs) %*% M %*% xs)
    trJ <- sum(diag(J)); detJ <- det(J)
    worst_tr <- max(worst_tr, abs(trJ + pibar)); worst_det <- max(worst_det, abs(detJ - S * prod(xs)))
    if (trJ < 0 && detJ > 0) n_sink <- n_sink + 1
    if (trJ^2 < 4 * detJ) n_focus <- n_focus + 1
    if (xs[3] >= cA / (cP + cA)) n_share_viol <- n_share_viol + 1
  } else {
    xd <- c(0, cP / (cP + cA), cA / (cP + cA)); pd <- as.vector(M %*% xd)
    lam <- pd[1] - sum(xd * pd)
    worst_lam <- max(worst_lam, abs(lam - H / (cP + cA)))
    if (bD * cA / (cA + aD) < cP && cP <= bD) n_band <- n_band + 1
  }
}
cat(sprintf("\nAccommodation regime, %d random parameter draws:\n", n_trials))
cat(sprintf("  interior equilibrium exists iff H > 0: verified; %d draws with H > 0\n", n_int))
cat(sprintf("  closed form (N_D, N_P, N_A)/S matches direct solve: max error %.1e (gate 1e-8)\n", worst_x))
cat(sprintf("  tr J = -pibar(x*): max error %.1e;  det J = S x_D x_P x_A: max error %.1e (gate 1e-6, central differences with h = 1e-6)\n", worst_tr, worst_det))
cat(sprintf("  all interior equilibria are sinks: %s (%d foci, %d nodes)\n", n_sink == n_int, n_focus, n_int - n_focus))
cat(sprintf("  disruptor share at coexistence below coalition value in all cases: %s\n", n_share_viol == 0))
cat(sprintf("  transverse exponent of D at coalition = H/(c_P + c_A): max error %.1e (gate 1e-10)\n", worst_lam))
cat(sprintf("  draws in the previously open band (H < 0 and c_P <= b_D): %d\n", n_band))
if (worst_x >= 1e-8) stop("closed-form equilibrium check failed")
if (worst_tr >= 1e-6 || worst_det >= 1e-6) stop("trace or determinant identity check failed")
if (n_sink != n_int) stop("an interior equilibrium is not a sink")
if (n_share_viol > 0) stop("disruptor-share inequality violated")
if (worst_lam >= 1e-10) stop("transverse exponent check failed")
cat("  all exact-identity gates passed\n")

## ---- 3. global convergence spot checks ----------------------------------------
cases <- list(
  "H > 0 (coexistence)"        = c(aD = 1, bD = 1, bP = 1, aA = 1, cP = 0.3, cA = 1),
  "band: H < 0, c_P <= b_D"    = c(aD = 1, bD = 1, bP = 1, aA = 1, cP = 0.7, cA = 1),
  "H = 0 exactly"              = c(aD = 1, bD = 1, bP = 1, aA = 1, cP = 0.5, cA = 1),
  "strong accommodation c_P > b_D" = c(aD = 1, bD = 1, bP = 1, aA = 1, cP = 2, cA = 1))
cat("\nGlobal convergence (6 random initial conditions per case; finite-time diagnostics, loosely gated):\n")
n_neg_all <- 0L; drift_all <- 0; clip_all <- 0L; min_raw_all <- Inf
for (nm in names(cases)) {
  p <- cases[[nm]]
  M <- matrix(c(0, -p["aD"], p["bD"],  p["bP"], 0, p["cP"],  -p["aA"], p["cA"], 0), 3, byrow = TRUE)
  H <- p["bD"] * p["cA"] - p["cP"] * (p["cA"] + p["aD"])
  target <- if (H > 0) interior_eq(M) else c(0, p["cP"] / (p["cP"] + p["cA"]), p["cA"] / (p["cP"] + p["cA"]))
  dev <- 0; T_end <- if (abs(H) > 1e-12) 600 else 3000; rate_ratio <- numeric(0)
  for (k in 1:6) {
    x0 <- rgamma(3, 1); x0 <- x0 / sum(x0)
    r <- integrate_end(M, x0, T = T_end); xe <- r$x
    n_neg_all <- n_neg_all + r$n_neg; drift_all <- max(drift_all, r$max_drift); clip_all <- clip_all + r$n_clip; min_raw_all <- min(min_raw_all, r$min_raw)
    dev <- max(dev, max(abs(xe - target)))
    if (abs(H) < 1e-12) {
      ## algebraic rate: x_D(T) * T * S0/L0 should approach 1
      L0 <- p["cP"] + p["cA"]
      S0 <- (p["bD"] * (p["aA"] + p["bP"]) - p["aA"] * p["cP"]) + (p["aD"] * (p["aA"] + p["bP"]) + p["bP"] * p["cA"])
      rate_ratio <- c(rate_ratio, xe[1] * T_end * S0 / L0)
    }
  }
  cat(sprintf("  %-32s H = %+.2f  target = (%.3f, %.3f, %.3f)  max deviation = %.1e%s\n", nm, H, target[1], target[2], target[3], dev,
              if (abs(H) < 1e-12) "  (non-hyperbolic: convergence is algebraic, so slower; gate 1e-2)" else "  (gate 1e-6)"))
  if (abs(H) < 1e-12) {
    cat(sprintf("  %-32s algebraic rate check: x_D(T) T S0/L0 at T = %d is %s (1 expected; gate 0.8 to 1.05)\n", "", T_end,
                paste(formatC(rate_ratio, format = "f", digits = 3), collapse = ", ")))
    if (dev >= 1e-2) stop("convergence diagnostic failed in the H = 0 case")
    if (any(rate_ratio < 0.8 | rate_ratio > 1.05)) stop("algebraic-rate diagnostic failed in the H = 0 case")
  } else if (dev >= 1e-6) stop("convergence diagnostic failed in case: ", nm)
}
cat(sprintf("  RK4 step check over all integrations: %d steps with a negative proposal (0 required), max drift of the coordinate sum before renormalization %.1e (gate 1e-9), smallest proposed component %.1e, components clipped at 1e-300: %d\n",
            n_neg_all, drift_all, min_raw_all, clip_all))
if (n_neg_all > 0 || drift_all >= 1e-9) stop("RK4 step check failed")
cat("All checks passed.\n")

}, finally = { while (sink.number() > .sink_depth0) sink() })
