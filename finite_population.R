## finite_population.R
## ------------------------------------------------------------------------------
## Finite-population check of the resurgence mechanism (Section 4.3, Figure 4).
##
## Process: the "local update" (linear pairwise-comparison) process of Traulsen,
## Claussen & Hauert (2005, Phys. Rev. Lett. 95, 238701).  In each elementary step
## a focal individual a and a role model b are drawn at random from a population of
## N actors; a adopts b's strategy with probability
##     p_{a->b} = 1/2 + w (pi_b - pi_a) / (2 dpi_max),
## where pi are the mean-field payoffs (M x) and dpi_max is the largest payoff
## difference that can occur, so that p lies in [0,1].  With probability mu the focal
## individual instead adopts a uniformly random strategy (reintroduction/mutation).
## One generation is N elementary steps.  As N -> infinity with mu = 0 the process
## converges to the replicator equation with time rescaled by w/dpi_max, so
##     replicator time  =  (w/dpi_max) * generations.
##
## Per elementary step, the probability that a type-i individual becomes type j is
##     r_ij = x_i [ (1 - mu) x_j p_{i->j} + mu/3 ].
## Because p_{i->j} + p_{j->i} = 1, the probability of a change within the unordered
## pair {i, j} is  R_ij = (1 - mu) x_i x_j + (mu/3)(x_i + x_j),  and the total
## probability of a change is
##     q = (1 - mu)(x_D x_P + x_D x_A + x_P x_A) + 2 mu / 3,
## which needs no payoffs.  The simulation is exact and event-driven: the holding
## time until the next change is geometric with parameter q (a long sojourn near a
## vertex costs one step), the pair is drawn with probability R_ij / q, and the
## direction i -> j with probability r_ij / R_ij, for which the single payoff
## difference (M_j - M_i) x is evaluated.  Before the experiments start the script
## (i) checks this decomposition against the six ordered rates computed directly
## from the definition at 307 states and four values of mu, with a tolerance gate;
## (ii) checks the sampler in use by Monte Carlo, comparing the first transition
## and the mean holding time at one interior state with their exact values; and
## (iii) runs the C++ and R versions from the same seed and requires that they agree
## on the sampled trajectory, the number of changes, and the number of episodes in
## that test.
##
## Experiment A also counts the disruptor-majority episodes (x_A above 1/2 after having
## been below 1/10) that occur before the population freezes.
##
## Implementation: the single-chain simulator is written in C++ (Rcpp), compiled once
## by this session from finitepop_sim.cpp (written next to this script) into a cache
## directory that the parallel workers then load without recompiling.  Every run
## (each replicate of each parameter set) is one task; the tasks are ordered longest
## first, grouped into bundles of a few seconds' work, and handed one bundle at a
## time to detectCores() - 1 workers with furrr, so all replicates of a parameter
## set run in parallel and no worker is tied to one parameter set.  A C++ compiler
## is required: if compilation fails, or a worker cannot load the compiled library,
## the script stops rather than substituting another simulator or design.  A pure-R
## simulator with the same semantics is kept only for the same-seed check (iii).
##
## Gates: the script stops if (i) or (iii) fails, or if any experiment-A run did not
## end with one strategy holding the whole population (so that the fixation and
## first-loss percentages quoted as shares of all runs are exactly that).  The Monte
## Carlo check (ii) is statistical and only warns.
##
## Packages (all on CRAN): Rcpp, future, furrr.  Missing ones are installed.
##
## Outputs (written next to this script):
##   fig_finitepop.pdf      Figure 4
##   results_finitepop.tex  LaTeX macros with the numbers quoted in the text
##   finitepop_summary.csv  per-experiment summary table
##   finitepop_runs.csv     one row per run (all raw results; at_vertex says whether the
##                          run ended with one strategy holding the whole population,
##                          absorbed is at_vertex restricted to mu = 0, where a vertex is
##                          absorbing)
##   finitepop_log.txt      copy of the console output
##   finitepop_sim.cpp      the C++ simulator (written by this script)
##
## Figures are written with cairo_pdf() when available so that the fonts are embedded.
##
## Usage:  source("finite_population.R")
## ------------------------------------------------------------------------------

set.seed(20260905)
## Output location: the folder containing this script when it is source()d, otherwise the
## working directory.  Run from the "npj Revision 1" folder so that the two coincide.
out_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(out_dir) || length(out_dir) != 1 || is.na(out_dir) || !nzchar(out_dir)) out_dir <- getwd()
cat("Writing outputs to:", out_dir, "\n")
## Everything printed below is also written to finitepop_log.txt in that folder.  The
## sink and any graphics device this script opens are closed on exit, whether the
## script finishes, stops at a gate, or is interrupted; sinks and devices that were
## open before are left alone.
.sink_depth0 <- sink.number(); .dev0 <- dev.cur()
sink(file.path(out_dir, "finitepop_log.txt"), split = TRUE)
.cleanup <- function() {
  if (dev.cur() != .dev0 && dev.cur() > 1) dev.off()
  if (requireNamespace("future", quietly = TRUE)) try(future::plan(future::sequential), silent = TRUE)
  while (sink.number() > .sink_depth0) sink()
}
tryCatch({
t_start <- proc.time()[["elapsed"]]
cat(R.version.string, "\n")

## ---- packages -----------------------------------------------------------------
pkgs <- c("Rcpp", "future", "furrr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({ library(future); library(furrr) })
cat(sprintf("Rcpp %s, future %s, furrr %s\n", packageVersion("Rcpp"), packageVersion("future"), packageVersion("furrr")))
n_workers <- max(1L, parallel::detectCores() - 1L)
cat("Parallel workers:", n_workers, "\n")

## ---- exploitation-regime parameters of Figure 2a ---------------------------
aD <- 1; bD <- 0.5; aP <- 3; bP <- 1; aA <- 1; bA <- 1
M1 <- matrix(c( 0,  -aD,  bD,
               bP,    0, -aP,
              -aA,   bA,   0), nrow = 3, byrow = TRUE)   # rows/cols ordered (D, P, A)
w <- 1
dpi_max <- max(apply(M1, 2, function(col) max(col) - min(col)))   # = 3.5 here
scale <- w / dpi_max                                              # replicator time per generation
stopifnot(aP * aA * aD > bP * bA * bD)                            # attracting-cycle condition (11)

## interior equilibrium of the replicator dynamics (equal payoffs)
interior_eq <- function(M) {
  A <- rbind(M[1, ] - M[2, ], M[2, ] - M[3, ], c(1, 1, 1))
  as.vector(solve(A, c(0, 0, 1)))
}
xstar <- interior_eq(M1)
cat(sprintf("Interior equilibrium x* = (%.3f, %.3f, %.3f); dpi_max = %.2f; replicator time per generation = %.4f\n",
            xstar[1], xstar[2], xstar[3], dpi_max, scale))

## initial counts: the interior equilibrium rounded to integers summing to N
init_counts <- function(N) {
  n0 <- floor(N * xstar)
  n0[3] <- N - n0[1] - n0[2]
  as.integer(n0)
}

## ordered transitions (i -> j), table `trans`: a type-i individual adopts strategy j
trans <- rbind(c(1, 2), c(1, 3), c(2, 1), c(2, 3), c(3, 1), c(3, 2))

## ==============================================================================
## Single-chain simulator, C++ version.  Signature (identical to the R version):
##   sim(M, N, mu, w, dpi_max, T_gen, n0, stop_at_vertex, track_episodes, lo, hi,
##       record_dt, max_events)
## T_gen is the horizon in generations; stop_at_vertex ends the run when one strategy
## holds the whole population; track_episodes switches the resurgence bookkeeping on;
## record_dt > 0 records the composition every record_dt generations; max_events > 0
## ends the run after that many changes (used only by the checks below).
## Returns a list: t_ext (generations at which a strategy was first lost, NA if never),
##   lost_first (1 = D, 2 = P, 3 = A, NA if none lost), at_vertex (TRUE if the run ended
##   with one strategy holding the whole population; with mu = 0 and stop_at_vertex this
##   is absorption, with mu > 0 it is only the state at the end), t_vertex (generations
##   at which the run ended if at_vertex, NA otherwise), final (counts), n_ep
##   (disruptor resurgence episodes: x_A rises above hi after having fallen below lo),
##   time_hi (elementary steps with x_A > hi), t_run (generations simulated), n_events
##   (changes applied), rec (recorded series with columns gen, xD, xP, xA when
##   record_dt > 0, else NULL).
## event_probs_*(M, N, mu, w, dpi_max, n) returns the six per-step transition
## probabilities in the row order of `trans`, computed through the pair decomposition
## exactly as the sampler uses it (backward direction as R_ij minus the forward rate).
## ==============================================================================
cpp_file <- file.path(out_dir, "finitepop_sim.cpp")
writeLines(c(
"#include <Rcpp.h>",
"#include <cmath>",
"#include <vector>",
"using namespace Rcpp;",
"",
"// Unordered pairs k = 0, 1, 2 are (i, j) = (D, P), (D, A), (P, A) with i < j.",
"static const int PI_[3] = {0, 0, 1};",
"static const int PJ_[3] = {1, 2, 2};",
"// Position of the ordered transitions i -> j and j -> i in the R table `trans`",
"// (rows 1->2, 1->3, 2->1, 2->3, 3->1, 3->2; 0-based here).",
"static const int FWD[3] = {0, 1, 3};",
"static const int BWD[3] = {2, 4, 5};",
"",
"struct Rates {",
"  double dM[3][3];   // dM[k][m] = M(j, m) - M(i, m), so that pi_j - pi_i = dM[k] . x",
"  double om, mu3, sel;",
"  Rates(const NumericMatrix& M, double mu, double w, double dpi_max) {",
"    for (int k = 0; k < 3; k++) for (int m = 0; m < 3; m++) dM[k][m] = M(PJ_[k], m) - M(PI_[k], m);",
"    om = 1.0 - mu; mu3 = mu / 3.0; sel = 0.5 * w / dpi_max;",
"  }",
"  // probability per step of a change within pair k",
"  inline double pair_weight(int k, const double* x) const {",
"    return om * x[PI_[k]] * x[PJ_[k]] + mu3 * (x[PI_[k]] + x[PJ_[k]]);",
"  }",
"  // probability per step that a type-i individual becomes type j (pair k, i < j)",
"  inline double forward_rate(int k, const double* x) const {",
"    double dpi = dM[k][0] * x[0] + dM[k][1] * x[1] + dM[k][2] * x[2];",
"    double p = 0.5 + sel * dpi;",
"    return x[PI_[k]] * (om * x[PJ_[k]] * p + mu3);",
"  }",
"};",
"",
"// [[Rcpp::export]]",
"NumericVector event_probs_cpp(NumericMatrix M, int N, double mu, double w, double dpi_max, IntegerVector n) {",
"  Rates rt(M, mu, w, dpi_max);",
"  double x[3]; for (int i = 0; i < 3; i++) x[i] = (double)n[i] / (double)N;",
"  NumericVector out(6);",
"  for (int k = 0; k < 3; k++) {",
"    double R = rt.pair_weight(k, x), f = rt.forward_rate(k, x);",
"    out[FWD[k]] = f; out[BWD[k]] = R - f;",
"  }",
"  return out;",
"}",
"",
"// Event-driven local update process with reintroduction for a three-strategy game.",
"// [[Rcpp::export]]",
"List simulate_chain_cpp(NumericMatrix M, int N, double mu, double w, double dpi_max,",
"                        double T_gen, IntegerVector n0, bool stop_at_vertex, bool track_episodes,",
"                        double lo, double hi, double record_dt, double max_events) {",
"  Rates rt(M, mu, w, dpi_max);",
"  int n[3] = {(int)n0[0], (int)n0[1], (int)n0[2]};",
"  const double Nd = (double)N;",
"  const double max_steps = T_gen * Nd;         // horizon in elementary steps",
"  double t_steps = 0.0, n_events = 0.0;",
"  double t_ext = NA_REAL; int lost_first = NA_INTEGER;",
"  bool armed = false; int n_ep = 0; double time_hi = 0.0;",
"  const bool recording = record_dt > 0;",
"  const double rec_step = record_dt * Nd;",
"  double next_rec = 0.0;",
"  std::vector<double> rec_t, rec_D, rec_P, rec_A;",
"  double x[3], R[3];",
"  while (true) {",
"    for (int i = 0; i < 3; i++) x[i] = (double)n[i] / Nd;",
"    for (int k = 0; k < 3; k++) R[k] = rt.pair_weight(k, x);",
"    double q = R[0] + R[1] + R[2];",
"    double hold;",
"    if (q <= 0.0) hold = R_PosInf;",
"    else { double u = R::runif(0.0, 1.0); hold = std::floor(std::log(u) / std::log1p(-q)) + 1.0; }",
"    double t_next = t_steps + hold;            // step at which the next change occurs",
"    if (recording) {                           // the state is constant on [t_steps, t_next)",
"      while (next_rec < t_next && next_rec <= max_steps) {",
"        rec_t.push_back(next_rec / Nd); rec_D.push_back(x[0]); rec_P.push_back(x[1]); rec_A.push_back(x[2]);",
"        next_rec += rec_step;",
"      }",
"    }",
"    if (track_episodes && x[2] > hi) time_hi += std::min(hold, max_steps - t_steps);",
"    if (t_next > max_steps) { t_steps = max_steps; break; }   // no further change within the horizon",
"    t_steps = t_next; n_events += 1.0;",
"    // draw the pair with probability R_k / q, then the direction; v is uniform on [0, q)",
"    double v = R::runif(0.0, 1.0) * q;",
"    int k; double cum;",
"    if (v < R[0]) { k = 0; cum = 0.0; }",
"    else if (v < R[0] + R[1]) { k = 1; cum = R[0]; }",
"    else { k = 2; cum = R[0] + R[1]; }",
"    double v_res = v - cum;",
"    if (R[k] <= 0.0) {                         // reachable only through rounding at the top of [0, q)",
"      k = 0; if (R[1] > R[k]) k = 1; if (R[2] > R[k]) k = 2;",
"      v_res = 0.0;",
"    }",
"    if (v_res < 0.0) v_res = 0.0;",
"    const int i = PI_[k], j = PJ_[k];",
"    bool forward = v_res < rt.forward_rate(k, x);",
"    if (forward && n[i] == 0) forward = false;         // guards against rounding at the boundary",
"    else if (!forward && n[j] == 0) forward = true;",
"    if (forward) { n[i] -= 1; n[j] += 1; } else { n[j] -= 1; n[i] += 1; }",
"    if (n[0] < 0 || n[1] < 0 || n[2] < 0 || n[0] > N || n[1] > N || n[2] > N) stop(\"count out of range\");",
"    if (lost_first == NA_INTEGER) {",
"      for (int m = 0; m < 3; m++) if (n[m] == 0) { t_ext = t_steps / Nd; lost_first = m + 1; break; }",
"    }",
"    if (track_episodes) {",
"      double xA = (double)n[2] / Nd;",
"      if (!armed && xA < lo) armed = true;",
"      else if (armed && xA > hi) { armed = false; n_ep += 1; }",
"    }",
"    if (stop_at_vertex && (n[0] == N || n[1] == N || n[2] == N)) break;",
"    if (max_events > 0 && n_events >= max_events) break;",
"  }",
"  const bool at_vertex = (n[0] == N || n[1] == N || n[2] == N);",
"  IntegerVector fin = IntegerVector::create(n[0], n[1], n[2]);",
"  RObject rec = R_NilValue;",
"  if (recording) {",
"    int m = rec_t.size(); NumericMatrix Rm(m, 4);",
"    for (int r = 0; r < m; r++) { Rm(r, 0) = rec_t[r]; Rm(r, 1) = rec_D[r]; Rm(r, 2) = rec_P[r]; Rm(r, 3) = rec_A[r]; }",
"    colnames(Rm) = CharacterVector::create(\"gen\", \"xD\", \"xP\", \"xA\");",
"    rec = Rm;",
"  }",
"  return List::create(Named(\"t_ext\") = t_ext, Named(\"lost_first\") = lost_first,",
"                      Named(\"at_vertex\") = at_vertex, Named(\"t_vertex\") = at_vertex ? t_steps / Nd : NA_REAL,",
"                      Named(\"final\") = fin, Named(\"n_ep\") = n_ep, Named(\"time_hi\") = time_hi,",
"                      Named(\"t_run\") = t_steps / Nd, Named(\"n_events\") = n_events, Named(\"rec\") = rec);",
"}"), cpp_file)

## Compile once, into a cache directory that the workers load from.  Fused
## multiply-add is switched off so that the C++ and R simulators round identically
## (this only matters for the same-seed check below).
cache_dir <- file.path(tempdir(), "finitepop_cache")
dir.create(cache_dir, showWarnings = FALSE)
old_cxxflags <- Sys.getenv("PKG_CXXFLAGS", unset = NA)
Sys.setenv(PKG_CXXFLAGS = paste(if (!is.na(old_cxxflags)) old_cxxflags else "", "-ffp-contract=off"))
t0 <- proc.time()[["elapsed"]]
compile_error <- tryCatch({ Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir); NULL }, error = function(e) conditionMessage(e))
if (is.na(old_cxxflags)) Sys.unsetenv("PKG_CXXFLAGS") else Sys.setenv(PKG_CXXFLAGS = old_cxxflags)
if (!is.null(compile_error))
  stop("C++ compilation failed (", compile_error, "). A C++ compiler is required; on macOS, install the Xcode command line tools (xcode-select --install).")
cat(sprintf("Simulator: C++ (Rcpp), compiled in %.0f s\n", proc.time()[["elapsed"]] - t0))

## ==============================================================================
## Single-chain simulator, pure-R version with the same signature, semantics, and
## order of floating-point operations.  Used only by the same-seed check (iii).
## ==============================================================================
simulate_chain_r <- function(M, N, mu, w, dpi_max, T_gen, n0, stop_at_vertex, track_episodes, lo, hi, record_dt, max_events) {
  PI_ <- c(1L, 1L, 2L); PJ_ <- c(2L, 3L, 3L)
  dM <- M[PJ_, ] - M[PI_, ]                     # dM[k, ] = M[j, ] - M[i, ]
  om <- 1 - mu; mu3 <- mu / 3; sel <- 0.5 * w / dpi_max
  n <- as.integer(n0); N <- as.integer(N); t_steps <- 0; max_steps <- T_gen * N; n_events <- 0
  t_ext <- NA_real_; lost_first <- NA_integer_; armed <- FALSE; n_ep <- 0L; time_hi <- 0
  recording <- record_dt > 0; rec <- NULL
  if (recording) { n_rec <- floor(T_gen / record_dt) + 1; rec <- matrix(NA_real_, n_rec, 4); kr <- 1; next_rec <- 0; rec_step <- record_dt * N }
  repeat {
    x <- n / N
    R1 <- om * x[1] * x[2] + mu3 * (x[1] + x[2])
    R2 <- om * x[1] * x[3] + mu3 * (x[1] + x[3])
    R3 <- om * x[2] * x[3] + mu3 * (x[2] + x[3])
    q <- R1 + R2 + R3
    hold <- if (q <= 0) Inf else floor(log(runif(1)) / log1p(-q)) + 1
    t_next <- t_steps + hold
    if (recording) while (kr <= n_rec && next_rec < t_next) { rec[kr, ] <- c(next_rec / N, x); kr <- kr + 1; next_rec <- next_rec + rec_step }
    if (track_episodes && x[3] > hi) time_hi <- time_hi + min(hold, max_steps - t_steps)
    if (t_next > max_steps) { t_steps <- max_steps; break }
    t_steps <- t_next; n_events <- n_events + 1
    v <- runif(1) * q
    if (v < R1) { k <- 1L; cum <- 0; Rk <- R1 } else if (v < R1 + R2) { k <- 2L; cum <- R1; Rk <- R2 } else { k <- 3L; cum <- R1 + R2; Rk <- R3 }
    v_res <- v - cum
    if (Rk <= 0) { k <- which.max(c(R1, R2, R3)); v_res <- 0 }
    if (v_res < 0) v_res <- 0
    i <- PI_[k]; j <- PJ_[k]
    dpi <- dM[k, 1] * x[1] + dM[k, 2] * x[2] + dM[k, 3] * x[3]
    forward <- v_res < x[i] * (om * x[j] * (0.5 + sel * dpi) + mu3)
    if (forward && n[i] == 0L) forward <- FALSE else if (!forward && n[j] == 0L) forward <- TRUE
    if (forward) { n[i] <- n[i] - 1L; n[j] <- n[j] + 1L } else { n[j] <- n[j] - 1L; n[i] <- n[i] + 1L }
    if (any(n < 0L) || any(n > N)) stop("count out of range")
    if (is.na(lost_first) && any(n == 0L)) { t_ext <- t_steps / N; lost_first <- which(n == 0L)[1] }
    if (track_episodes) { xA <- n[3] / N; if (!armed && xA < lo) armed <- TRUE else if (armed && xA > hi) { armed <- FALSE; n_ep <- n_ep + 1L } }
    if (stop_at_vertex && any(n == N)) break
    if (max_events > 0 && n_events >= max_events) break
  }
  if (recording) { rec <- rec[seq_len(kr - 1), , drop = FALSE]; colnames(rec) <- c("gen", "xD", "xP", "xA") }
  at_vertex <- any(n == N)
  list(t_ext = t_ext, lost_first = lost_first, at_vertex = at_vertex, t_vertex = if (at_vertex) t_steps / N else NA_real_,
       final = n, n_ep = n_ep, time_hi = time_hi, t_run = t_steps / N, n_events = n_events, rec = rec)
}
event_probs_r <- function(M, N, mu, w, dpi_max, n) {
  PI_ <- c(1L, 1L, 2L); PJ_ <- c(2L, 3L, 3L); FWD <- c(1L, 2L, 4L); BWD <- c(3L, 5L, 6L)
  dM <- M[PJ_, ] - M[PI_, ]; om <- 1 - mu; mu3 <- mu / 3; sel <- 0.5 * w / dpi_max
  x <- n / N; out <- numeric(6)
  for (k in 1:3) {
    i <- PI_[k]; j <- PJ_[k]
    Rk <- om * x[i] * x[j] + mu3 * (x[i] + x[j])
    dpi <- dM[k, 1] * x[1] + dM[k, 2] * x[2] + dM[k, 3] * x[3]
    f <- x[i] * (om * x[j] * (0.5 + sel * dpi) + mu3)
    out[FWD[k]] <- f; out[BWD[k]] <- Rk - f
  }
  out
}
## the six ordered transition probabilities straight from the definition of the process
event_probs_def <- function(M, N, mu, w, dpi_max, n) {
  x <- n / N; pi <- as.vector(M %*% x)
  p <- 0.5 + 0.5 * w * (pi[trans[, 2]] - pi[trans[, 1]]) / dpi_max
  x[trans[, 1]] * ((1 - mu) * x[trans[, 2]] * p + mu / 3)
}

## The compiled simulator in the current R session, as list(f = function, name = "C++").
## A worker copies the cache index written by the main session into its own temporary
## directory and loads the shared library the main session built; nothing is
## recompiled, and no two processes write to the same cache.  A worker that cannot
## load the library stops with an error, which future propagates to the main session.
get_sim <- function() {
  if (!exists("simulate_chain_cpp", envir = globalenv(), inherits = FALSE)) {
    local_cache <- file.path(tempdir(), "finitepop_cache")
    if (normalizePath(local_cache, mustWork = FALSE) != normalizePath(cache_dir, mustWork = FALSE)) {
      dir.create(local_cache, showWarnings = FALSE)
      file.copy(list.files(cache_dir, full.names = TRUE), local_cache, recursive = TRUE, overwrite = TRUE)
    }
    Rcpp::sourceCpp(cpp_file, cacheDir = local_cache)
  }
  list(f = get("simulate_chain_cpp", envir = globalenv()), name = "C++")
}

## ==============================================================================
## Checks of the event sampler against the definition of the process
## ==============================================================================
cat("\nChecks of the event sampler:\n")
## (i) the pair decomposition reproduces the six ordered transition probabilities
set.seed(1)
states <- list(c(10L, 0L, 0L), c(0L, 10L, 0L), c(0L, 0L, 10L), c(9L, 1L, 0L), c(0L, 1L, 9L), c(1L, 0L, 9L), c(5L, 5L, 0L))
for (t in 1:300) states[[length(states) + 1]] <- as.integer(rmultinom(1, sample(c(5L, 10L, 50L, 1000L), 1), runif(3)))
err_r <- 0; err_cpp <- 0
for (nv in states) for (mu_v in c(0, 1e-3, 0.05, 0.5)) {
  Nv <- sum(nv); ref <- event_probs_def(M1, Nv, mu_v, w, dpi_max, nv)
  err_r <- max(err_r, max(abs(event_probs_r(M1, Nv, mu_v, w, dpi_max, nv) - ref)))
  err_cpp <- max(err_cpp, max(abs(event_probs_cpp(M1, Nv, mu_v, w, dpi_max, nv) - ref)))
}
cat(sprintf("  six transition probabilities, pair decomposition vs definition, %d states x 4 values of mu: max error R %.1e, C++ %.1e\n",
            length(states), err_r, err_cpp))
if (err_r >= 1e-12 || err_cpp >= 1e-12) stop("event sampler check (i) failed: the pair decomposition does not reproduce the transition probabilities")
## (ii) Monte Carlo check of the sampler actually used: first transition and holding time
sim <- get_sim()
Nv <- 10L; nv <- c(3L, 4L, 3L); mu_v <- 0.05; n_mc <- 20000L
ref <- event_probs_def(M1, Nv, mu_v, w, dpi_max, nv); q_ref <- sum(ref); p_ref <- ref / q_ref
set.seed(2)
first <- t(vapply(seq_len(n_mc), function(s) {
  r <- sim$f(M1, Nv, mu_v, w, dpi_max, 1e6, nv, FALSE, FALSE, 0.1, 0.5, 0, 1)
  c(r$final - nv, r$t_run * Nv)
}, numeric(4)))
from <- apply(first[, 1:3], 1, function(d) which(d == -1L)[1]); to <- apply(first[, 1:3], 1, function(d) which(d == 1L)[1])
k_obs <- match(paste(from, to), paste(trans[, 1], trans[, 2]))
freq <- tabulate(k_obs, 6) / n_mc
z <- (freq - p_ref) / sqrt(p_ref * (1 - p_ref) / n_mc)
z_hold <- (mean(first[, 4]) - 1 / q_ref) / (sd(first[, 4]) / sqrt(n_mc))
cat(sprintf("  Monte Carlo, %s simulator, %d single-change runs at N = %d, n = (%d, %d, %d), mu = %g:\n", sim$name, n_mc, Nv, nv[1], nv[2], nv[3], mu_v))
cat(sprintf("    transition frequencies vs r_ij/q: max |z| = %.2f;  mean holding time %.3f steps vs 1/q = %.3f steps (z = %.2f)\n",
            max(abs(z)), mean(first[, 4]), 1 / q_ref, z_hold))
if (max(abs(z)) > 4.5 || abs(z_hold) > 4.5) cat("    WARNING: a deviation above 4.5 standard errors; the sampler should be inspected\n")
## (iii) the C++ and R simulators follow the same trajectory from the same seed (a
## deterministic identity, so any disagreement stops the script)
args <- list(M1, 60L, 1e-2, w, dpi_max, 300, init_counts(60L), FALSE, TRUE, 0.1, 0.5, 1, 0)
set.seed(3); a <- do.call(simulate_chain_cpp, args)
set.seed(3); b <- do.call(simulate_chain_r, args)
same_dim <- identical(dim(a$rec), dim(b$rec))
max_diff <- if (same_dim) max(abs(a$rec - b$rec)) else NA
cat(sprintf("  same seed, C++ vs R (N = 60, mu = 0.01, 300 generations, shares sampled once per generation): %d vs %d changes, %d vs %d episodes, max |difference| in sampled shares %s (0 expected)\n",
            as.integer(a$n_events), as.integer(b$n_events), a$n_ep, b$n_ep,
            if (same_dim) formatC(max_diff, format = "e", digits = 1) else "n/a (different lengths)"))
if (!same_dim || max_diff != 0 || a$n_events != b$n_events || a$n_ep != b$n_ep)
  stop("event sampler check (iii) failed: the C++ and R simulators disagree from the same seed")

## ==============================================================================
## Design.  Experiment A: no reintroduction (mu = 0), time to the first loss of a
## strategy and the state in which the population freezes, as a function of N.
## Experiment B: with reintroduction (mu > 0), pooled simulated time per resurgence
## (total simulated time in a cell divided by the number of resurgences; an inverse
## resurgence rate, stored under the name `period`, that estimates the mean interval
## between resurgences only when the process is stationary) and share of time with a
## disruptor majority over long runs.
## ==============================================================================
N_A <- c(100, 300, 1000, 3000, 10000, 30000, 100000); reps_A <- 1000
N_B <- c(100, 1000, 10000, 100000); mu_B <- c(1e-5, 1e-4, 1e-3, 1e-2); reps_B <- 20; T_B <- 2000
T_A <- 1e7                                       # horizon in generations for experiment A (never reached in practice)
tasks <- rbind(
  expand.grid(exp = "A", N = N_A, mu = 0,    rep = seq_len(reps_A), stringsAsFactors = FALSE),
  expand.grid(exp = "B", N = N_B, mu = mu_B, rep = seq_len(reps_B), stringsAsFactors = FALSE))
tasks$T_gen <- ifelse(tasks$exp == "A", T_A, T_B / scale)
## Rough cost proxy (expected number of changes) used only to order the tasks, longest
## first, and to group them: experiment A runs last about 50 + 80 log10(N/100)
## generations and about a fifth of their steps are changes; in experiment B the share of
## steps that are changes falls with mu.
tasks$cost <- ifelse(tasks$exp == "A", 0.21 * tasks$N * (50 + 80 * log10(tasks$N / 100)),
                     0.24 * tasks$N * tasks$T_gen * (tasks$mu / 1e-2)^0.17)
tasks <- tasks[order(-tasks$cost, tasks$exp, -tasks$N, -tasks$mu, tasks$rep), ]
rownames(tasks) <- NULL
tasks$id <- seq_len(nrow(tasks))
## One future per bundle of consecutive tasks worth about target_changes changes of work
## (a few seconds), so that the heaviest runs are dispatched singly while the many light
## runs do not each pay the dispatch overhead of a future.
target_changes <- 3e7
bundles <- list(); cur <- list(); cur_cost <- 0
for (i in seq_len(nrow(tasks))) {
  cur[[length(cur) + 1]] <- tasks[i, ]; cur_cost <- cur_cost + tasks$cost[i]
  if (cur_cost >= target_changes) { bundles[[length(bundles) + 1]] <- cur; cur <- list(); cur_cost <- 0 }
}
if (length(cur) > 0) bundles[[length(bundles) + 1]] <- cur
cat(sprintf("\nExperiment A: N in {%s}, %d runs each (mu = 0), run until one strategy holds the whole population\n",
            paste(N_A, collapse = ", "), reps_A))
cat(sprintf("Experiment B: N in {%s}, mu in {%s}, %d runs each of %d replicator time units (%.0f generations)\n",
            paste(N_B, collapse = ", "), paste(mu_B, collapse = ", "), reps_B, T_B, T_B / scale))
cat(sprintf("%d runs in total, longest first, in %d bundles of roughly equal work; one bundle at a time per worker\n",
            nrow(tasks), length(bundles)))

run_task <- function(task) {
  sim <- get_sim()
  N <- as.integer(task$N)
  t0 <- proc.time()[["elapsed"]]
  res <- sim$f(M1, N, task$mu, w, dpi_max, task$T_gen, init_counts(N),
               task$exp == "A", TRUE, 0.1, 0.5, 0, 0)
  list(id = task$id, exp = task$exp, N = task$N, mu = task$mu, rep = task$rep,
       t_ext = res$t_ext, lost_first = res$lost_first, at_vertex = res$at_vertex, t_vertex = res$t_vertex,
       final = res$final, n_ep = res$n_ep, time_hi = res$time_hi, t_run = res$t_run,
       n_events = res$n_events, elapsed = proc.time()[["elapsed"]] - t0, sim_used = sim$name)
}

plan(multisession, workers = n_workers)
## warm-up: every worker loads the simulator before the clock starts
invisible(future_map(seq_len(n_workers), function(i) { Sys.sleep(1); get_sim()$name },
                     .options = furrr_options(seed = TRUE, chunk_size = 1L)))
run_bundle <- function(bundle) lapply(bundle, run_task)
set.seed(20260905)                                # per-bundle seeds are derived from this state
t0 <- proc.time()[["elapsed"]]
results <- future_map(bundles, run_bundle,
                      .options = furrr_options(seed = TRUE, chunk_size = 1L),
                      .progress = FALSE)
results <- unlist(results, recursive = FALSE)     # one entry per run
plan(sequential)
wall <- proc.time()[["elapsed"]] - t0
cat(sprintf("Simulations finished in %.0f s\n", wall))

## ---- assemble results -------------------------------------------------------
res_df <- data.frame(
  id = vapply(results, function(r) r$id, integer(1)),
  exp = vapply(results, function(r) r$exp, character(1)),
  N = vapply(results, function(r) r$N, numeric(1)),
  mu = vapply(results, function(r) r$mu, numeric(1)),
  rep = vapply(results, function(r) r$rep, integer(1)),
  t_ext = vapply(results, function(r) as.numeric(r$t_ext), numeric(1)),
  lost_first = vapply(results, function(r) as.integer(r$lost_first), integer(1)),
  at_vertex = vapply(results, function(r) r$at_vertex, logical(1)),
  t_vertex = vapply(results, function(r) as.numeric(r$t_vertex), numeric(1)),
  final_vertex = vapply(results, function(r) if (r$at_vertex) which(r$final == r$N)[1] else NA_integer_, integer(1)),
  final_D = vapply(results, function(r) r$final[1], integer(1)),
  final_P = vapply(results, function(r) r$final[2], integer(1)),
  final_A = vapply(results, function(r) r$final[3], integer(1)),
  n_ep = vapply(results, function(r) as.integer(r$n_ep), integer(1)),
  time_hi = vapply(results, function(r) r$time_hi, numeric(1)),
  t_run = vapply(results, function(r) r$t_run, numeric(1)),
  n_events = vapply(results, function(r) r$n_events, numeric(1)),
  elapsed = vapply(results, function(r) r$elapsed, numeric(1)),
  sim_used = vapply(results, function(r) r$sim_used, character(1)))
## a vertex is absorbing only without reintroduction; with mu > 0, at_vertex is the state at the horizon
res_df$absorbed <- res_df$at_vertex & res_df$mu == 0
sim_tab <- table(res_df$sim_used)
cat("Simulator used by the workers:", paste(sprintf("%s in %d runs", names(sim_tab), sim_tab), collapse = ", "), "\n")
if (any(res_df$sim_used != "C++")) stop("a worker did not use the compiled simulator")

## Experiment A summary (replicator time units)
A <- res_df[res_df$exp == "A", ]
censored_A <- sum(!A$absorbed); noloss_A <- sum(is.na(A$t_ext))
cat(sprintf("\nExperiment A (mu = 0): %d of %d runs ended with one strategy holding the whole population; %d ended before any strategy was lost\n",
            sum(A$absorbed), nrow(A), noloss_A))
## Gate: the manuscript quotes the first-loss times and the fixation shares as summaries of
## all runs, so every experiment-A run must have lost a strategy and reached a vertex.
## If not (a changed design or horizon), the raw runs are saved for inspection under a
## separate name and nothing else is written.
if (censored_A > 0 || noloss_A > 0) {
  write.csv(res_df[order(res_df$id), ], file.path(out_dir, "finitepop_runs_INCOMPLETE.csv"), row.names = FALSE)
  stop(sprintf("%d experiment-A run(s) did not reach a vertex and %d had no strategy loss within the horizon; no outputs written (raw runs saved to finitepop_runs_INCOMPLETE.csv)", censored_A, noloss_A))
}
resA <- data.frame(N = N_A, ext_median = NA, ext_q25 = NA, ext_q75 = NA, ext_mean = NA,
                   frac_D = NA, frac_P = NA, frac_A = NA, lost_D = NA, lost_P = NA, lost_A = NA, censored = NA,
                   ep_median = NA, ep_mean = NA, ep_le1 = NA, ep_zero = NA)
for (r in seq_along(N_A)) {
  a <- A[A$N == N_A[r], ]; te <- a$t_ext * scale; fv <- a$final_vertex[a$absorbed]
  resA$ext_median[r] <- median(te, na.rm = TRUE); resA$ext_mean[r] <- mean(te, na.rm = TRUE)
  resA$ext_q25[r] <- quantile(te, 0.25, na.rm = TRUE); resA$ext_q75[r] <- quantile(te, 0.75, na.rm = TRUE)
  resA$frac_D[r] <- mean(fv == 1); resA$frac_P[r] <- mean(fv == 2); resA$frac_A[r] <- mean(fv == 3)
  resA$lost_D[r] <- mean(a$lost_first == 1, na.rm = TRUE); resA$lost_P[r] <- mean(a$lost_first == 2, na.rm = TRUE)
  resA$lost_A[r] <- mean(a$lost_first == 3, na.rm = TRUE); resA$censored[r] <- sum(!a$absorbed)
  ## disruptor-majority episodes (x_A above 1/2 after having been below 1/10) before the population freezes
  resA$ep_median[r] <- median(a$n_ep); resA$ep_mean[r] <- mean(a$n_ep)
  resA$ep_le1[r] <- mean(a$n_ep <= 1); resA$ep_zero[r] <- mean(a$n_ep == 0)
  cat(sprintf("  N = %6d: median first loss %6.1f (IQR %5.1f-%5.1f) replicator time; lost first D/P/A = %.2f/%.2f/%.2f; freezes at D/P/A = %.2f/%.2f/%.2f; not absorbed: %d\n",
              N_A[r], resA$ext_median[r], resA$ext_q25[r], resA$ext_q75[r], resA$lost_D[r], resA$lost_P[r], resA$lost_A[r],
              resA$frac_D[r], resA$frac_P[r], resA$frac_A[r], resA$censored[r]))
  cat(sprintf("              disruptor-majority episodes before freezing: median %g, mean %.2f, none in %.0f%% of runs, at most one in %.0f%% of runs\n",
              resA$ep_median[r], resA$ep_mean[r], 100 * resA$ep_zero[r], 100 * resA$ep_le1[r]))
}
ext_slope <- coef(lm(ext_median ~ log10(N), data = resA))[2]

## Experiment B summary
B <- res_df[res_df$exp == "B", ]
resB <- expand.grid(N = N_B, mu = mu_B)
resB$episodes <- NA; resB$period <- NA; resB$frac_hi <- NA
cat("\nExperiment B (mu > 0), T =", T_B, "replicator time units per run:\n")
for (row in seq_len(nrow(resB))) {
  b <- B[B$N == resB$N[row] & B$mu == resB$mu[row], ]
  total_ep <- sum(b$n_ep); total_T <- sum(b$t_run) * scale
  resB$episodes[row] <- total_ep
  resB$period[row]   <- if (total_ep > 0) total_T / total_ep else NA   # pooled time per resurgence (inverse rate), not a measured mean interval
  resB$frac_hi[row]  <- sum(b$time_hi) / (sum(b$t_run) * resB$N[row])   # steps with x_A > 1/2 over all steps
  cat(sprintf("  N = %6d, mu = %g: %5d episodes in %d runs; pooled time per resurgence %7.1f; share of time with x_A > 1/2: %.2f\n",
              resB$N[row], resB$mu[row], total_ep, nrow(b), resB$period[row], resB$frac_hi[row]))
}

sparse_threshold <- 50                            # cells with fewer resurgences in total are drawn with open symbols
i_min <- which.min(resB$episodes)
zero_runs_min <- sum(B$N == resB$N[i_min] & B$mu == resB$mu[i_min] & B$n_ep == 0)
cat(sprintf("  fewest resurgences: %d in %d runs at N = %d, mu = %g (%d of those runs had none); cells below %d are drawn with open symbols\n",
            resB$episodes[i_min], reps_B, resB$N[i_min], resB$mu[i_min], zero_runs_min, sparse_threshold))

## Run sizes and timings
cat("\nChanges simulated and wall time per run, by parameter set:\n")
sets <- unique(tasks[, c("exp", "N", "mu")])
for (s in seq_len(nrow(sets))) {
  d <- res_df[res_df$exp == sets$exp[s] & res_df$N == sets$N[s] & res_df$mu == sets$mu[s], ]
  cat(sprintf("  %s N = %6d mu = %-6g: %9.3g changes, %10.1f generations, %7.2f s per run (max %7.2f s)\n",
              sets$exp[s], sets$N[s], sets$mu[s], mean(d$n_events), mean(d$t_run), mean(d$elapsed), max(d$elapsed)))
}
cat(sprintf("Total: %.3g changes in %.0f s of worker time; wall time %.0f s on %d workers (parallel efficiency %.0f%%); %.3g changes per worker-second\n",
            sum(res_df$n_events), sum(res_df$elapsed), wall, n_workers, 100 * sum(res_df$elapsed) / (wall * n_workers),
            sum(res_df$n_events) / sum(res_df$elapsed)))

## ---- deterministic replicator trajectory (RK4), time in replicator units --------
## The clip at 1e-300 protects against underflow and the renormalization removes
## rounding drift; the raw RK4 proposal is checked before both (no negative proposal,
## drift of the coordinate sum below 1e-9, else the script stops), and the number of
## components the clip changed is reported.
replicator_traj <- function(M, x0, T, dt = 0.01) {
  rhs <- function(x) { p <- as.vector(M %*% x); x * (p - sum(x * p)) }
  n_steps <- ceiling(T / dt); out <- matrix(NA_real_, n_steps + 1, 4)
  x <- x0 / sum(x0); out[1, ] <- c(0, x)
  min_raw <- Inf; n_neg <- 0L; max_drift <- 0; n_clip <- 0L
  for (s in 1:n_steps) {
    k1 <- rhs(x); k2 <- rhs(x + 0.5 * dt * k1); k3 <- rhs(x + 0.5 * dt * k2); k4 <- rhs(x + dt * k3)
    x_raw <- x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    min_raw <- min(min_raw, min(x_raw)); if (any(x_raw < 0)) n_neg <- n_neg + 1L
    max_drift <- max(max_drift, abs(sum(x_raw) - 1)); n_clip <- n_clip + sum(x_raw < 1e-300)
    x <- pmax(x_raw, 1e-300); x <- x / sum(x)
    out[s + 1, ] <- c(s * dt, x)
  }
  cat(sprintf("Deterministic curve, RK4 step check (raw proposals): %d negative proposals (0 required); max drift of the coordinate sum %.1e (gate 1e-9); smallest proposed component %.1e; components clipped at 1e-300: %d\n",
              n_neg, max_drift, min_raw, n_clip))
  if (n_neg > 0 || max_drift >= 1e-9) stop("RK4 step check failed for the deterministic curve")
  colnames(out) <- c("t", "xD", "xP", "xA"); out
}

## ---- time-series panels: one run each at N = 1000, mu = 0 and mu = 1e-3 ---------
sim <- get_sim()
N_ts <- 1000L
## Panel (a) shows a run with the typical outcome at N = 1000 (the public lost first, the
## population frozen at the all-defender state within the window); seeds are tried in order
## until one gives it, and the seed used is reported.
panelA <- NULL
for (seed_a in 1:200) {
  set.seed(seed_a)
  r <- sim$f(M1, N_ts, 0, w, dpi_max, 150 / scale, init_counts(N_ts), FALSE, TRUE, 0.1, 0.5, 0.25, 0)
  if (!is.na(r$lost_first) && r$lost_first == 2 && r$at_vertex && r$final[1] == N_ts) { panelA <- r; break }
}
if (is.null(panelA)) { set.seed(1); panelA <- sim$f(M1, N_ts, 0, w, dpi_max, 150 / scale, init_counts(N_ts), FALSE, TRUE, 0.1, 0.5, 0.25, 0); seed_a <- NA }
rec0 <- panelA$rec
## time at which the population freezes: first record at which one share equals 1
frozen_rows <- which(apply(rec0[, 2:4], 1, max) >= 1)
t_freeze_a <- if (length(frozen_rows) > 0) rec0[frozen_rows[1], "gen"] * scale else NA
cat(sprintf("\nPanel (a): seed %s; public lost at %.1f replicator units; frozen at the all-defender state from %.1f; %d disruptor-majority episode(s) before freezing\n",
            as.character(seed_a), panelA$t_ext * scale, t_freeze_a, panelA$n_ep))
set.seed(4)
rec1 <- sim$f(M1, N_ts, 1e-3, w, dpi_max, 400 / scale, init_counts(N_ts), FALSE, TRUE,  0.1, 0.5, 0.25, 0)$rec
det  <- replicator_traj(M1, init_counts(N_ts) / N_ts, T = 150)

## ==============================================================================
## Figure 4
## ==============================================================================
## Drawn at 7.6 x 7.1 in and included at the text width (6.5 in), so that the axis
## titles print at about 10 pt and the tick labels, legends and secondary headings at
## about 9 pt; the height is the most the page allows next to the caption.  Three
## rows: the two single-run time series (a, b); an enlargement of three successive
## cycles of run (b) across the full width (c); and the two aggregate panels (d, e).
## Headings sit at the top left of each panel, bold, with the run conditions on a
## smaller second line where there are conditions to state.  The only text inside a
## panel is the three direct curve labels in the enlargement; no statistics are
## written inside the panels.
## The key to the curves of the time-series panels is one shared legend in the outer
## margin above the top row, so that no legend sits on the data.  The curve colors are
## the Okabe-Ito blue, black, and vermilion, which remain distinct in grayscale and
## under color-vision deficiency; simulated shares are solid, the deterministic
## disruptor share dashed.
colD <- "#0072B2"; colP <- "#000000"; colA <- "#D55E00"
lwdD <- 1.6; lwdP <- 1.6; lwdA <- 1.8

## Interval enlarged in panel (c): from just before one disruptor resurgence of run (b)
## to just before the third resurgence after it, so that three full cycles are shown.
## A resurgence starts when the disruptor share first exceeds 1/2 after having been
## below 1/10 (the definition used for the counts); the first resurgence starting after
## replicator time 120 is used, past the initial transient.
tt1 <- rec1[, "gen"] * scale; xa1 <- rec1[, "xA"]
starts <- numeric(0); armed <- FALSE
for (i in seq_along(xa1)) {
  if (xa1[i] < 0.1) armed <- TRUE
  if (armed && xa1[i] > 0.5) { starts <- c(starts, tt1[i]); armed <- FALSE }
}
k <- which(starts >= 120)[1]
if (!is.na(k) && length(starts) >= k + 3) {
  pad  <- 0.15 * (starts[k + 1] - starts[k])
  zoom <- round(c(starts[k] - pad, starts[k + 3] - pad))
} else {
  zoom <- c(150, 225)
}
cat(sprintf("Panel (c): enlargement of run (b) over replicator time %d to %d; %d resurgences start inside it\n",
            zoom[1], zoom[2], sum(starts >= zoom[1] & starts <= zoom[2])))

if (isTRUE(capabilities("cairo"))) cairo_pdf(file.path(out_dir, "fig_finitepop.pdf"), width = 7.6, height = 7.1, family = "Helvetica") else pdf(file.path(out_dir, "fig_finitepop.pdf"), width = 7.6, height = 7.1)
layout(matrix(c(1, 2, 3, 3, 4, 5), 3, 2, byrow = TRUE), heights = c(1, 0.85, 1))
par(cex = 0.83)                                    # the text scale of a two-row layout, kept so the panel text prints as before
par(oma = c(0, 0, 2.2, 0), mar = c(3.8, 5.0, 2.4, 1), mgp = c(2.6, 0.65, 0), cex.lab = 1.15, cex.axis = 1.05)
cex_panel <- par("cex")                            # reused for the shared legend
## panel label only; the journal requires titles to appear in the figure legend rather
## than in the image, so the descriptive heading and the run conditions are not drawn
heading <- function(label, main = NULL, sub = NULL) {
  mtext(label, side = 3, adj = 0, line = 0.6, font = 2, cex = 1.15)
}

## time-series panel; shade = c(t1, t2) marks the interval enlarged in panel (c)
ts_panel <- function(rec, tmax, shade = NULL) {
  t <- rec[, "gen"] * scale
  plot(t, rec[, "xA"], type = "n", xlim = c(0, tmax), ylim = c(0, 1), xlab = "time (replicator units)",
       ylab = "population share", las = 1)
  if (!is.null(shade)) { u <- par("usr"); rect(shade[1], u[3], shade[2], u[4], col = "grey90", border = NA) }
  lines(t, rec[, "xD"], col = colD, lwd = lwdD)
  lines(t, rec[, "xP"], col = colP, lwd = lwdP)
  lines(t, rec[, "xA"], col = colA, lwd = lwdA)
  box()                                            # redraw the frame over the shading
}
## (a) no reintroduction; the deterministic disruptor share from the same start is dashed
## (the shared legend above the row is drawn after panel (e))
ts_panel(rec0, 150)
lines(det[, "t"], det[, "xA"], col = colA, lwd = lwdA, lty = 2)
heading("(a)", "Without reintroduction", expression(paste(italic(N) == 1000, ",  ", mu == 0)))
## (b) with reintroduction; the shaded interval is enlarged in (c)
ts_panel(rec1, 400, shade = zoom)
heading("(b)", "With reintroduction", expression(paste(italic(N) == 1000, ",  ", mu == 10^-3)))

## (c) enlargement of the shaded interval of (b), full width, same share scale (the
## axis stops at 1; the small headroom above it holds the curve labels).  Each curve is
## labeled once, above its highest peak in the interval, where the other two shares
## are near zero.
sel <- tt1 >= zoom[1] & tt1 <= zoom[2]
par(mar = c(3.8, 5.0, 2.3, 1))                     # one heading line, so a shorter top margin
plot(tt1[sel], xa1[sel], type = "n", xlim = zoom, ylim = c(0, 1.12), yaxt = "n",
     xlab = "time (replicator units)", ylab = "population share", las = 1)
axis(2, at = seq(0, 1, 0.2), las = 1)
lines(tt1[sel], rec1[sel, "xD"], col = colD, lwd = lwdD)
lines(tt1[sel], rec1[sel, "xP"], col = colP, lwd = lwdP)
lines(tt1[sel], rec1[sel, "xA"], col = colA, lwd = lwdA)
label_peak <- function(col_name, label, col) {
  y <- rec1[sel, col_name]; i <- which.max(y)
  text(tt1[sel][i], y[i] + 0.075, label, col = col, font = 2, cex = 1.1, xpd = NA)
}
label_peak("xD", "D", colD); label_peak("xP", "P", colP); label_peak("xA", "A", colA)
heading("(c)", "Detail of (b), the shaded interval")
par(mar = c(3.8, 5.0, 2.4, 1))                     # restore the panel-label margin

## (d) first-loss time vs N (mu = 0): median with interquartile range.  The share of runs
## that freeze at the all-defender state is reported in the text, not on the panel.
plot(resA$N, resA$ext_median, log = "x", type = "n", ylim = c(0, max(resA$ext_q75) * 1.15),
     xlab = expression(paste("population size ", italic(N))), ylab = "time to first loss of a strategy", las = 1, xaxt = "n")
axis(1, at = c(100, 1000, 10000, 1e5), labels = expression(10^2, 10^3, 10^4, 10^5))
arrows(resA$N, resA$ext_q25, resA$N, resA$ext_q75, angle = 90, code = 3, length = 0.03, col = colD)
lines(resA$N, resA$ext_median, col = colD, lwd = 1.5); points(resA$N, resA$ext_median, pch = 16, col = colD, cex = 1.1)
heading("(d)", "Time to first strategy loss", expression(paste(mu == 0, ",  median and interquartile range")))

## (e) pooled time per resurgence vs mu (mu > 0).  One marker shape per N (square, circle,
## triangle, diamond) on a light-to-dark blue ramp; a cell with fewer than sparse_threshold
## resurgences in total keeps its shape but is drawn open.
cols_all <- c("#9ecae1", "#6baed6", "#3182bd", "#08519c"); pch_all <- c(22, 21, 24, 23)
cols_N <- cols_all[seq_along(N_B)]; pch_N <- pch_all[seq_along(N_B)]
plot(resB$mu, resB$period, log = "xy", type = "n", xlab = expression(paste("reintroduction rate ", mu)),
     ylab = "time per resurgence", las = 1, xaxt = "n",
     ylim = range(resB$period, na.rm = TRUE) * c(0.7, 1.4))
axis(1, at = mu_B, labels = expression(10^-5, 10^-4, 10^-3, 10^-2))
for (i in seq_along(N_B)) {
  s <- resB$N == N_B[i]
  lines(resB$mu[s], resB$period[s], col = cols_N[i], lwd = 1.5)
  points(resB$mu[s], resB$period[s], col = cols_N[i], pch = pch_N[i],
         bg = ifelse(resB$episodes[s] >= sparse_threshold, cols_N[i], "white"), cex = 1.25, lwd = 1.3)
}
legend("topright", legend = parse(text = sprintf("italic(N) == 10^%d", round(log10(N_B)))),
       col = cols_N, pt.bg = cols_N, pch = pch_N, lwd = 1.5, pt.cex = 1.25, bty = "n", cex = 1.05)
heading("(e)", "Time per resurgence", expression(paste(mu > 0, ",  pooled over runs")))

## shared legend for the time-series panels (a)-(c), centered in the outer margin above
## the top row: leave the layout for a full-device plot region without clearing the
## page, then place the legend at its top edge
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("top", inset = 0.006, horiz = TRUE, bty = "n", cex = 1.0 * cex_panel / par("cex"),
       seg.len = 2.4, x.intersp = 0.8, text.width = NA,   # column-wise widths (R >= 4.1), so the long last entry does not widen every column
       legend = c("defenders D", "public P", "disruptors A", "deterministic A (panel a)"),
       col = c(colD, colP, colA, colA), lwd = c(lwdD, lwdP, lwdA, lwdA), lty = c(1, 1, 1, 2))
dev.off()

## ==============================================================================
## Numbers quoted in the manuscript and the response letter
## ==============================================================================
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
rowN <- function(N) which(resA$N == N)
rowB <- function(N, m) which(resB$N == N & resB$mu == m)
pow10 <- function(N) sprintf("10^{%d}", round(log10(N)))
macros <- c(
  "% Generated by finite_population.R -- do not edit by hand",
  "\\newcommand{\\FPsimulator}{C++}",
  sprintf("\\newcommand{\\FPscale}{%s}", fmt(scale, 3)),
  sprintf("\\newcommand{\\FPrepsA}{%d}", reps_A),
  sprintf("\\newcommand{\\FPrepsB}{%d}", reps_B),
  sprintf("\\newcommand{\\FPTB}{%s}", fmt(T_B)),
  sprintf("\\newcommand{\\FPNAmax}{%s}", pow10(max(N_A))),
  sprintf("\\newcommand{\\FPNBmin}{%s}", pow10(min(N_B))),
  sprintf("\\newcommand{\\FPNBmax}{%s}", pow10(max(N_B))),
  sprintf("\\newcommand{\\FPcensoredA}{%d}", censored_A),
  sprintf("\\newcommand{\\FPsparseThreshold}{%d}", sparse_threshold),
  sprintf("\\newcommand{\\FPminEpisodes}{%d}", resB$episodes[i_min]),
  sprintf("\\newcommand{\\FPminEpisodesN}{%s}", pow10(resB$N[i_min])),
  sprintf("\\newcommand{\\FPminEpisodesMu}{10^{%d}}", round(log10(resB$mu[i_min]))),
  sprintf("\\newcommand{\\FPzeroEpisodeRuns}{%d}", zero_runs_min),
  sprintf("\\newcommand{\\FPpanelALoss}{%s}", fmt(panelA$t_ext * scale, 0)),
  sprintf("\\newcommand{\\FPpanelAFreeze}{%s}", fmt(t_freeze_a, 0)),
  sprintf("\\newcommand{\\FPepAtMostOneThousand}{%s}", fmt(100 * resA$ep_le1[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPepAtMostOneHundredThousand}{%s}", fmt(100 * resA$ep_le1[rowN(1e5)], 0)),
  sprintf("\\newcommand{\\FPepMedianThousand}{%s}", fmt(resA$ep_median[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPepMedianHundredThousand}{%s}", fmt(resA$ep_median[rowN(1e5)], 0)),
  sprintf("\\newcommand{\\FPepMeanThousand}{%s}", fmt(resA$ep_mean[rowN(1000)], 1)),
  sprintf("\\newcommand{\\FPepMeanHundredThousand}{%s}", fmt(resA$ep_mean[rowN(1e5)], 1)),
  sprintf("\\newcommand{\\FPextMedianHundred}{%s}", fmt(resA$ext_median[rowN(100)], 0)),
  sprintf("\\newcommand{\\FPextMedianThousand}{%s}", fmt(resA$ext_median[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPextMedianTenThousand}{%s}", fmt(resA$ext_median[rowN(10000)], 0)),
  sprintf("\\newcommand{\\FPextMedianHundredThousand}{%s}", fmt(resA$ext_median[rowN(1e5)], 0)),
  sprintf("\\newcommand{\\FPextSlope}{%s}", fmt(ext_slope, 1)),
  sprintf("\\newcommand{\\FPfreezeDThousand}{%s}", fmt(100 * resA$frac_D[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPfreezePThousand}{%s}", fmt(100 * resA$frac_P[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPfreezeAThousand}{%s}", fmt(100 * resA$frac_A[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPfreezeDHundredThousand}{%s}", fmt(100 * resA$frac_D[rowN(1e5)], 0)),
  sprintf("\\newcommand{\\FPlostPThousand}{%s}", fmt(100 * resA$lost_P[rowN(1000)], 0)),
  sprintf("\\newcommand{\\FPlostPHundredThousand}{%s}", fmt(100 * resA$lost_P[rowN(1e5)], 0)),
  sprintf("\\newcommand{\\FPperiodThousandEminusTwo}{%s}", fmt(resB$period[rowB(1000, 1e-2)], 0)),
  sprintf("\\newcommand{\\FPperiodThousandEminusThree}{%s}", fmt(resB$period[rowB(1000, 1e-3)], 0)),
  sprintf("\\newcommand{\\FPperiodThousandEminusFour}{%s}", fmt(resB$period[rowB(1000, 1e-4)], 0)),
  sprintf("\\newcommand{\\FPperiodThousandEminusFive}{%s}", fmt(resB$period[rowB(1000, 1e-5)], 0)),
  sprintf("\\newcommand{\\FPperiodTenThousandEminusFive}{%s}", fmt(resB$period[rowB(10000, 1e-5)], 0)),
  sprintf("\\newcommand{\\FPperiodTenThousandEminusTwo}{%s}", fmt(resB$period[rowB(10000, 1e-2)], 0)),
  sprintf("\\newcommand{\\FPfracHiThousandEminusThree}{%s}", fmt(100 * resB$frac_hi[rowB(1000, 1e-3)], 0)),
  sprintf("\\newcommand{\\FPfracHiMin}{%s}", fmt(100 * min(resB$frac_hi[resB$N >= 1000], na.rm = TRUE), 0)),
  sprintf("\\newcommand{\\FPfracHiMax}{%s}", fmt(100 * max(resB$frac_hi[resB$N >= 1000], na.rm = TRUE), 0))
)
writeLines(macros, file.path(out_dir, "results_finitepop.tex"))
write.csv(rbind(cbind(experiment = "A", N = resA$N, mu = 0, stat = "ext_median", value = resA$ext_median),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "ext_q25", value = resA$ext_q25),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "ext_q75", value = resA$ext_q75),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "frac_D", value = resA$frac_D),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "frac_P", value = resA$frac_P),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "frac_A", value = resA$frac_A),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "lost_first_P", value = resA$lost_P),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "not_absorbed", value = resA$censored),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "episodes_median", value = resA$ep_median),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "episodes_mean", value = resA$ep_mean),
                cbind(experiment = "A", N = resA$N, mu = 0, stat = "episodes_at_most_one", value = resA$ep_le1),
                cbind(experiment = "B", N = resB$N, mu = resB$mu, stat = "time_per_resurgence", value = resB$period),
                cbind(experiment = "B", N = resB$N, mu = resB$mu, stat = "episodes", value = resB$episodes),
                cbind(experiment = "B", N = resB$N, mu = resB$mu, stat = "frac_hi", value = resB$frac_hi)),
          file.path(out_dir, "finitepop_summary.csv"), row.names = FALSE)
write.csv(res_df[order(res_df$id), ], file.path(out_dir, "finitepop_runs.csv"), row.names = FALSE)
cat("\nWrote fig_finitepop.pdf, results_finitepop.tex, finitepop_summary.csv, finitepop_runs.csv to", out_dir, "\n")
cat(sprintf("Total time %.0f s\n", proc.time()[["elapsed"]] - t_start))
}, finally = .cleanup())
