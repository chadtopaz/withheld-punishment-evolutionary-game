# Code and simulation outputs for *How Withheld Punishment Sustains Anti-Democratic Behavior in an Evolutionary Game*

Chad M. Topaz (Williams College; QSIDE Institute; University of Colorado Boulder)

This repository contains everything needed to reproduce the figures, the numerical checks of the two propositions, and every number quoted in the finite-population section of the article, together with the outputs those scripts produced on the author's machine. The repository is archived on Zenodo at `https://doi.org/10.5281/zenodo.22735002`, a DOI that represents all versions and resolves to the latest one; please cite that DOI and the article.

The model is a three-strategy replicator system on the simplex, with strategies *defenders* (D), *non-punishing public* (P), and *disruptors* (A). Two payoff regimes are analyzed, exploitation and accommodation, and the resurgence cycle of the exploitation regime is also simulated in finite populations with and without reintroduction of lost strategies.

## Contents

| File | Role |
| --- | --- |
| `verify_propositions.R` | Numerical checks of Proposition 1 (the exploitation regime) and Proposition 2 (the accommodation regime): the closed-form interior equilibrium, the trace and determinant identities, the exclusion threshold, and convergence from random initial conditions. Writes `verify_report.txt`. |
| `make_figures.R` | Draws Figure 2 (phase portraits of the exploitation regime, `fig_regimes.pdf`), writes `figures_log.txt`, and sources `make_figure3.R`. |
| `make_figure3.R` | Draws Figure 3 (the accommodation regime: regime diagram, equilibrium branches, and phase portraits, `fig_accommodation.pdf`) and writes `figure3_checks.txt`. Sourced by `make_figures.R`; not run on its own. |
| `finite_population.R` | The finite-population simulations of the article's Results (the subsection on strategy loss and resurgence in finite populations). Writes the C++ simulator `finitepop_sim.cpp`, compiles it with Rcpp, runs the two experiments in parallel, and writes Figure 4 (`fig_finitepop.pdf`), the LaTeX macros `results_finitepop.tex`, the tables `finitepop_summary.csv` and `finitepop_runs.csv`, and `finitepop_log.txt`. |
| `session_info.R` | Records the software environment (`session_info.txt`). |
| `fig_schematic.tex` | Standalone LaTeX/TikZ source of Figure 1 (the invasion relations of the two regimes); `pdflatex fig_schematic.tex` produces `fig_schematic.pdf`. In the article this figure is drawn inside the manuscript source with the same code. |
| `fig_regimes.pdf`, `fig_accommodation.pdf`, `fig_finitepop.pdf`, `fig_schematic.pdf` | Figures 2, 3, 4, and 1 as produced by the author. |
| `results_finitepop.tex` | LaTeX macros holding every number from the simulations that is quoted in the article (the finite-population results and the Figure 4 legend). |
| `finitepop_summary.csv` | One row per experimental cell (population size and reintroduction rate): first-loss statistics, which strategy was lost first, where the population froze, resurgence counts, pooled time per resurgence, and the fraction of time with a disruptor majority. |
| `finitepop_runs.csv` | One row per simulation run (7,320 rows), the raw results behind the summary. |
| `finitepop_sim.cpp` | The C++ simulator as written by `finite_population.R`; supplied for reference and regenerated on each run. |
| `verify_report.txt`, `figures_log.txt`, `figure3_checks.txt`, `finitepop_log.txt`, `session_info.txt` | Console output of the author's runs, including every check the scripts perform. |
| `checksums.txt` | MD5 checksums of the deterministic outputs, for comparison after a rerun. |

All scripts write their outputs next to themselves, so the repository root is the working directory. Rerunning the scripts regenerates the supplied outputs in place; `git diff` or `checksums.txt` then shows whether anything changed.

## Requirements

The author's runs used R 4.5.3 on macOS with Rcpp 1.1.0, future 1.70.0, and furrr 0.4.0; the exact platform and package versions are recorded in `session_info.txt`.

- **R** 4.1 or later. (`finite_population.R` uses a legend option introduced in R 4.1.)
- **CRAN packages** `Rcpp`, `future`, and `furrr`, needed only by `finite_population.R`, which installs any that are missing. `verify_propositions.R`, `make_figures.R`, and `make_figure3.R` use base R only.
- **A C++ compiler** for `finite_population.R`. On macOS, install the Xcode command line tools (`xcode-select --install`); on Linux, `g++` or `clang++`; on Windows, Rtools matching the R version. The script stops with a message if compilation fails; it never substitutes another simulator or a reduced design.
- **Cairo graphics support** in R (`capabilities("cairo")`) is used when available so that the figure fonts are embedded; otherwise the standard `pdf()` device is used and the figures are the same apart from font embedding.
- **LaTeX** with TikZ (any TeX Live from the last several years) only for compiling `fig_schematic.tex`.

Nothing needs to be installed beyond R, the three packages, and a compiler. There is no installation step for the repository itself: clone or unzip it and run the scripts from its root.

## How to run

Start R (or RStudio) with the repository root as the working directory and run the scripts in this order:

```r
source("verify_propositions.R")   # numerical checks of Propositions 1 and 2; base R; about one minute
source("make_figures.R")          # Figures 2 and 3; base R; under one minute
source("finite_population.R")     # Figure 4 and the finite-population numbers; Rcpp, future, furrr, C++ compiler
source("session_info.R")          # records the software environment
```

Then, if wanted, `pdflatex fig_schematic.tex` for Figure 1.

**Quick demonstration.** The first two scripts run in about a minute together and exercise the model, the integrators, and the figure code without the parallel simulation. `verify_propositions.R` prints its checks to the console and ends with `All checks passed.`; `make_figures.R` writes the two figure PDFs and reports the integrator checks for each panel. Compare the console output with the supplied `verify_report.txt`, `figures_log.txt`, and `figure3_checks.txt`.

**Runtime of the simulation.** `finite_population.R` runs 7,320 simulation runs (7,000 without reintroduction and 320 with it) on `parallel::detectCores() - 1` workers. In the author's run on 23 workers the simulations took 52 seconds and the whole script 57 seconds, with about 1,100 seconds of total worker time, so on a four-core machine expect roughly five to seven minutes. The log lists the wall time per run for each parameter set.

## Random seeds and reproducibility

All seeds are fixed inside the scripts. In `finite_population.R`, `set.seed(20260905)` precedes the parallel run, and `furrr` derives one L'Ecuyer-CMRG stream per bundle of runs from that state (`furrr_options(seed = TRUE)`). The bundles are fixed by the experimental design, so the results do not depend on the number of workers. The sampler checks use `set.seed(1)` to `set.seed(3)`, the single run shown in Figure 4a is chosen by trying seeds in order until one shows the most common outcome at its population size (seed 2 in the author's run, reported in the log), and the run shown in Figure 4b,c uses `set.seed(4)`.

Rerunning `finite_population.R` reproduces `results_finitepop.tex` and `finitepop_summary.csv` byte for byte (see `checksums.txt`). `finitepop_runs.csv` reproduces byte for byte except for its `elapsed` column, the wall time of each run. The figure PDFs reproduce apart from the creation date stored in the PDF metadata; the drawn content is identical.

## What the scripts check

Every script stops rather than writing output when a check fails.

- `verify_propositions.R`: over 20,000 random parameter draws in the accommodation regime, an interior equilibrium exists exactly when defenders can invade the coalition (H > 0); the closed-form equilibrium matches a direct linear solve to 1e-8; the trace and determinant identities hold to 1e-6; every interior equilibrium is a sink; the disruptor share at coexistence lies below its coalition value; the transverse exponent at the coalition equals H/(c_P + c_A) to 1e-10. Convergence from random initial conditions is checked in the hyperbolic cases (deviation below 1e-6) and, with the algebraic rate predicted by the center-manifold expansion, at the threshold H = 0. The exploitation-regime equilibrium, trace, and determinant of Figure 2a are printed and compared with their formulas.
- `make_figures.R` and `make_figure3.R`: the fixed-step RK4 integrators clip each share at 1e-300 and renormalize after each step; each script inspects the raw RK4 proposals and stops if any is negative or the coordinate sum drifts from one by 1e-9 or more. On the attracting cycle of Figure 2a the public share falls far below 1e-300 late in the run, so a nonzero clip count is expected there and is reported, not an error. `make_figure3.R` also verifies the signs, the equilibrium residuals, and a direct linear solve for the panels of Figure 3.
- `finite_population.R`: before the experiments, the event-driven sampler's pair decomposition is checked against the six transition probabilities computed from the definition at 307 states and four reintroduction rates (the script stops above 1e-12; the author's run reports 7e-17); a Monte Carlo check of the sampler's first transition and mean holding time at one interior state warns if |z| is large; and the C++ and R simulators are run from the same seed and must agree exactly on the sampled trajectory, the number of changes, and the number of episodes. After the runs, the script stops before writing any output unless every run without reintroduction ended with one strategy holding the whole population, because the article quotes first-loss times and fixation shares as summaries of all such runs; in that case the raw runs are saved as `finitepop_runs_INCOMPLETE.csv` for inspection.

## Expected results

From `verify_report.txt`: `All checks passed.`, with the maximum errors of the exact identities near machine precision.

From `finitepop_log.txt` (the numbers quoted in the article; population size N, reintroduction rate mu, time in replicator units):

- Without reintroduction, the median time to the first loss of a strategy is 8.7 at N = 100, 27.5 at N = 1,000, 48.1 at N = 10,000, and 66.4 at N = 100,000, growing like log N.
- The public strategy is lost first in 74% of runs at N = 1,000 and 89% at N = 100,000, and the population freezes at the all-defender state in 73% and 89% of all runs at those sizes.
- The median number of disruptor resurgences before fixation is 0, 1, 1, 1, 2, 2, 2 for N = 100, 300, 1,000, 3,000, 10,000, 30,000, 100,000.
- With reintroduction, resurgences occur in every cell. At N = 1,000 the pooled time per resurgence is 13.8, 24.4, 91.5, and 519.5 replicator units at mu = 1e-2, 1e-3, 1e-4, and 1e-5. The sparsest cell, N = 100 and mu = 1e-5, has 15 resurgences in 20 runs.
- One generation is 0.286 replicator time units (w = 1, largest payoff difference 3.5).

## How the outputs enter the article

`results_finitepop.tex` defines LaTeX macros (`\FPextMedianThousand`, `\FPfreezeDThousand`, and so on) that the manuscript reads with `\input{results_finitepop}`, so the numbers in the finite-population results and the Figure 4 legend are the ones the simulation wrote. Figures 2, 3, and 4 are included as the three figure PDFs. Figure 1 is drawn in the manuscript source with the TikZ code in `fig_schematic.tex`.

## License and citation

The code and the generated outputs in this repository are released under the MIT License (see `LICENSE`). To cite, use the Zenodo DOI above together with the article; `CITATION.cff` carries the metadata in machine-readable form.
