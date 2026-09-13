## session_info.R
## ------------------------------------------------------------------------------
## Records the software environment in which the scripts were run: R version,
## platform, operating system, attached and loaded packages with their versions,
## and the C++ compiler that Rcpp uses.  Run it after the three main scripts,
## from the folder that contains them:
##
##     source("session_info.R")
##
## It writes session_info.txt next to itself.  The copy supplied with the package
## is the author's.
## ------------------------------------------------------------------------------

out_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(out_dir) || length(out_dir) != 1 || is.na(out_dir) || !nzchar(out_dir)) out_dir <- getwd()

for (pkg in c("Rcpp", "future", "furrr")) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
suppressPackageStartupMessages({ library(Rcpp); library(future); library(furrr) })

con <- file(file.path(out_dir, "session_info.txt"), open = "wt")
sink(con)
cat("Recorded:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
cat("R.version.string:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
si <- Sys.info()
cat("Operating system:", si[["sysname"]], si[["release"]], "(", si[["version"]], ")\n")
cat("Machine:", si[["machine"]], "\n")
cat("Cores available (parallel::detectCores()):", parallel::detectCores(), "\n")
cat("cairo graphics support:", isTRUE(capabilities("cairo")), "\n\n")
cat("C++ compiler used by Rcpp (R CMD config CXX):\n")
cxx <- tryCatch(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE), error = function(e) "unavailable")
cat(" ", paste(cxx, collapse = " "), "\n")
ver <- tryCatch(system2(strsplit(cxx, " ")[[1]][1], "--version", stdout = TRUE, stderr = TRUE)[1], error = function(e) "unavailable")
cat(" ", ver, "\n\n")
cat("Package versions:\n")
for (pkg in c("Rcpp", "future", "furrr", "parallel")) cat(sprintf("  %-9s %s\n", pkg, as.character(packageVersion(pkg))))
cat("\n")
print(sessionInfo())
sink()
close(con)
cat("Wrote", file.path(out_dir, "session_info.txt"), "\n")
