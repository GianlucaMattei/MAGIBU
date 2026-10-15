#!/usr/bin/env Rscript
# Installs all required packages into the library of the bundled R. Run it with the Rscript of build/R-portable.
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "packages.R"))
lib <- .Library
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
install.packages(c("BiocManager", magibu_cran), lib = lib)
BiocManager::install(magibu_bioc, lib = lib, update = FALSE, ask = FALSE)
miss <- setdiff(c(magibu_cran, magibu_bioc), rownames(installed.packages(lib.loc = lib)))
if (length(miss) > 0) stop("Packages not installed: ", paste(miss, collapse = ", "))
cat("All packages installed in", lib, "\n")
