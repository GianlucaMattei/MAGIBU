#!/usr/bin/env Rscript
# Writes THIRD_PARTY_LICENSES.txt from the packages actually present in an R library.
# Usage: Rscript make_licenses.R [library] [output file]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else .Library
out <- if (length(args) >= 2) args[2] else "THIRD_PARTY_LICENSES.txt"
ip <- installed.packages(lib.loc = lib, fields = c("License", "URL"))
ip <- ip[order(tolower(ip[, "Package"])), , drop = FALSE]
txt <- c(
  "MAGIBU - THIRD PARTY SOFTWARE", "=============================", "",
  "Magibu is distributed together with the following third party software.",
  "Each component remains under its own licence; the full licence texts are in the",
  "R installation (R/share/licenses) and in each package folder (R/library/<package>).", "",
  "R - The R Project for Statistical Computing. Licence: GPL-2 | GPL-3.",
  "  Source code: https://cran.r-project.org/src/base/", "",
  "Electron. Licence: MIT. https://github.com/electron/electron (includes Chromium and Node.js, see LICENSES.chromium.html).", "",
  "magibu-meth (app/bin, BAM reader) is statically linked with htslib 1.21 (MIT/Expat), zlib 1.3.1 (zlib licence)",
  "  and, on Windows, winpthreads from mingw-w64 (MIT). Full texts: app/bin/LICENSES_magibu-meth.txt", "",
  "R packages (source code: https://cran.r-project.org and https://bioconductor.org):", "",
  sprintf("  %-48s %-12s %s", ip[, "Package"], ip[, "Version"], gsub("\\s+", " ", ip[, "License"])), "",
  "Reference model 'capper_v2': derived from the DNA methylation reference cohort of",
  "Capper et al., Nature 2018 (GEO GSE90496). Research use only.")
writeLines(txt, out)
cat("Written", out, "with", nrow(ip), "packages\n")
