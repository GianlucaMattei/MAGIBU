# Downloads the Windows binaries (R 4.4 / Bioconductor 3.20) of the required packages and unpacks them into the
# library of R-extract. Called by prepare_win_r.sh; paths are relative to this folder.
args <- commandArgs(trailingOnly = FALSE)
here <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "packages.R"))
options(timeout = 3600)
lib <- file.path(here, "R-extract", "app", "library")
zipdir <- file.path(here, "zips"); dir.create(zipdir, showWarnings = FALSE)
urls <- c(CRAN = "https://cran.r-project.org/bin/windows/contrib/4.4",
          BioC = "https://bioconductor.org/packages/3.20/bioc/bin/windows/contrib/4.4")
av <- do.call(rbind, lapply(names(urls), function(n) { a <- available.packages(contriburl = urls[[n]], type = "win.binary", filters = list()); cbind(a, Repo = n) }))
av <- av[!duplicated(av[, "Package"]), ]
cat("available:", nrow(av), "\n")
data_pkgs <- c("IlluminaHumanMethylationEPICv2anno.20a1.hg38", "IlluminaHumanMethylationEPICv2manifest",
               "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylation450kmanifest", "GenomeInfoDbData")
top <- setdiff(c(magibu_cran, magibu_bioc), data_pkgs)
base <- rownames(installed.packages(lib.loc = lib))
deps <- unique(c(top, unlist(tools::package_dependencies(top, db = av, which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
deps <- setdiff(deps, c(base, "R"))
need_data <- intersect(deps, data_pkgs)
missing <- setdiff(deps, c(rownames(av), data_pkgs))
cat("to download:", length(setdiff(deps, data_pkgs)), " | not found:", paste(missing, collapse = " "), "\n")
get <- setdiff(deps, c(data_pkgs, missing))
dl <- download.packages(get, destdir = zipdir, available = av, contriburl = urls, type = "win.binary", quiet = TRUE)
for (i in seq_len(nrow(dl))) utils::unzip(dl[i, 2], exdir = lib)
cat("unpacked:", nrow(dl), "\n")
writeLines(c(data_pkgs), file.path(zipdir, "data_pkgs.txt"))
