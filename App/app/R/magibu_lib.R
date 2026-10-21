# ==============================================================================
# Magibu App - engine (used by the Shiny interface and from the command line)
# Wraps the functions in R/core (copies of Magibu/Scripts/FinalScripts) into three
# operations: build a model, predict, plot (class zoom and global map).
# ==============================================================================

MAGIBU_VERSION <- "0.7.1"

magibu_init <- function(app_dir) {
  suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(ggplot2)
  })
  core <- file.path(app_dir, "R", "core")
  for (f in list.files(core, pattern = "\\.[Rr]$", full.names = TRUE)) source(f, local = FALSE)
  options(magibu.app_dir = app_dir)
  invisible(TRUE)
}

magibu_resource <- function(f) file.path(getOption("magibu.app_dir"), "resources", f)

magibu_default_cores <- function() max(1L, min(4L, parallel::detectCores() - 1L, na.rm = TRUE))

# ------------------------------------------------------------------------------
# Model: one folder with <name>.umap_model.components.rds, <name>.umap_model.uwot,
# <name>.umap_metadata_result.tsv (class column already included)
# ------------------------------------------------------------------------------
magibu_model_files <- function(model_dir) {
  if (!nzchar(model_dir) || !dir.exists(model_dir)) stop("Model folder not found: ", model_dir)
  pick <- function(pattern, what) {
    f <- list.files(model_dir, pattern = pattern, full.names = TRUE)
    if (length(f) == 0) stop("Model file missing: ", what, " (", model_dir, ")")
    if (length(f) > 1) stop("More than one model in the folder: ", what, " (keep one model per folder)")
    f
  }
  list(components = pick("components\\.rds$", "*.umap_model.components.rds"),
       uwot       = pick("\\.uwot$", "*.umap_model.uwot"),
       metadata   = pick("metadata_result\\.tsv$", "*.umap_metadata_result.tsv"))
}

magibu_class_col <- function(dt) {
  for (col in c("sample.title", "methylation_class", "general_class", "diagnosis", "subclass")) {
    if (col %in% names(dt)) return(col)
  }
  NULL
}

# Metadata in standard form: Sample_ID_Full, UMAP1, UMAP2, sample.title, Batch
magibu_standard_meta <- function(meta) {
  meta <- data.table::as.data.table(meta)
  cc <- magibu_class_col(meta)
  if (is.null(cc)) meta[, sample.title := NA_character_]
  else if (cc != "sample.title") meta[, sample.title := as.character(get(cc))]
  meta[, sample.title := as.character(sample.title)]
  meta[!is.na(sample.title) & !nzchar(sample.title), sample.title := NA_character_]
  if (!"Batch" %in% names(meta)) meta[, Batch := NA_character_]
  meta[]
}

magibu_read_model_meta <- function(model_dir) {
  magibu_standard_meta(data.table::fread(magibu_model_files(model_dir)$metadata))
}

magibu_is_new <- function(meta) !is.na(meta$Batch) & meta$Batch == "New_Projected"

# ------------------------------------------------------------------------------
# Input: single files and/or folders -> list of samples (idat Grn/Red pairs, or tsv)
# ------------------------------------------------------------------------------
.mg_idat_re <- "_(Grn|Red)\\.idat(\\.gz)?$"
.mg_tsv_re  <- "\\.tsv(\\.gz)?$"
.mg_bam_re  <- "\\.bam$"

# What a tsv file contains, from its first line:
#   "magibu"    Magibu tsv (Probe_ID, Beta, chr_hg38, pos_hg38), any platform
#   "entropy"   magibu-meth or modkit extract + ModkitResorter (Poremeth2) output: chr, pos, -, -, beta (0-1), -
magibu_file_kind <- function(f) {
  l <- readLines(f, n = 5, warn = FALSE)
  l <- l[nzchar(l) & !startsWith(l, "#") & !startsWith(l, "track")]
  if (length(l) == 0) stop("Empty file: ", basename(f))
  if (grepl("Probe_ID", l[1], fixed = TRUE)) return("magibu")
  x <- strsplit(l[min(2, length(l))], "\t", fixed = TRUE)[[1]]
  num <- function(v) suppressWarnings(as.numeric(v))
  if (length(x) == 6 && !is.na(num(x[2])) && !is.na(num(x[5])) && num(x[5]) >= 0 && num(x[5]) <= 1) return("entropy")
  stop("Unrecognised methylation file (expected Magibu tsv or entropy tsv): ", basename(f))
}

.mg_sample_name <- function(f) sub("\\.(tsv|bam)$", "", sub("\\.gz$", "", basename(f)), ignore.case = TRUE)

magibu_collect_inputs <- function(paths) {
  paths <- unique(trimws(unlist(strsplit(paste(paths, collapse = "\n"), "\n", fixed = TRUE))))
  paths <- paths[nzchar(paths)]
  if (length(paths) == 0) stop("No file or folder given.")
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) stop("Path not found: ", paste(missing, collapse = "; "))
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, pattern = paste(.mg_idat_re, .mg_tsv_re, .mg_bam_re, sep = "|"), full.names = TRUE, ignore.case = TRUE)
    else p
  }))
  is_idat <- grepl(.mg_idat_re, files, ignore.case = TRUE)
  is_meth <- grepl(.mg_tsv_re, files, ignore.case = TRUE)
  is_bam <- grepl(.mg_bam_re, files, ignore.case = TRUE)
  bad <- files[!is_idat & !is_meth & !is_bam]
  if (length(bad) > 0) stop("Unrecognised file (expected *_Grn.idat / *_Red.idat, also .gz, nanopore BAM or .tsv): ",
                            paste(basename(bad), collapse = "; "))
  idat_base <- unique(sub(.mg_idat_re, "", files[is_idat], ignore.case = TRUE))
  meth <- files[is_meth]
  # Magibu tsv are used as they are; entropy tsv are converted first
  meth_kind <- vapply(meth, function(f) if (magibu_file_kind(f) == "magibu") "tsv" else "nanopore", "")
  out <- rbind(
    data.frame(sample = basename(idat_base), kind = rep("idat", length(idat_base)), path = idat_base, stringsAsFactors = FALSE),
    data.frame(sample = vapply(meth, .mg_sample_name, ""), kind = unname(meth_kind), path = meth, stringsAsFactors = FALSE),
    data.frame(sample = vapply(files[is_bam], .mg_sample_name, ""), kind = rep("bam", sum(is_bam)), path = files[is_bam], stringsAsFactors = FALSE))
  rownames(out) <- NULL
  if (nrow(out) == 0) stop("No idat, BAM or tsv file found in the given paths.")
  dup <- unique(out$sample[duplicated(out$sample)])
  if (length(dup) > 0) stop("Same sample given more than once: ", paste(dup, collapse = "; "))
  out
}

# Platform from the idat (number of probes read): 450k ~622k, EPICv1 ~1.05M, EPICv2 ~1.1M
magibu_detect_array <- function(grn_idat) {
  n <- illuminaio::readIDAT(grn_idat)$nSNPsRead
  if (n < 800000) "Infinium" else if (n < 1080000) "EpicV1" else "EpicV2"
}

.mg_copy_idat <- function(base, channel, dest_dir) {
  src <- paste0(base, "_", channel, ".idat")
  dest <- file.path(dest_dir, basename(src))
  if (file.exists(src)) {
    if (!file.copy(src, dest, overwrite = TRUE)) stop("Copy failed: ", src)
  } else if (file.exists(paste0(src, ".gz"))) {
    con <- gzfile(paste0(src, ".gz"), "rb"); on.exit(close(con))
    out <- file(dest, "wb")
    repeat { buf <- readBin(con, "raw", 4e6); if (length(buf) == 0) break; writeBin(buf, out) }
    close(out)
  } else {
    stop("Missing file (both channels _Grn and _Red are needed): ", basename(src))
  }
  dest
}

# idat -> Magibu tsv (Probe_ID, Beta, chr_hg38, pos_hg38). Returns the tsv paths.
magibu_idat_to_tsv <- function(inputs, work_dir, tsv_dir, array_type = "auto", n_cores = magibu_default_cores()) {
  idat <- inputs[inputs$kind == "idat", , drop = FALSE]
  if (nrow(idat) == 0) return(character(0))
  dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
  stage <- file.path(work_dir, "idat_staging")
  unlink(stage, recursive = TRUE); dir.create(stage, recursive = TRUE)
  message("Preparing ", nrow(idat), " idat samples...")
  types <- character(nrow(idat))
  for (i in seq_len(nrow(idat))) {
    tmp <- file.path(stage, "_in"); dir.create(tmp, showWarnings = FALSE)
    grn <- .mg_copy_idat(idat$path[i], "Grn", tmp)
    red <- .mg_copy_idat(idat$path[i], "Red", tmp)
    types[i] <- if (identical(array_type, "auto")) magibu_detect_array(grn) else array_type
    d <- file.path(stage, types[i]); dir.create(d, showWarnings = FALSE)
    file.rename(c(grn, red), file.path(d, basename(c(grn, red))))
  }
  unlink(file.path(stage, "_in"), recursive = TRUE)
  message("Platforms: ", paste(names(table(types)), table(types), sep = " x", collapse = ", "))
  for (ty in unique(types)) {
    # one tsv folder per platform: the EpicV1 branch re-annotates EVERY tsv it finds
    tsv_ty <- file.path(stage, paste0("tsv_", ty))
    annot <- switch(ty, EpicV1 = magibu_resource("epic_v1_manifest_hg38.csv"),
                        Infinium = magibu_resource("annot_450k_hg38_from_epicv2.tsv"), NULL)
    process_methylation_data(path.to.idat = file.path(stage, ty), save.path = tsv_ty,
                             array.type = ty, annot.file.path = annot, n.cores = n_cores)
    made <- list.files(tsv_ty, pattern = "\\.tsv$", full.names = TRUE)
    file.copy(made, tsv_dir, overwrite = TRUE)
  }
  unlink(stage, recursive = TRUE)
  out <- file.path(tsv_dir, paste0(idat$sample, ".tsv"))
  if (any(!file.exists(out))) stop("idat conversion failed for: ", paste(idat$sample[!file.exists(out)], collapse = "; "))
  out
}

# ------------------------------------------------------------------------------
# Nanopore entropy tsv -> Magibu tsv. Keeps only CpGs within `window` bp of a model CpG:
# the projection uses nothing else (exact match, or nearest CpG within 500 bp), so the
# result is unchanged while the file shrinks from ~30M to ~1M rows.
# Entropy tsv (magibu-meth from a BAM, or modkit extract + ModkitResorter): beta = column 5,
# as in the original Magibu nanopore workflow.
# ------------------------------------------------------------------------------
.mg_model_positions <- function(mf) {
  p <- readRDS(mf$components)$common_positions
  p <- unique(as.character(p[!is.na(p)]))
  data.table::data.table(chr = sub(":.*", "", p), pos = as.integer(sub(".*:", "", p)))
}

.mg_near_model <- function(dt, model, window) {
  keep <- logical(nrow(dt))
  for (ch in intersect(unique(dt$chr), unique(model$chr))) {
    ref <- sort(model$pos[model$chr == ch]); i <- which(dt$chr == ch); x <- dt$pos[i]
    j <- findInterval(x, ref)
    d_lo <- ifelse(j >= 1, x - ref[pmax(j, 1)], Inf)
    d_hi <- ifelse(j < length(ref), ref[pmin(j + 1, length(ref))] - x, Inf)
    keep[i] <- pmin(abs(d_lo), abs(d_hi)) <= window
  }
  dt[keep]
}

magibu_read_nanopore <- function(f) {
  if (magibu_file_kind(f) != "entropy") stop("Not a nanopore file: ", basename(f))
  dt <- data.table::fread(f, sep = "\t", header = FALSE, select = c(1, 2, 5), col.names = c("chr", "pos", "beta"))
  dt[, chr := as.character(chr)]
  dt[!startsWith(chr, "chr"), chr := paste0("chr", chr)]
  dt[!is.na(beta) & beta >= 0 & beta <= 1]
}

# Share of model CpGs measured in a Magibu tsv: exact position (+-1 bp) or a CpG within `window` bp,
# i.e. what the projection can use instead of the training mean. Low values = unreliable position.
magibu_cpg_coverage <- function(tsv, model, window = 502L) {
  dt <- data.table::fread(tsv, select = c("chr_hg38", "pos_hg38"), col.names = c("chr", "pos"))
  dt[, chr := as.character(chr)]
  hit <- .mg_near_model(model, dt, window)      # model CpGs with a sample CpG nearby
  round(100 * nrow(hit) / nrow(model), 1)
}

# ------------------------------------------------------------------------------
# Nanopore BAM (dorado, MM/ML methylation tags) -> beta per CpG with magibu-meth --entropy --no-entropy
# (C, htslib, build/native/prog/magibu_meth.c). Beta as in the original Magibu nanopore workflow
# (modkit extract full + ModkitResorter.sh / ParseModkit.pl of PoreMeth2: 5mC if p > 0.85, unmethylated
# if p5mC + p5hmC < 0.3; with include_5hmc methylated if p5mC + p5hmC > 0.85), without the entropy,
# so no 3-CpG window and no 3-read pattern rule: every CpG with at least min_cov confident calls. Through the BAM index only the regions around the model CpGs
# (+-flank bp) are read, so the regions follow the chosen model. No reference FASTA needed.
# The file (entropy tsv layout, entropy = NA) is kept in <out>/Bam_beta: only the regions of this model.
# ------------------------------------------------------------------------------
magibu_meth_bin <- function() {
  sys <- Sys.info()[["sysname"]]
  f <- switch(sys, Windows = "win-x64/magibu-meth.exe", Linux = "linux-x64/magibu-meth", "")
  p <- if (nzchar(f)) file.path(getOption("magibu.app_dir"), "bin", f) else ""
  if (!nzchar(p) || !file.exists(p))
    stop("Reading BAM files is not available on this system (Windows and Linux only): ", sys)
  p
}

magibu_bam_to_beta <- function(inputs, out_dir, mf, n_cores = magibu_default_cores(), min_cov = 1,
                               include_5hmc = FALSE, flank = 600L) {
  is_bam <- inputs$kind == "bam"
  if (!any(is_bam)) return(inputs)
  bin <- magibu_meth_bin()
  for (f in inputs$path[is_bam]) {
    idx <- c(paste0(f, ".bai"), sub("\\.bam$", ".bai", f, ignore.case = TRUE), paste0(f, ".csi"))
    if (!any(file.exists(idx))) stop("BAM index (.bai) not found, create it with samtools index: ", basename(f))
  }
  ent_dir <- file.path(out_dir, "Bam_beta")
  dir.create(ent_dir, recursive = TRUE, showWarnings = FALSE)
  regions <- file.path(ent_dir, "model_regions.bed")
  model <- .mg_model_positions(mf)
  data.table::fwrite(model[order(chr, pos), .(chr, start = pmax(0L, pos - flank), end = pos + flank)],
                     regions, sep = "\t", col.names = FALSE)
  for (i in which(is_bam)) {
    f <- inputs$path[i]
    out <- file.path(ent_dir, paste0(inputs$sample[i], ".beta.tsv"))
    message("BAM ", sum(is_bam[seq_len(i)]), "/", sum(is_bam), ": ", basename(f), " (only the regions of the model are read, ",
            if (include_5hmc) "5mC + 5hmC" else "5mC", ")")
    last <- -1
    r <- processx::run(bin, c("--entropy", "--no-entropy", "--min-cov", as.character(max(1L, as.integer(min_cov))),
                                if (include_5hmc) "--sum-5hmc", "-r", regions, "-t", as.character(n_cores), "-o", out, f), error_on_status = FALSE,
                       stderr_line_callback = function(l, proc) {
                         if (startsWith(l, "PROGRESS ")) {
                           x <- as.numeric(strsplit(substring(l, 10), "/", fixed = TRUE)[[1]])
                           pct <- 10 * floor(10 * x[1] / x[2])
                           if (pct > last) { last <<- pct; message("  ", pct, "%") }
                         } else if (nzchar(l)) message("  ", l)
                       })
    if (r$status != 0 || !file.exists(out)) {
      err <- grep("^ERROR", strsplit(r$stderr, "\n", fixed = TRUE)[[1]], value = TRUE)
      stop("BAM reading failed: ", basename(f), if (length(err)) paste0(" (", sub("^ERROR: ", "", err[1]), ")"))
    }
    if (file.info(out)$size == 0) stop("BAM reading failed: ", basename(f), " (no CpG with enough reads in the model regions)")
    inputs$kind[i] <- "nanopore"; inputs$path[i] <- out
  }
  inputs
}

magibu_nanopore_to_tsv <- function(inputs, tsv_dir, mf, window = 502L) {
  nano <- inputs[inputs$kind == "nanopore", , drop = FALSE]
  if (nrow(nano) == 0) return(character(0))
  dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
  model <- .mg_model_positions(mf)
  out <- file.path(tsv_dir, paste0(nano$sample, ".tsv"))
  for (i in seq_len(nrow(nano))) {
    message("Nanopore ", i, "/", nrow(nano), ": ", basename(nano$path[i]))
    dt <- magibu_read_nanopore(nano$path[i])
    n_all <- nrow(dt)
    dt <- .mg_near_model(dt, model, window)
    message("  CpGs read: ", n_all, ", near the model CpGs: ", nrow(dt))
    if (nrow(dt) < 0.05 * nrow(model))
      warning("Few CpGs near the model positions in ", basename(nano$path[i]), " (", nrow(dt), "): check that coordinates are hg38.", call. = FALSE)
    if (nrow(dt) == 0) stop("No CpG near the model positions: ", basename(nano$path[i]))
    data.table::fwrite(dt[, .(Probe_ID = paste0(chr, ":", pos), Beta = round(beta, 4), chr_hg38 = chr, pos_hg38 = pos)],
                       out[i], sep = "\t")
  }
  out
}

# ------------------------------------------------------------------------------
# Sample -> reference distances
# ------------------------------------------------------------------------------
magibu_distances <- function(meta, sample_ids, knn = 5, k_majority = c(5, 7, 9, 11)) {
  meta <- magibu_standard_meta(meta)
  ref <- meta[!magibu_is_new(meta) & !is.na(sample.title)]
  if (nrow(ref) == 0) stop("The model has no annotated classes: distances cannot be computed.")
  long <- list(); summ <- list()
  for (s in sample_ids) {
    t <- meta[Sample_ID_Full == s][1]
    d <- sqrt((ref$UMAP1 - t$UMAP1)^2 + (ref$UMAP2 - t$UMAP2)^2)
    dd <- data.table::data.table(Class = ref$sample.title, Ref_ID = ref$Sample_ID_Full, Dist = d)[order(Dist)]
    per <- dd[, .(MeanDist_kNN = mean(Dist[seq_len(min(.N, knn))]), MinDist = Dist[1], N_ref = .N), by = Class][order(MeanDist_kNN)]
    per[, `:=`(Sample = s, Rank = seq_len(.N))]
    long[[s]] <- per[, .(Sample, Rank, Class, MeanDist_kNN = round(MeanDist_kNN, 4), MinDist = round(MinDist, 4), N_ref)]
    row <- data.table::data.table(Sample = s, UMAP1 = round(t$UMAP1, 4), UMAP2 = round(t$UMAP2, 4))
    for (j in 1:3) {
      row[[paste0("NN", j, "_Class")]] <- dd$Class[j]; row[[paste0("NN", j, "_Dist")]] <- round(dd$Dist[j], 4)
      row[[paste0("NN", j, "_ID")]] <- dd$Ref_ID[j]
    }
    for (j in 1:3) {
      row[[paste0("Group", j, "_Class")]] <- per$Class[j]; row[[paste0("Group", j, "_MeanDist")]] <- round(per$MeanDist_kNN[j], 4)
    }
    for (k in k_majority) {
      tb <- sort(table(dd$Class[seq_len(min(k, nrow(dd)))]), decreasing = TRUE)
      row[[paste0("kNN", k, "_Class")]] <- names(tb)[1]; row[[paste0("kNN", k, "_Purity")]] <- round(tb[[1]] / sum(tb), 2)
    }
    summ[[s]] <- row
  }
  long <- data.table::rbindlist(long); summ <- data.table::rbindlist(summ)
  wide <- data.table::dcast(long, Class ~ Sample, value.var = "MeanDist_kNN")
  data.table::setnames(long, "MeanDist_kNN", paste0("MeanDist_kNN", knn))
  list(Summary = summ, Distances = long, MeanDist_wide = wide)
}

# ------------------------------------------------------------------------------
# Sample names (display names chosen by the user; the ID stays the file name)
# ------------------------------------------------------------------------------
# Named vector id -> display name (the id itself when no name was given)
magibu_sample_names <- function(res) {
  s <- res$samples
  out <- stats::setNames(s, s)
  l <- res$labels
  if (!is.null(l)) {
    l <- l[names(l) %in% s & !is.na(l) & nzchar(trimws(l))]
    out[names(l)] <- trimws(l)
  }
  out
}

# Names must be unique: a repeated name gets its id appended
magibu_clean_labels <- function(labels, samples) {
  l <- stats::setNames(rep("", length(samples)), samples)
  if (!is.null(labels)) {
    keep <- intersect(names(labels), samples)
    l[keep] <- trimws(as.character(labels[keep]))
  }
  l[is.na(l)] <- ""
  shown <- ifelse(nzchar(l), l, samples)
  dup <- shown %in% shown[duplicated(shown)] & nzchar(l)
  l[dup] <- paste0(l[dup], " (", samples[dup], ")")
  l
}

.mg_with_names <- function(dt, nm) {
  dt <- data.table::copy(dt)
  dt[, Name := unname(nm[Sample])]
  data.table::setcolorder(dt, c("Sample", "Name"))
  dt[]
}

# Classic plot (local neighbourhood + bar plot): one page block per sample, header with its name
magibu_classic_pdf <- function(meta, samples, names, file, n_top_groups = 10, ncol = 2) {
  meta <- as.data.frame(meta)
  plots <- lapply(seq_along(samples), function(k) {
    sid <- samples[k]
    message("Plot: ", sid)
    p <- suppressWarnings(plot_umap_zoom(meta, sid, cases_files = samples, circle_radius = 0.3,
                                         distance_multiplier = 1, n_top_groups = n_top_groups))
    head <- if (names[[sid]] == sid) sid else paste0(names[[sid]], "  (", sid, ")")
    if (is.null(p)) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", 0, 0, label = paste("No data:", head)))
    hdr <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = head, size = 6, fontface = "bold") +
      ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(b = 5, t = 10))
    patchwork::wrap_plots(hdr, p, ncol = 1, heights = c(1, 15))
  })
  height <- max(9, 5 * length(samples) / 2)
  width <- max(16, min(ncol, length(samples)) * 10)
  ggplot2::ggsave(file, patchwork::wrap_plots(plots, ncol = ncol), width = width, height = height, limitsize = FALSE)
  file
}

# Writes Excel + classic PDF from a result (with the current names)
magibu_write_outputs <- function(res, pdf = TRUE) {
  nm <- magibu_sample_names(res)
  d <- res$distances
  wide <- data.table::copy(d$MeanDist_wide)
  data.table::setnames(wide, res$samples, unname(nm[res$samples]), skip_absent = TRUE)
  writexl::write_xlsx(list(Summary = .mg_with_names(d$Summary, nm), Distances = .mg_with_names(d$Distances, nm),
                           MeanDist_wide = wide), res$files$xlsx)
  if (pdf) {
    message("Classic plot (UMAP + bar plot)...")
    magibu_classic_pdf(res$metadata, res$samples, nm, res$files$pdf)
  }
  invisible(res)
}

magibu_result_rds <- function(res) file.path(res$out_dir, paste0(res$name, ".magibu_result.rds"))

# Saves new names: updates Excel, classic PDF and the result file
magibu_apply_names <- function(res, labels) {
  res$labels <- magibu_clean_labels(labels, res$samples)
  magibu_write_outputs(res)
  saveRDS(res, magibu_result_rds(res))
  res
}

# ------------------------------------------------------------------------------
# PREDICTION
# ------------------------------------------------------------------------------
# The UMAP projection (uwot::umap_transform) is stochastic: without a seed the same sample moves
# by ~0.2 units between runs. seed = NULL (default) keeps the original Magibu behaviour (one batch,
# no seed). With a number every sample is projected ALONE with that seed: same file -> same point.
magibu_project <- function(tsv_files, meta_train, mf, save_dir, name, seed = NULL) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  comp <- readRDS(mf$components)
  uw <- uwot::load_uwot(mf$uwot)
  tmp <- file.path(save_dir, "_tmp")
  if (is.null(seed)) {
    meta <- project_new_samples(meta_train, comp, uw, new.samples = tsv_files, save.path = tmp, name = name)$metadata
  } else {
    new_rows <- lapply(seq_along(tsv_files), function(i) {
      message("--- Sample ", i, "/", length(tsv_files), ": ", basename(tsv_files[i]))
      set.seed(seed)
      m <- data.table::as.data.table(project_new_samples(meta_train, comp, uw, new.samples = tsv_files[i],
                                                         save.path = tmp, name = name)$metadata)
      m[!is.na(Batch) & Batch == "New_Projected"]
    })
    meta <- data.table::rbindlist(c(list(meta_train), new_rows), use.names = TRUE, fill = TRUE)
  }
  unlink(tmp, recursive = TRUE)
  meta <- magibu_standard_meta(meta)
  data.table::fwrite(meta, file.path(save_dir, paste0(name, "_metadata_combined.tsv")), sep = "\t")
  meta
}

magibu_run_prediction <- function(paths, model_dir, out_dir, name = "magibu", array_type = "auto",
                                  n_cores = magibu_default_cores(), knn = 5, seed = NULL, labels = NULL,
                                  min_cov = 1, include_5hmc = FALSE) {
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)
  mf <- magibu_model_files(model_dir)
  inputs <- magibu_collect_inputs(paths)
  message("Samples: ", nrow(inputs), " (", sum(inputs$kind == "idat"), " idat, ", sum(inputs$kind == "bam"), " BAM, ",
          sum(inputs$kind == "nanopore"), " entropy tsv, ", sum(inputs$kind == "tsv"), " tsv)")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tsv_dir <- file.path(out_dir, "Tsv")
  tsv_new <- magibu_idat_to_tsv(inputs, out_dir, tsv_dir, array_type, n_cores)
  inputs <- magibu_bam_to_beta(inputs, out_dir, mf, n_cores, min_cov, include_5hmc)
  tsv_nano <- magibu_nanopore_to_tsv(inputs, tsv_dir, mf)
  tsv_all <- c(tsv_new, tsv_nano, inputs$path[inputs$kind == "tsv"])

  meta_train <- magibu_read_model_meta(model_dir)
  meta <- magibu_project(tsv_all, meta_train, mf, file.path(out_dir, "Projection"), name, seed)
  samples <- meta$Sample_ID_Full[magibu_is_new(meta)]

  model <- .mg_model_positions(mf)
  cover <- vapply(tsv_all, magibu_cpg_coverage, 0, model = model)
  names(cover) <- tools::file_path_sans_ext(basename(tsv_all))
  for (k in names(cover)) message("Model CpGs measured in ", k, ": ", cover[[k]], "%",
                                  if (cover[[k]] < 50) "  (LOW: the position on the map is unreliable)" else "")

  message("Computing distances...")
  dist <- magibu_distances(meta, samples, knn = knn)
  dist$Summary[, Model_CpG_coverage_pct := unname(cover[Sample])]
  data.table::setcolorder(dist$Summary, c("Sample", "Model_CpG_coverage_pct"))
  res <- list(version = MAGIBU_VERSION, name = name, model_dir = model_dir, out_dir = out_dir, date = Sys.time(),
              seed = seed, metadata = meta, samples = samples, distances = dist,
              labels = magibu_clean_labels(labels, samples),
              files = list(xlsx = file.path(out_dir, paste0(name, ".distances.xlsx")),
                           pdf = file.path(out_dir, paste0(name, ".results.plots.pdf"))))
  magibu_write_outputs(res)
  rds <- magibu_result_rds(res)
  saveRDS(res, rds)
  message("=== PREDICTION DONE === ", out_dir)
  rds
}

magibu_find_result <- function(path) {
  if (file.exists(path) && !dir.exists(path)) return(path)
  f <- list.files(path, pattern = "\\.magibu_result\\.rds$", full.names = TRUE)
  if (length(f) == 0) stop("No Magibu result (*.magibu_result.rds) in: ", path)
  f[which.max(file.mtime(f))]
}

# ------------------------------------------------------------------------------
# BUILD A MODEL
# ------------------------------------------------------------------------------
magibu_read_annotation <- function(annotation_file, id_col = NULL, class_col = NULL) {
  if (!file.exists(annotation_file)) stop("Annotation file not found: ", annotation_file)
  ann <- if (grepl("\\.xlsx?$", annotation_file, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Reading xlsx needs the readxl package: save the annotation as tsv/csv.")
    data.table::as.data.table(readxl::read_excel(annotation_file))
  } else data.table::fread(annotation_file)
  if (ncol(ann) < 2) stop("The annotation needs at least two columns: sample identifier and class.")
  if (is.null(id_col) || !nzchar(id_col)) id_col <- names(ann)[1]
  if (is.null(class_col) || !nzchar(class_col)) class_col <- names(ann)[2]
  if (!all(c(id_col, class_col) %in% names(ann))) stop("Columns not found in the annotation: ", id_col, ", ", class_col)
  out <- data.table::data.table(id = as.character(ann[[id_col]]), class = as.character(ann[[class_col]]))
  out <- out[!is.na(id) & nzchar(id)]
  out[!duplicated(id)]
}

# exact id, otherwise the prefix before the first "_" (e.g. GSM123_sentrix)
magibu_match_annotation <- function(sample_ids, ann) {
  cls <- ann$class[match(sample_ids, ann$id)]
  miss <- is.na(cls)
  cls[miss] <- ann$class[match(sub("_.*", "", sample_ids[miss]), ann$id)]
  cls
}

magibu_slim_components <- function(rds) {
  m <- readRDS(rds)
  k <- if (!is.null(m$pca_components)) m$pca_components else 50L
  p <- m$pca_model
  k <- min(k, ncol(p$rotation))
  slim <- list(sdev = p$sdev[seq_len(k)], rotation = p$rotation[, seq_len(k), drop = FALSE],
               center = p$center, scale = p$scale)
  class(slim) <- "prcomp"
  m$pca_model <- slim
  saveRDS(m, rds)
  invisible(rds)
}

magibu_reference_map <- function(meta, title = "Reference map") {
  meta <- magibu_standard_meta(meta)
  ref <- meta[!magibu_is_new(meta)]
  n <- data.table::uniqueN(ref$sample.title)
  ggplot2::ggplot(ref, ggplot2::aes(UMAP1, UMAP2, color = sample.title)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::coord_fixed() + ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(0.35, "cm")) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = max(1, ceiling(n / 35)), override.aes = list(size = 2.5, alpha = 1))) +
    ggplot2::labs(title = title, subtitle = paste0(nrow(ref), " samples, ", n, " classes"),
                  color = "Class", x = "UMAP 1", y = "UMAP 2")
}

magibu_run_training <- function(dataset_dir, annotation_file, model_dir, name = "model",
                                id_col = NULL, class_col = NULL, array_type = "auto",
                                n_top_vars = 50000, pca_components = 50, n_neighbors = 30, min_dist = 0.3,
                                n_cores = magibu_default_cores()) {
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)
  if (!dir.exists(dataset_dir)) stop("Dataset folder not found: ", dataset_dir)
  if (dir.exists(model_dir) && length(list.files(model_dir, pattern = "components\\.rds$|\\.uwot$")) > 0)
    stop("The destination folder already contains a model, choose a new or empty one: ", model_dir)
  inputs <- magibu_collect_inputs(dataset_dir)
  if (any(inputs$kind %in% c("nanopore", "bam")))
    stop("Building a model from modkit files is not supported, convert them to Magibu tsv first: ",
         paste(head(basename(inputs$path[inputs$kind %in% c("nanopore", "bam")]), 3), collapse = "; "))
  n <- nrow(inputs)
  message("Dataset: ", n, " samples (", sum(inputs$kind == "idat"), " idat, ", sum(inputs$kind == "tsv"), " tsv)")

  # annotation checked BEFORE the long steps
  ann <- magibu_read_annotation(annotation_file, id_col, class_col)
  cls <- magibu_match_annotation(inputs$sample, ann)
  n_ok <- sum(!is.na(cls))
  message("Annotation: ", n_ok, "/", n, " samples with a class, ", data.table::uniqueN(cls[!is.na(cls)]), " classes")
  if (n_ok == 0) stop("No dataset sample matches the annotation identifiers. ",
                      "Example sample: '", inputs$sample[1], "', example annotation id: '", ann$id[1], "'")
  if (n < 15) stop("At least 15 samples are needed to build a model (found ", n, ").")
  pca_components <- as.integer(min(pca_components, n - 1))
  n_neighbors <- as.integer(min(n_neighbors, n - 1))

  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  work <- file.path(model_dir, "_work")
  if (all(inputs$kind == "tsv")) {
    tsv_dir <- dataset_dir
  } else {
    tsv_dir <- file.path(work, "Tsv")
    magibu_idat_to_tsv(inputs, work, tsv_dir, array_type, n_cores)
    file.copy(inputs$path[inputs$kind == "tsv"], tsv_dir, overwrite = TRUE)
  }

  generate_umap_training_model(tsv.path = tsv_dir, save.path = model_dir, metadata_file = NULL, name = name,
                               n_top_vars = n_top_vars, pca_components = pca_components,
                               n_neighbors = n_neighbors, min_dist = min_dist, n_cores = n_cores)
  mf <- magibu_model_files(model_dir)
  meta <- data.table::fread(mf$metadata)
  meta[, sample.title := magibu_match_annotation(Sample_ID_Full, ann)]
  data.table::fwrite(meta, mf$metadata, sep = "\t")
  message("Compacting the model...")
  magibu_slim_components(mf$components)

  ggplot2::ggsave(file.path(model_dir, paste0(name, ".reference_map.pdf")),
                  magibu_reference_map(meta, paste0("Reference map: ", name)), width = 16, height = 11, limitsize = FALSE)
  writeLines(c(paste("Magibu model:", name), paste("Created:", format(Sys.time())), paste("Magibu App:", MAGIBU_VERSION),
               paste("Dataset:", dataset_dir), paste("Annotation:", annotation_file),
               paste("Samples:", nrow(meta)), paste("Annotated samples:", sum(!is.na(meta$sample.title))),
               paste("Classes:", data.table::uniqueN(meta$sample.title[!is.na(meta$sample.title)])),
               paste("n_top_vars:", n_top_vars), paste("pca_components:", pca_components),
               paste("n_neighbors:", n_neighbors), paste("min_dist:", min_dist)),
             file.path(model_dir, paste0(name, ".model_info.txt")))
  message("=== MODEL BUILT === ", model_dir)
  model_dir
}

# ------------------------------------------------------------------------------
# Colours: one per class (or family). Hues spread so that alphabetical neighbours
# (often neighbours on the map too) never get similar colours.
# ------------------------------------------------------------------------------
magibu_family <- function(x) trimws(sub(",.*", "", x))

magibu_class_colors <- function(classes) {
  classes <- sort(unique(classes[!is.na(classes)]))
  n <- length(classes)
  if (n == 0) return(character(0))
  base <- grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = c(80, 55, 95)[(seq_len(n) - 1) %% 3 + 1],
                         l = c(55, 72, 40)[(seq_len(n) - 1) %% 3 + 1])
  gcd <- function(a, b) { while (b > 0) { t <- b; b <- a %% b; a <- t }; a }
  step <- max(1L, as.integer(round(n * 0.38)))
  while (n > 2 && gcd(step, n) != 1) step <- step + 1L
  stats::setNames(base[((seq_len(n) - 1) * step) %% n + 1], classes)
}

# ------------------------------------------------------------------------------
# GLOBAL MAP: every reference, samples as small asterisks. xlim/ylim = current view.
# ------------------------------------------------------------------------------
magibu_global_data <- function(meta, color_by = "class") {
  meta <- magibu_standard_meta(meta)
  ref <- meta[!magibu_is_new(meta)]
  ref[, group := if (color_by == "family") magibu_family(sample.title) else sample.title]
  ref[is.na(group), group := "Unlabelled"]
  ref
}

magibu_global_plot <- function(meta, samples = character(0), names = NULL, color_by = "class",
                               xlim = NULL, ylim = NULL, label_samples = TRUE, legend = TRUE, title = NULL) {
  ref <- magibu_global_data(meta, color_by)
  meta <- magibu_standard_meta(meta)
  fg <- meta[magibu_is_new(meta) & Sample_ID_Full %in% samples]
  pal <- magibu_class_colors(ref$group)
  if ("Unlabelled" %in% ref$group) pal["Unlabelled"] <- "grey80"
  zoomed <- !is.null(xlim)
  if (zoomed) ref <- ref[UMAP1 >= xlim[1] & UMAP1 <= xlim[2] & UMAP2 >= ylim[1] & UMAP2 <= ylim[2]]
  n <- data.table::uniqueN(ref$group)
  p <- ggplot2::ggplot(ref, ggplot2::aes(UMAP1, UMAP2, color = group)) +
    ggplot2::geom_point(size = if (zoomed) 1.6 else 0.9, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal, name = if (color_by == "family") "Family" else "Class") +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = max(1, ceiling(n / 40)), override.aes = list(size = 2.5, alpha = 1)))
  if (nrow(fg) > 0) {
    fg[, lab := if (is.null(names)) Sample_ID_Full else unname(names[Sample_ID_Full])]
    p <- p + ggplot2::geom_point(data = fg, ggplot2::aes(UMAP1, UMAP2), inherit.aes = FALSE,
                                 shape = 8, size = 2.2, stroke = 0.9, color = "black")
    if (label_samples) p <- p + ggrepel::geom_text_repel(data = fg, ggplot2::aes(UMAP1, UMAP2, label = lab), inherit.aes = FALSE,
                                                         size = 3, fontface = "bold", color = "black", min.segment.length = 0,
                                                         box.padding = 0.4, max.overlaps = Inf, seed = 1)
  }
  sub <- paste0(nrow(ref), " reference samples, ", n, if (color_by == "family") " families" else " classes",
                if (nrow(fg) > 0) paste0("; asterisk = projected sample", if (nrow(fg) > 1) "s" else "") else "")
  p + ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = !zoomed) + ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(0.35, "cm"), legend.position = if (legend) "right" else "none") +
    ggplot2::labs(title = title, subtitle = sub, x = "UMAP 1", y = "UMAP 2")
}

# ------------------------------------------------------------------------------
# CLASS ZOOM
# Window = box containing ALL references of the chosen classes + the chosen samples.
# Class (text before the comma, e.g. "LGG") -> shape; subclass (full label) -> colour.
# ------------------------------------------------------------------------------
.mg_kelly <- c("#BE0032", "#0067A5", "#008856", "#F38400", "#875692", "#F3C300", "#E68FAC", "#A1CAF1", "#882D17",
               "#8DB600", "#604E97", "#F99379", "#2B3D26", "#DCD300", "#B3446C", "#654522", "#E25822", "#C2B280", "#F6A600")
.mg_shapes <- c(21, 22, 24, 23, 25)

magibu_zoom_plot <- function(meta, classes = character(0), samples = character(0), labels = NULL,
                             show_other = TRUE, label_classes = FALSE, label_samples = TRUE,
                             margin = 0.08, point_size = 3, title = NULL) {
  meta <- magibu_standard_meta(meta)
  ref <- meta[!magibu_is_new(meta)]
  fg  <- meta[magibu_is_new(meta) & Sample_ID_Full %in% samples]
  classes <- intersect(classes, unique(ref$sample.title))
  sel <- ref[sample.title %in% classes]
  box <- rbind(sel[, .(UMAP1, UMAP2)], fg[, .(UMAP1, UMAP2)])
  if (nrow(box) == 0) box <- ref[, .(UMAP1, UMAP2)]
  xr <- range(box$UMAP1); yr <- range(box$UMAP2)
  side <- max(diff(xr), diff(yr), 0.5)                      # never a degenerate window
  pad <- side * margin
  xlim <- c(xr[1] - pad, xr[2] + pad); ylim <- c(yr[1] - pad, yr[2] + pad)
  if (diff(xlim) < side) xlim <- mean(xr) + c(-1, 1) * (side / 2 + pad)
  if (diff(ylim) < side) ylim <- mean(yr) + c(-1, 1) * (side / 2 + pad)
  in_win <- function(d) d[UMAP1 >= xlim[1] & UMAP1 <= xlim[2] & UMAP2 >= ylim[1] & UMAP2 <= ylim[2]]
  # window wider than half the map: samples as small asterisks, no ring
  wide <- side > 0.5 * max(diff(range(ref$UMAP1)), diff(range(ref$UMAP2)))

  p <- ggplot2::ggplot()
  if (show_other) {
    oth <- in_win(ref[!sample.title %in% classes | is.na(sample.title)])
    p <- p + ggplot2::geom_point(data = oth, ggplot2::aes(UMAP1, UMAP2), color = "grey82", size = point_size * if (wide) 0.3 else 0.55)
  }
  if (nrow(sel) > 0) {
    lev <- data.table::data.table(cls = classes, fam = magibu_family(classes))
    lev <- lev[order(fam, cls)]
    fams <- unique(lev$fam)
    lev[, shape := .mg_shapes[(match(fam, fams) - 1) %% length(.mg_shapes) + 1]]
    lev[, color := .mg_kelly[(seq_len(.N) - 1) %% length(.mg_kelly) + 1]]
    lev[, n := vapply(cls, function(x) sum(sel$sample.title == x), integer(1))]
    lev[, lab := paste0(cls, " (", n, ")")]
    sel[, cls := factor(sample.title, levels = lev$cls)]
    pal <- stats::setNames(lev$color, lev$cls); shp <- stats::setNames(lev$shape, lev$cls); lb <- stats::setNames(lev$lab, lev$cls)
    p <- p +
      ggplot2::geom_point(data = sel, ggplot2::aes(UMAP1, UMAP2, fill = cls, shape = cls),
                          size = point_size * if (wide) 0.6 else 1, stroke = 0.35, color = "grey15", alpha = 0.9) +
      ggplot2::scale_fill_manual(values = pal, labels = lb, name = "Class (n)", drop = FALSE) +
      ggplot2::scale_shape_manual(values = shp, labels = lb, name = "Class (n)", drop = FALSE) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = if (nrow(lev) > 22) 2 else 1,
                                                   override.aes = list(size = 4, alpha = 1)))
    if (label_classes) {
      anchor <- function(x, y) { d <- as.matrix(stats::dist(cbind(x, y))); k <- max(1, ceiling(length(x) / 2))
        i <- which.min(apply(d, 1, function(r) sum(sort(r)[seq_len(k)]))); list(UMAP1 = x[i], UMAP2 = y[i]) }
      cen <- sel[, anchor(UMAP1, UMAP2), by = cls]
      p <- p + ggrepel::geom_label_repel(data = cen, ggplot2::aes(UMAP1, UMAP2, label = cls, color = cls),
                                         fill = ggplot2::alpha("white", 0.85), fontface = "bold", size = 3.4,
                                         label.size = 0.4, min.segment.length = 0, box.padding = 0.5,
                                         max.overlaps = Inf, seed = 1, show.legend = FALSE) +
        ggplot2::scale_color_manual(values = pal, guide = "none")
    }
  }
  if (nrow(fg) > 0) {
    fg[, lab := Sample_ID_Full]
    if (!is.null(labels)) fg[Sample_ID_Full %in% names(labels), lab := labels[Sample_ID_Full]]
    if (!wide) p <- p + ggplot2::geom_point(data = fg, ggplot2::aes(UMAP1, UMAP2), shape = 21, size = point_size * 2.6,
                                            stroke = 1.5, fill = NA, color = "black")
    p <- p + ggplot2::geom_point(data = fg, ggplot2::aes(UMAP1, UMAP2), shape = 8, size = point_size * if (wide) 0.75 else 1.2,
                                 stroke = if (wide) 0.9 else 1.3, color = "black")
    if (label_samples) p <- p + ggrepel::geom_label_repel(data = fg, ggplot2::aes(UMAP1, UMAP2, label = lab), fill = "black",
                                         color = "white", fontface = "bold", size = 3.6, point.padding = 1,
                                         box.padding = 0.8, min.segment.length = 0, segment.color = "black",
                                         max.overlaps = Inf, seed = 1)
  }
  sub <- paste0(length(classes), " classes highlighted: shape = class, colour = subclass",
                if (show_other) "; grey = other reference samples" else "",
                if (nrow(fg) > 0) paste0("\nAsterisk = projected sample", if (nrow(fg) > 1) "s" else "", " (", nrow(fg), ")") else "")
  p + ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.5, "cm"),
                   legend.text = ggplot2::element_text(size = 10)) +
    ggplot2::labs(title = title, subtitle = sub, x = "UMAP 1", y = "UMAP 2")
}

# Classes closest to the samples (to propose a first selection in the zoom)
magibu_nearest_classes <- function(dist, samples, n = 4) {
  d <- dist$Distances[Sample %in% samples & Rank <= n]
  unique(d[order(Rank)]$Class)
}
