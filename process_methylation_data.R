#' Process Methylation Data (Illumina IDAT or Nanopore)
#'
#' Processes raw methylation data from various platforms (EPICv2, EPICv1, 450K,
#' Nanopore) into standardised TSV files with a consistent schema:
#' \code{Probe_ID}, \code{Beta}, \code{chr_hg38}, \code{pos_hg38}.
#'
#' For Illumina arrays, IDAT files are read via \pkg{minfi}, beta values are
#' extracted, and genomic coordinates are added from annotation tables.
#' For Nanopore data, raw BED-like files are parsed using user-defined column
#' mappings.
#'
#' @param path.to.idat Character vector. One or more directory paths containing
#'   \code{.idat} files (for Illumina arrays), or a vector of specific file
#'   paths (for Nanopore). For Nanopore, can also be a single directory.
#' @param save.path Character. Directory where output TSV files will be saved.
#' @param array.type Character. One of \code{"EpicV2"}, \code{"Infinium"},
#'   \code{"EpicV1"}, or \code{"Nanopore"}.
#' @param annot.file.path Character or \code{NULL}. Path to an external
#'   annotation file (TSV/CSV with columns \code{Name}, \code{chr}, \code{pos}).
#'   Required for \code{array.type = "Infinium"} or \code{"EpicV1"}.
#' @param nanopore.settings Named character vector or \code{NULL}. Maps the
#'   column names \code{"chr"}, \code{"pos"}, \code{"beta"} to the actual column
#'   identifiers in the Nanopore files (e.g.
#'   \code{c(chr = "V1", pos = "V2", beta = "V5")}).
#'   Required for \code{array.type = "Nanopore"}.
#' @param file.pattern Character or \code{NULL}. A regex pattern to filter
#'   Nanopore files when \code{path.to.idat} is a directory.
#' @param n.cores Integer. Number of cores for parallel processing. Default is
#'   \code{min(100, parallel::detectCores() - 1)}.
#'
#' @return Called for its side effect of writing TSV files to \code{save.path}.
#'   Returns \code{NULL} invisibly.
#'
#' @details
#' \strong{EPICv2}: Uses \code{minfi::preprocessRaw} and the
#' \code{IlluminaHumanMethylationEPICv2anno.20a1.hg38} annotation package.
#' Processes samples in batches of 10.
#'
#' \strong{Infinium / EpicV1}: Uses \code{minfi::preprocessIllumina} with
#' background correction and control normalisation. Requires an external
#' annotation file via \code{annot.file.path}.
#'
#' \strong{Nanopore}: Reads raw BED-like files using \code{data.table::fread}
#' with column selection driven by \code{nanopore.settings}. Each file is saved
#' individually inside the parallel loop to minimise RAM usage.
#'
#' @importFrom data.table data.table fread fwrite setnames
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom tools file_path_sans_ext
#' @importFrom utils data
#' @export
process_methylation_data <- function(
    path.to.idat,
    save.path,
    array.type = "EpicV2",
    annot.file.path = NULL,
    nanopore.settings = c(chr = "V1", pos = "V2", beta = "V5"),
    file.pattern = NULL,
    n.cores = min(100, parallel::detectCores() - 1)
) {

  message("=== Starting Processing: ", array.type, " ===")

  if (n.cores < 1) n.cores <- 1
  doParallel::registerDoParallel(cores = n.cores)
  suppressPackageStartupMessages(require(minfi))
  message("Using ", n.cores, " cores.")

  if (!dir.exists(save.path)) dir.create(save.path, recursive = TRUE)

  final_results_list <- list()

  # ============================================================================
  # BLOCK: EPICv2
  # ============================================================================
  if (array.type == "EpicV2") {

    if (!requireNamespace("IlluminaHumanMethylationEPICv2anno.20a1.hg38", quietly = TRUE)) {
      stop("Package 'IlluminaHumanMethylationEPICv2anno.20a1.hg38' is required but not installed.")
    }

    message("Loading EPICv2 annotation...")

    # NOTE: minfi::getAnnotation() requires the annotation package to be
    # attached to the search path. loadNamespace() and utils::data(envir=)
    # are NOT sufficient — minfi internally dispatches via the search list.
    # library() is the only reliable approach here.
    suppressPackageStartupMessages(
      library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    )
    data("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
    anno_data <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

    annot_df <- data.table::data.table(
      Name = rownames(anno_data),
      chr  = anno_data$chr,
      pos  = anno_data$pos
    )

    targets <- .read_targets_safe(path.to.idat)
    batch.size <- 10
    target_batches <- split(targets, ceiling(seq_len(nrow(targets)) / batch.size))

    `%dopar%` <- foreach::`%dopar%`
    final_results_list <- foreach::foreach(
      batch = target_batches,
      .packages = c("minfi", "data.table")
    ) %dopar% {
      .process_batch_epicv2(batch, annot_df)
    }
    final_results_list <- unlist(final_results_list, recursive = FALSE)

  # ============================================================================
  # BLOCK: INFINIUM / EPICv1
  # ============================================================================
  } else if (array.type == "Infinium" || array.type == "EpicV1") {

    if (is.null(annot.file.path)) {
      stop("ERROR: 'annot.file.path' is required for array.type 'Infinium' or 'EpicV1'.")
    }

    message("Loading external annotation...")
    annot_ext <- data.table::fread(annot.file.path)

    targets <- .read_targets_safe(path.to.idat)
    batch.size <- 10
    target_batches <- split(targets, ceiling(seq_len(nrow(targets)) / batch.size))

    final_results_list <- list()
    for (i in seq_along(target_batches)) {
      message("Batch ", i, "/", length(target_batches))
      batch_res <- .process_batch_infinium_legacy(target_batches[[i]], annot_ext)
      final_results_list <- c(final_results_list, batch_res)
    }
  # ============================================================================
  # BLOCK: NANOPORE
  # ============================================================================
  } else if (array.type == "Nanopore") {

    if (is.null(nanopore.settings)) stop("'nanopore.settings' must be specified!")

    if (dir.exists(path.to.idat[1])) {
      # Case A: directory — apply file.pattern if provided
      files <- list.files(path.to.idat, pattern = file.pattern,
                          full.names = TRUE, recursive = TRUE)
    } else {
      # Case B: vector of specific file paths
      files <- path.to.idat
      files <- files[file.exists(files)]
    }

    if (length(files) == 0) stop("No files found matching the specified criteria!")
    message("Found ", length(files), " Nanopore files. Starting stream processing...")

    # Save directly inside the parallel loop to avoid RAM saturation
    `%dopar%` <- foreach::`%dopar%`
    foreach::foreach(fp = files, .packages = c("data.table"), .export = ".process_nanopore_files") %dopar% {      
      res_list <- .process_nanopore_files(c(fp), nanopore.settings)

      sample_name <- names(res_list)[1]
      dt <- res_list[[1]]

      if (!is.null(dt) && nrow(dt) > 0) {
        clean_name <- tools::file_path_sans_ext(sample_name)
        out_file <- file.path(save.path, paste0(clean_name, ".tsv"))
        data.table::fwrite(dt, file = out_file, sep = "\t")
      }
      return(NULL)
    }

    message("Nanopore processing completed. Files saved to: ", save.path)
    final_results_list <- list()

  } else {
    stop("Unrecognised array type: '", array.type, "'")
  }

  doParallel::stopImplicitCluster()

  # ============================================================================
  # FINAL SAVE (Only for Epic/Infinium — Nanopore saves inside the loop)
  # ============================================================================
  if (length(final_results_list) > 0) {
    message("Saving TSV files to: ", save.path)
    count <- 0
    for (sample_name in names(final_results_list)) {
      dt_res <- final_results_list[[sample_name]]
      if (!is.null(dt_res) && nrow(dt_res) > 0) {
        clean_name <- gsub("_Grn$|_Red$", "", sample_name)
        out_file <- file.path(save.path, paste0(clean_name, ".tsv"))
        data.table::fwrite(dt_res, file = out_file, sep = "\t")
        count <- count + 1
      }
    }
    message("Done! Generated ", count, " TSV files.")
  }

  invisible(NULL)
}


# ==============================================================================
# Internal Helper Functions (not exported)
# ==============================================================================

#' Read IDAT Targets from Directories
#'
#' Scans one or more directories for \code{_Grn.idat} files and builds a
#' targets data.frame suitable for \code{minfi::read.metharray.exp}.
#'
#' @param base_dirs Character vector of directory paths to search.
#' @return A data.frame with columns \code{Basename} and \code{Sample_Name}.
#' @keywords internal
.read_targets_safe <- function(base_dirs) {
  valid_paths <- base_dirs[dir.exists(base_dirs)]

  if (length(valid_paths) == 0) {
    stop("None of the specified directories exist: ",
         paste(base_dirs, collapse = ", "))
  }

  grn_files <- list.files(valid_paths, pattern = "_Grn.idat$",
                          recursive = TRUE, full.names = TRUE)

  if (length(grn_files) == 0) {
    stop("No IDAT files found in: ", paste(valid_paths, collapse = ", "))
  }

  basenames <- sub("_Grn.idat$", "", grn_files)

  targets <- data.frame(
    Basename    = basenames,
    Sample_Name = basename(basenames),
    stringsAsFactors = FALSE
  )
  return(targets)
}


#' Process a Batch of EPICv2 Samples
#'
#' Reads a batch of IDAT files, extracts beta values via
#' \code{minfi::preprocessRaw}, and merges with EPICv2 annotation.
#'
#' @param batch.targets A data.frame of targets (subset of full targets).
#' @param annot.epicv2 A \code{data.table} with columns \code{Name},
#'   \code{chr}, \code{pos}.
#' @return A named list of \code{data.table} objects, one per sample.
#' @keywords internal
.process_batch_epicv2 <- function(batch.targets, annot.epicv2) {
  rgSet <- minfi::read.metharray.exp(targets = batch.targets,
                                     force = TRUE)

  if (minfi::annotation(rgSet)["array"] == "Unknown") {
    minfi::annotation(rgSet) <- c(array = "IlluminaHumanMethylationEPICv2",
                                  annotation = "20a1.hg38")
  }

  mSet <- minfi::preprocessRaw(rgSet)
  beta_vals <- minfi::getBeta(mSet)

  dt_list <- list()
  for (j in seq_len(ncol(beta_vals))) {
    sample_name <- colnames(beta_vals)[j]
    sample_beta <- beta_vals[, j, drop = FALSE]

    dt <- data.table::data.table(
      Probe_ID = rownames(sample_beta),
      Beta     = as.numeric(sample_beta)
    )
    dt <- merge(dt, annot.epicv2, by.x = "Probe_ID", by.y = "Name", all.x = TRUE)

    data.table::setnames(dt, old = c("chr", "pos"),
                         new = c("chr_hg38", "pos_hg38"),
                         skip_absent = TRUE)

    cols_needed <- c("Probe_ID", "Beta", "chr_hg38", "pos_hg38")
    dt <- dt[, ..cols_needed]

    dt_list[[sample_name]] <- dt
  }
  return(dt_list)
}


#' Process a Batch of Infinium/EPICv1 Samples
#'
#' Reads a batch of IDAT files, applies Illumina preprocessing (background
#' correction + control normalisation), and merges with an external annotation.
#'
#' @param batch.targets A data.frame of targets (subset of full targets).
#' @param annot.external A \code{data.table} with columns \code{Name},
#'   \code{chr}, \code{pos}.
#' @return A named list of \code{data.table} objects, one per sample.
#' @keywords internal
.process_batch_infinium_legacy <- function(batch.targets, annot.external) {
  rgSet <- minfi::read.metharray.exp(targets = batch.targets,
                                     force = TRUE)
  mSet <- minfi::preprocessIllumina(rgSet, bg.correct = TRUE,
                                    normalize = "controls")
  beta_values <- minfi::getBeta(mSet)

  dt_list <- list()
  for (j in seq_len(ncol(beta_values))) {
    sample_name <- colnames(beta_values)[j]
    sample_beta <- beta_values[, j, drop = FALSE]

    dt <- data.table::data.table(
      Probe_ID = rownames(sample_beta),
      Beta     = as.numeric(sample_beta)
    )
    dt <- merge(dt, annot.external, by.x = "Probe_ID", by.y = "Name", all.x = TRUE)

    if ("chr" %in% names(dt)) data.table::setnames(dt, "chr", "chr_hg38")
    if ("pos" %in% names(dt)) data.table::setnames(dt, "pos", "pos_hg38")

    cols_needed <- c("Probe_ID", "Beta", "chr_hg38", "pos_hg38")
    missing_cols <- setdiff(cols_needed, names(dt))
    if (length(missing_cols) > 0) dt[, (missing_cols) := NA]

    dt <- dt[, ..cols_needed]
    dt_list[[sample_name]] <- dt
  }
  return(dt_list)
}


#' Process Nanopore Methylation Files
#'
#' Reads raw Nanopore BED-like files using user-defined column mappings,
#' rounds beta values, and constructs a standardised output table.
#'
#' @param file_paths Character vector of file paths.
#' @param col_settings Named character vector mapping \code{"chr"}, \code{"pos"},
#'   \code{"beta"} to column identifiers in the files (e.g. \code{"V1"},
#'   \code{"V2"}, \code{"V5"}).
#' @return A named list of \code{data.table} objects, one per file.
#' @keywords internal
.process_nanopore_files <- function(file_paths, col_settings) {
  dt_list <- list()

  cols_to_keep <- as.character(col_settings)

  for (fp in file_paths) {
    sample_name <- tools::file_path_sans_ext(basename(fp))

    # Optimised read: select only the needed columns.
    # Note: tryCatch returns the value so raw_dt is properly assigned.
    raw_dt <- tryCatch({
      data.table::fread(fp, select = cols_to_keep, header = FALSE)
    }, error = function(e) {
      warning("Optimised read failed, falling back to full read: ", fp)
      data.table::fread(fp, header = FALSE)
    })

    available_cols <- names(raw_dt)
    if (!all(cols_to_keep %in% available_cols)) {
      # Retry with header = TRUE
      raw_dt <- data.table::fread(fp, select = cols_to_keep, header = TRUE)
      if (!all(cols_to_keep %in% names(raw_dt))) {
        warning(paste("File skipped (missing columns):", fp))
        next
      }
    }

    dt <- data.table::data.table(
      chr_hg38 = raw_dt[[col_settings["chr"]]],
      pos_hg38 = raw_dt[[col_settings["pos"]]],
      Beta     = raw_dt[[col_settings["beta"]]]
    )
    
    dt[, Beta := suppressWarnings(as.numeric(Beta))]
    dt <- dt[!is.na(Beta)]
    
    # Round beta values (nanopore has excessive decimal places)
    dt[, Beta := round(Beta, 4)]

    # Create Probe_ID
    dt[, Probe_ID := paste0(chr_hg38, ":", pos_hg38)]

    # Final column selection
    cols_needed <- c("Probe_ID", "Beta", "chr_hg38", "pos_hg38")
    dt <- dt[, ..cols_needed]

    dt_list[[sample_name]] <- dt

    rm(raw_dt)
  }
  return(dt_list)
}