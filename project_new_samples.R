#' Project New Samples onto Reference UMAP
#'
#' Projects new methylation samples (array or Nanopore) onto an existing PCA + UMAP
#' space. Auto-detects whether new samples use probe names (\code{cg*}) or genomic
#' positions (\code{chr:pos}), and always matches via genomic coordinates.
#'
#' For low-coverage Nanopore samples, missing CpGs can be recovered by proximity
#' imputation: if a reference CpG is not directly observed, the beta value of the
#' nearest observed CpG on the same chromosome (within \code{max_distance} bp) is used.
#'
#' @param training.metadata Path to training metadata TSV, or a data.frame/data.table.
#' @param model.components Path to \code{.rds} model file, or the list object.
#' @param uwot.model Path to \code{.uwot} model file, or the uwot object.
#' @param new.samples File paths, a single directory, or a single TSV path.
#' @param save.path Directory for output files.
#' @param name Prefix for output filenames. Default: \code{"projected"}.
#' @param path.new.tsv (Deprecated) Use \code{new.samples} instead.
#' @param annotation_file Optional path to clinical annotation TSV for metadata merge.
#' @param max_distance Integer. Maximum distance (bp) for proximity-based CpG
#'   imputation. Set to 0 to disable. Default: 500.
#'
#' @return A list with \code{metadata} (combined training + projected) and
#'   \code{new_umap_coords} (matrix of UMAP coordinates for new samples).
#'
#' @export
project_new_samples2 <- function(
    training.metadata,
    model.components,
    uwot.model,
    new.samples = NULL,
    save.path,
    name = "projected",
    path.new.tsv = NULL,
    annotation_file = NULL,
    max_distance = 25L
) {
 message("=== STARTING PROJECTION OF NEW SAMPLES ===")

  # --- 0. Backwards Compatibility ---
  if (is.null(new.samples)) {
    if (!is.null(path.new.tsv)) {
      new.samples <- path.new.tsv
    } else {
      stop("You must specify the file or folder of new samples in 'new.samples'.")
    }
  }

  if (!dir.exists(save.path)) dir.create(save.path, recursive = TRUE)

  max_distance <- as.integer(max_distance)

  # =========================================================================
  # 1. Load Model
  # =========================================================================
  message("Loading model components...")
  if (is.character(model.components)) {
    if (!file.exists(model.components)) stop("model.components file not found.")
    model_obj <- readRDS(model.components)
  } else {
    model_obj <- model.components
  }

  message("Loading UMAP model (uwot)...")
  if (is.character(uwot.model)) {
    if (!file.exists(uwot.model)) stop("uwot model file not found.")
    umap_model_loaded <- uwot::load_uwot(uwot.model)
  } else {
    umap_model_loaded <- uwot.model
  }

  # --- Validate model contents ---
  if (is.null(model_obj$pca_model)) {
    # Fallback for old format (rotation/center/scale separate)
    if (!is.null(model_obj$rotation) && !is.null(model_obj$center)) {
      message("Detected OLD model format (rotation + center). Wrapping into prcomp-like object.")
      model_obj$pca_model <- list(
        rotation = model_obj$rotation,
        center   = model_obj$center,
        scale    = model_obj$scale
      )
      class(model_obj$pca_model) <- "prcomp"
    } else {
      stop("Model components missing 'pca_model'. Regenerate model with updated generate_umap_training_model().")
    }
  }

  # =========================================================================
  # 2. Input File Handling
  # =========================================================================
  files <- NULL
  if (length(new.samples) > 1) {
    files <- new.samples
  } else {
    if (grepl("\\.tsv$", new.samples, ignore.case = TRUE)) {
      files <- new.samples
    } else {
      if (!dir.exists(new.samples)) stop("Directory not found.")
      files <- list.files(new.samples, pattern = "\\.tsv$", full.names = TRUE)
    }
  }

  if (length(files) == 0) stop("No files found.")
  message("Found ", length(files), " files to project.")

  # =========================================================================
  # 3. Prepare Reference IDs (always chr:pos)
  # =========================================================================
  selected_cpgs    <- model_obj$selected_cpgs
  common_positions <- model_obj$common_positions
  probe_to_pos     <- model_obj$probe_to_position  # data.table or NULL

  # Build the reference ID vector
  if (!is.null(common_positions) && sum(!is.na(common_positions)) > 0) {
    ref_pos_ids <- common_positions[!is.na(common_positions)]
    pos_to_cpg_idx <- stats::setNames(
      match(names(ref_pos_ids), selected_cpgs),
      as.character(ref_pos_ids)
    )
    n_total_features <- length(selected_cpgs)
    message(paste("Reference: ", length(ref_pos_ids), "/", n_total_features,
                  " probes have genomic positions for matching."))
  } else if (!is.null(probe_to_pos)) {
    pos_map <- stats::setNames(probe_to_pos$pos_id, probe_to_pos$Probe_ID)
    ref_pos_ids <- pos_map[selected_cpgs]
    ref_pos_ids <- ref_pos_ids[!is.na(ref_pos_ids)]
    pos_to_cpg_idx <- stats::setNames(
      match(names(ref_pos_ids), selected_cpgs),
      as.character(ref_pos_ids)
    )
    n_total_features <- length(selected_cpgs)
    message(paste("Reference (from lookup):", length(ref_pos_ids), "/",
                  n_total_features, "with positions."))
  } else {
    stop("Model has no position information (common_positions or probe_to_position). ",
         "Regenerate the model with the updated generate_umap_training_model().")
  }

  ref_pos_ids_vec <- as.character(ref_pos_ids)

  message("\n--- ID DIAGNOSTICS ---")
  message("Example model position ID: ", ref_pos_ids_vec[1])
  message("Total reference positions: ", length(ref_pos_ids_vec))
  if (max_distance > 0) {
    message("Proximity imputation enabled: max_distance = ", max_distance, " bp")
  } else {
    message("Proximity imputation DISABLED (max_distance = 0)")
  }

  # =========================================================================
  # 4. Sample Processing
  # =========================================================================

  read_and_align_single <- function(fp) {
    dt <- data.table::fread(fp)

    # --- AUTO-DETECT: Array (cg*) vs Nanopore (chr:pos) ---
    has_probe_id <- "Probe_ID" %in% names(dt)
    has_positions <- all(c("chr_hg38", "pos_hg38") %in% names(dt))

    is_array <- FALSE
    if (has_probe_id) {
      sample_ids <- head(dt$Probe_ID[!is.na(dt$Probe_ID)], 20)
      if (any(grepl("^cg", sample_ids))) {
        is_array <- TRUE
      }
    }

    # Will hold the sample's full observed CpGs (for proximity imputation)
    sample_all_positions <- NULL

    if (is_array && !has_positions && !is.null(probe_to_pos)) {
      # =============================================================
      # PATH A: ARRAY without positions — convert cg* -> chr:pos
      # =============================================================

      dt_merged <- merge(dt, probe_to_pos[, .(Probe_ID, pos_id, chr_hg38, pos_hg38)],
                         by = "Probe_ID", all.x = FALSE)

      if (nrow(dt_merged) == 0) {
        warning(paste("No probe names matched for", basename(fp), "- check Probe_ID format."))
        return(NULL)
      }

      dt_merged[, ID := pos_id]

      beta_col <- intersect(c("Beta", "beta_value"), names(dt_merged))
      if (length(beta_col) == 0) {
        beta_candidates <- grep("beta", names(dt_merged), ignore.case = TRUE, value = TRUE)
        if (length(beta_candidates) > 0) beta_col <- beta_candidates[1]
        else { warning(paste("No Beta column in", basename(fp))); return(NULL) }
      } else {
        beta_col <- beta_col[1]
      }

      dt_filtered <- dt_merged[ID %in% ref_pos_ids_vec, .(ID, Beta = get(beta_col))]

      # For proximity: all observed positions from this sample
      sample_all_positions <- dt_merged[!is.na(chr_hg38) & !is.na(pos_hg38),
                                        .(chr = chr_hg38, pos = pos_hg38, Beta = get(beta_col))]

      n_matched <- nrow(dt_filtered)
      pct_cov <- (n_matched / length(ref_pos_ids_vec)) * 100
      message(paste0(basename(fp), " [ARRAY] -> Matched: ", n_matched,
                     " probes (", round(pct_cov, 1), "%)"))

    } else if (has_positions) {
      # =============================================================
      # PATH B: NANOPORE / POSITION-BASED — auto-align (+1/-1)
      # =============================================================

      if (!"Beta" %in% names(dt)) {
        beta_candidates <- grep("beta", names(dt), ignore.case = TRUE, value = TRUE)
        if (length(beta_candidates) > 0) {
          data.table::setnames(dt, beta_candidates[1], "Beta")
        } else {
          warning(paste("No Beta column in", basename(fp)))
          return(NULL)
        }
      }

      # Try 3 alignment strategies
      dt[, ID_raw    := paste0(chr_hg38, ":", pos_hg38)]
      dt[, ID_plus1  := paste0(chr_hg38, ":", pos_hg38 + 1)]
      dt[, ID_minus1 := paste0(chr_hg38, ":", pos_hg38 - 1)]

      matches_raw    <- sum(dt$ID_raw    %in% ref_pos_ids_vec)
      matches_plus1  <- sum(dt$ID_plus1  %in% ref_pos_ids_vec)
      matches_minus1 <- sum(dt$ID_minus1 %in% ref_pos_ids_vec)

      scores <- c("Raw" = matches_raw, "Shift+1" = matches_plus1, "Shift-1" = matches_minus1)
      best_strat <- names(which.max(scores))
      best_score <- max(scores)
      pct_cov <- (best_score / length(ref_pos_ids_vec)) * 100

      shift <- 0L
      if (best_strat == "Shift+1") {
        dt[, ID := ID_plus1]
        shift <- 1L
        msg_suffix <- " (SHIFT +1)"
      } else if (best_strat == "Shift-1") {
        dt[, ID := ID_minus1]
        shift <- -1L
        msg_suffix <- " (SHIFT -1)"
      } else {
        dt[, ID := ID_raw]
        msg_suffix <- ""
      }

      data_type <- if (is_array) "ARRAY+POS" else "NANOPORE"
      message(paste0(basename(fp), " [", data_type, "] -> Overlap: ", best_score,
                     " (", round(pct_cov, 1), "%)", msg_suffix))

      if (pct_cov < 1) {
        warning(paste("CRITICAL LOW COVERAGE for", basename(fp),
                      "-> CHECK GENOME REFERENCE (hg19 vs hg38)!"))
      }

      dt_filtered <- dt[ID %in% ref_pos_ids_vec, .(ID, Beta)]

      # For proximity: all observed positions (with shift applied)
      sample_all_positions <- dt[!is.na(chr_hg38) & !is.na(pos_hg38),
                                 .(chr = chr_hg38, pos = pos_hg38 + shift, Beta)]

    } else {
      warning(paste("File", basename(fp),
                    "has neither cg* Probe_IDs (without position cols) nor chr_hg38/pos_hg38. Skipping."))
      return(NULL)
    }

    # =================================================================
    # PROXIMITY IMPUTATION: recover missing CpGs from nearby observed
    # =================================================================
    if (max_distance > 0 && !is.null(sample_all_positions) && nrow(sample_all_positions) > 0) {
      matched_ids <- dt_filtered$ID
      missing_ids <- setdiff(ref_pos_ids_vec, matched_ids)

      if (length(missing_ids) > 0) {
        imputed <- impute_cpg_by_proximity(
          sample_dt     = sample_all_positions,
          ref_positions = missing_ids,
          max_distance  = max_distance,
          verbose       = TRUE
        )
        if (nrow(imputed) > 0) {
          dt_filtered <- data.table::rbindlist(
            list(dt_filtered, imputed[, .(ID, Beta)]),
            use.names = TRUE
          )
        }
      }
    }

    # --- Build ordered beta vector (aligned to model features) ---
    result_vec <- rep(NA_real_, n_total_features)
    names(result_vec) <- selected_cpgs

    # Handle duplicates (multiple reads per site -> take mean)
    if (anyDuplicated(dt_filtered$ID)) {
      dt_filtered <- dt_filtered[, .(Beta = mean(Beta, na.rm = TRUE)), by = ID]
    }

    # Map position IDs back to feature indices
    matched_idx <- pos_to_cpg_idx[dt_filtered$ID]
    valid <- !is.na(matched_idx)
    result_vec[matched_idx[valid]] <- dt_filtered$Beta[valid]

    return(result_vec)
  }

  message("\nProcessing and Aligning...")
  if (.Platform$OS.type == "unix") {
    n_cores <- min(10, parallel::detectCores() - 1)
    res_list <- parallel::mclapply(files, read_and_align_single, mc.cores = n_cores)
  } else {
    res_list <- lapply(files, read_and_align_single)
  }

  names(res_list) <- tools::file_path_sans_ext(basename(files))
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]

  if (length(res_list) == 0) stop("No samples were successfully processed.")

  new_matrix <- do.call(rbind, res_list)

  # =========================================================================
  # 5. Mean Imputation (remaining NAs after proximity)
  # =========================================================================
  na_indices <- which(is.na(new_matrix), arr.ind = TRUE)
  if (nrow(na_indices) > 0) {
    pct_na <- (nrow(na_indices) / prod(dim(new_matrix))) * 100
    message(paste0("\nMEAN IMPUTATION: ", round(pct_na, 1), "% still missing after proximity."))

    if ("col_means_training" %in% names(model_obj)) {
      message("Using training column means for imputation.")
      col_means <- model_obj$col_means_training
      col_indices <- na_indices[, 2]
      new_matrix[na_indices] <- col_means[col_indices]
    } else {
      message("Using PCA center for imputation (missing CpGs contribute zero to projection).")
      pca_center <- model_obj$pca_model$center
      col_indices <- na_indices[, 2]
      new_matrix[na_indices] <- pca_center[col_indices]
    }
  }

  # =========================================================================
  # 6. Projection
  # =========================================================================
  message("Projecting onto PCA space...")
  pca_pred <- stats::predict(model_obj$pca_model, new_matrix)

  # Subset to the number of PCs used during UMAP training
  n_pcs <- if (!is.null(model_obj$pca_components)) {
    model_obj$pca_components
  } else {
    50L  # default fallback
  }
  if (ncol(pca_pred) > n_pcs) {
    pca_pred <- pca_pred[, 1:n_pcs, drop = FALSE]
  }

  message("Projecting onto UMAP space...")
  umap_pred <- uwot::umap_transform(pca_pred, umap_model_loaded)

  # =========================================================================
  # 7. Output
  # =========================================================================
  message("Preparing output...")

  if (is.character(training.metadata)) {
    if (!file.exists(training.metadata)) stop("Training metadata file not found.")
    meta.training <- data.table::fread(training.metadata)
  } else {
    meta.training <- training.metadata
  }

  meta_new <- data.frame(
    Sample = rownames(new_matrix),
    Sample_ID = rownames(new_matrix),
    Sample_ID_Full = rownames(new_matrix),
    UMAP1 = umap_pred[, 1],
    UMAP2 = umap_pred[, 2],
    Batch = "New_Projected",
    stringsAsFactors = FALSE
  )

  # Carry over class columns from training
  candidate_cols <- c("sample.title", "general_class", "subclass",
                       "diagnosis", "methylation_class")
  found_class_cols <- intersect(candidate_cols, names(meta.training))
  for (col in found_class_cols) {
    if (!col %in% names(meta_new)) meta_new[[col]] <- NA_character_
  }

  meta.final <- dplyr::bind_rows(meta.training, meta_new)

  # --- Optional Clinical Annotation Merge ---
  if (!is.null(annotation_file)) {
    message("Merging clinical annotations from: ", annotation_file)
    if (!file.exists(annotation_file)) stop("Annotation file not found: ", annotation_file)
    ann <- data.table::fread(annotation_file)
    meta.final <- data.table::as.data.table(meta.final)
    meta.final[, join_key := sub("_.*", "", Sample_ID_Full)]
    if ("sample.title" %in% names(meta.final)) meta.final[, sample.title := NULL]
    meta.final <- merge(meta.final, ann, by.x = "join_key", by.y = names(ann)[1], all.x = TRUE)
    meta.final[, join_key := NULL]
    message("Annotation merge complete. Added ", ncol(ann) - 1, " columns.")
  }

  file_meta <- file.path(save.path, paste0(name, "_metadata_combined.tsv"))
  data.table::fwrite(meta.final, file_meta, sep = "\t")

  message("\n=== PROJECTION COMPLETED ===")
  message("Results saved to: ", save.path)

  return(list(metadata = meta.final, new_umap_coords = umap_pred))
}
