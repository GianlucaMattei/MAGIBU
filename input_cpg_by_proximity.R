#' Impute Missing CpG Beta Values by Genomic Proximity
#'
#' For Nanopore samples with low read coverage, many CpG sites from the
#' reference (array-based training set) will not be directly observed. This
#' function recovers a subset of those missing sites by finding the nearest
#' observed CpG on the same chromosome within a user-defined distance
#' threshold. When found, the beta value of the nearest observed CpG is
#' assigned to the missing reference site.
#'
#' The rationale is that CpG methylation is spatially correlated along the
#' genome, especially within CpG islands and shores (< 1–2 kb). Using the
#' sample's own observed beta from a nearby CpG is typically more informative
#' than falling back to the population mean from the training set.
#'
#' @param sample_dt A \code{data.table} with the sample's observed CpG data
#'   after alignment. Required columns:
#'   \describe{
#'     \item{chr}{Character. Chromosome (e.g., \code{"chr1"}).}
#'     \item{pos}{Integer. Genomic position (already shift-corrected).}
#'     \item{Beta}{Numeric. Methylation beta value [0, 1].}
#'   }
#'   This is the sample's full set of observed CpGs (not just those matching
#'   the reference).
#' @param ref_positions A character vector of reference CpG position IDs in
#'   \code{"chr:pos"} format (e.g., \code{"chr1:10524"}). These are the
#'   model's \code{common_positions}.
#' @param max_distance Integer. Maximum genomic distance (in bp) to consider
#'   for proximity imputation. CpGs further than this threshold are left as
#'   \code{NA} (to be handled by downstream mean imputation). Default: 500.
#' @param verbose Logical. Print summary statistics? Default: \code{TRUE}.
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{ID}{Character. Position ID in \code{"chr:pos"} format.}
#'     \item{Beta}{Numeric. Imputed beta value from the nearest observed CpG.}
#'     \item{distance}{Integer. Genomic distance (bp) to the donor CpG.}
#'     \item{donor_pos}{Character. Position ID of the donor CpG.}
#'   }
#'   Only successfully imputed positions are returned (i.e., those within
#'   \code{max_distance}). Positions with no nearby CpG are excluded.
#'
#' @details
#' The algorithm works per-chromosome for efficiency:
#' \enumerate{
#'   \item Split reference positions and sample positions by chromosome.
#'   \item For each chromosome, identify reference positions NOT observed in
#'     the sample.
#'   \item For each missing reference position, perform a binary search on
#'     the sorted sample positions to find the nearest left and right
#'     observed CpGs.
#'   \item Select the closer of the two. If it is within \code{max_distance},
#'     assign its beta value.
#' }
#'
#' The function does NOT modify the original sample data. Its output is
#' designed to be \code{rbind}-ed with the direct-match results before
#' building the final aligned beta vector in \code{project_new_samples()}.
#'
#' @section Performance:
#' Uses \code{findInterval()} (binary search, O(log n) per query) rather
#' than brute-force distance calculation. Typical runtime is <1 second for
#' 50K reference sites even with 1M+ observed CpGs.
#'
#' @section Integration:
#' Call this function inside \code{read_and_align_single()} in
#' \code{project_new_samples.R}, after the direct position matching and
#' before building the final beta vector. Example:
#'
#' \preformatted{
#'   # After dt_filtered <- dt[ID %in% ref_pos_ids_vec, .(ID, Beta)]
#'   matched_ids   <- dt_filtered$ID
#'   missing_ids   <- setdiff(ref_pos_ids_vec, matched_ids)
#'
#'   if (length(missing_ids) > 0 && max_distance > 0) {
#'     imputed <- impute_cpg_by_proximity(
#'       sample_dt     = dt[, .(chr = chr_hg38, pos = pos_hg38, Beta)],
#'       ref_positions = missing_ids,
#'       max_distance  = max_distance
#'     )
#'     if (nrow(imputed) > 0) {
#'       dt_filtered <- rbind(dt_filtered, imputed[, .(ID, Beta)])
#'     }
#'   }
#' }
#'
#' @family projection
#' @export
impute_cpg_by_proximity <- function(sample_dt,
                                    ref_positions,
                                    max_distance = 500L,
                                    verbose = TRUE) {

  # =========================================================================
  # Input Validation
  # =========================================================================
  if (!data.table::is.data.table(sample_dt)) {
    if (is.data.frame(sample_dt)) {
      sample_dt <- data.table::as.data.table(sample_dt)
    } else {
      stop("`sample_dt` must be a data.table or data.frame, got: ",
           paste(class(sample_dt), collapse = ", "))
    }
  }

  required_cols <- c("chr", "pos", "Beta")
  missing_cols  <- setdiff(required_cols, names(sample_dt))
  if (length(missing_cols) > 0) {
    stop("`sample_dt` is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Expected: chr, pos, Beta.")
  }

  if (!is.character(ref_positions) || length(ref_positions) == 0) {
    stop("`ref_positions` must be a non-empty character vector of 'chr:pos' IDs.")
  }

  max_distance <- as.integer(max_distance)
  if (is.na(max_distance) || max_distance < 1) {
    stop("`max_distance` must be a positive integer (bp). Got: ", max_distance)
  }

  # =========================================================================
  # Preparation
  # =========================================================================

  # Clean sample data: remove NAs in position or Beta
  sample_clean <- sample_dt[!is.na(chr) & !is.na(pos) & !is.na(Beta)]

  if (nrow(sample_clean) == 0) {
    if (verbose) message("  [proximity] No valid sample CpGs. Nothing to impute.")
    return(data.table::data.table(
      ID        = character(0),
      Beta      = numeric(0),
      distance  = integer(0),
      donor_pos = character(0)
    ))
  }

  # Parse reference positions: "chr1:10524" -> chr + pos
  ref_dt <- data.table::data.table(ID = ref_positions)
  ref_dt[, c("chr", "pos") := data.table::tstrsplit(ID, ":", fixed = TRUE)]
  ref_dt[, pos := as.integer(pos)]

  # Remove any that failed to parse
  ref_dt <- ref_dt[!is.na(pos)]

  if (nrow(ref_dt) == 0) {
    if (verbose) message("  [proximity] No parseable reference positions.")
    return(data.table::data.table(
      ID        = character(0),
      Beta      = numeric(0),
      distance  = integer(0),
      donor_pos = character(0)
    ))
  }

  # =========================================================================
  # Per-Chromosome Proximity Search
  # =========================================================================

  # Split by chromosome
  sample_by_chr <- split(sample_clean, by = "chr", sorted = FALSE)
  ref_by_chr    <- split(ref_dt,       by = "chr", sorted = FALSE)

  # Only process chromosomes present in both sample and reference
  common_chrs <- intersect(names(ref_by_chr), names(sample_by_chr))

  if (length(common_chrs) == 0) {
    if (verbose) message("  [proximity] No shared chromosomes between sample and reference.")
    return(data.table::data.table(
      ID        = character(0),
      Beta      = numeric(0),
      distance  = integer(0),
      donor_pos = character(0)
    ))
  }

  result_list <- vector("list", length(common_chrs))

  for (i in seq_along(common_chrs)) {
    chr_name <- common_chrs[i]

    samp_chr <- sample_by_chr[[chr_name]]
    ref_chr  <- ref_by_chr[[chr_name]]

    # Sort sample positions for binary search
    data.table::setorder(samp_chr, pos)
    samp_pos  <- samp_chr$pos
    samp_beta <- samp_chr$Beta
    n_samp    <- length(samp_pos)

    # Reference missing positions on this chromosome
    ref_pos_vec <- ref_chr$pos
    ref_id_vec  <- ref_chr$ID
    n_ref       <- length(ref_pos_vec)

    # findInterval: for each ref position, find the index of the largest
    # sample position <= ref position. This gives us the left neighbor.
    left_idx <- findInterval(ref_pos_vec, samp_pos)
    # Right neighbor is left_idx + 1

    # Compute distances to left and right neighbors
    dist_left  <- rep(NA_integer_, n_ref)
    dist_right <- rep(NA_integer_, n_ref)
    beta_left  <- rep(NA_real_,    n_ref)
    beta_right <- rep(NA_real_,    n_ref)
    pos_left   <- rep(NA_integer_, n_ref)
    pos_right  <- rep(NA_integer_, n_ref)

    # Left neighbor (index > 0 means there is a CpG at or before this position)
    has_left <- left_idx > 0
    if (any(has_left)) {
      li <- left_idx[has_left]
      dist_left[has_left]  <- abs(ref_pos_vec[has_left] - samp_pos[li])
      beta_left[has_left]  <- samp_beta[li]
      pos_left[has_left]   <- samp_pos[li]
    }

    # Right neighbor (index < n_samp means there is a CpG after this position)
    right_idx <- left_idx + 1L
    has_right <- right_idx <= n_samp
    if (any(has_right)) {
      ri <- right_idx[has_right]
      dist_right[has_right] <- abs(samp_pos[ri] - ref_pos_vec[has_right])
      beta_right[has_right] <- samp_beta[ri]
      pos_right[has_right]  <- samp_pos[ri]
    }

    # Pick the closer neighbor
    # If tied, prefer left (arbitrary but deterministic)
    use_right <- !is.na(dist_right) &
      (is.na(dist_left) | dist_right < dist_left)

    best_dist <- ifelse(use_right, dist_right, dist_left)
    best_beta <- ifelse(use_right, beta_right, beta_left)
    best_pos  <- ifelse(use_right, pos_right,  pos_left)

    # Apply distance threshold
    keep <- !is.na(best_dist) & best_dist <= max_distance

    if (any(keep)) {
      result_list[[i]] <- data.table::data.table(
        ID        = ref_id_vec[keep],
        Beta      = best_beta[keep],
        distance  = as.integer(best_dist[keep]),
        donor_pos = paste0(chr_name, ":", best_pos[keep])
      )
    }
  }

  # =========================================================================
  # Combine and Report
  # =========================================================================
  imputed <- data.table::rbindlist(result_list, use.names = TRUE)

  if (verbose) {
    n_missing   <- nrow(ref_dt)
    n_imputed   <- nrow(imputed)
    n_remaining <- n_missing - n_imputed
    pct         <- round(n_imputed / n_missing * 100, 1)

    if (n_imputed > 0) {
      median_d <- median(imputed$distance)
      max_d    <- max(imputed$distance)
      message(sprintf(
        "  [proximity] Imputed %d / %d missing CpGs (%.1f%%) within %d bp | median dist: %d bp, max: %d bp",
        n_imputed, n_missing, pct, max_distance, median_d, max_d
      ))
    } else {
      message(sprintf(
        "  [proximity] 0 / %d missing CpGs had a neighbor within %d bp.",
        n_missing, max_distance
      ))
    }

    if (n_remaining > 0) {
      message(sprintf(
        "  [proximity] %d CpGs still missing -> will use training mean imputation.",
        n_remaining
      ))
    }
  }

  return(imputed)
}