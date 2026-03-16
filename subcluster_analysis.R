# ==============================================================================
# SUBCLUSTER-BASED DISTANCE ANALYSIS FOR METHYLATION UMAP
# ==============================================================================
#
# This module implements a subcluster-aware distance calculation strategy.
# Instead of computing KNN distances across an entire class (which can be
# misleading when a class is split into spatially distinct groups on the UMAP),
# it first identifies spatial subclusters within each reference class using
# HDBSCAN, then computes mean distances from a target sample to each subcluster.
#
# Dependencies: dbscan, data.table, dplyr, ggplot2, patchwork, ggrepel, scales
# ==============================================================================

# --- Helper: Identify Class Column ---
# Shared logic used across all functions in this module.
# Returns the name of the first matching class column found in the data.
.find_class_column <- function(dt) {
    candidates <- c("methylation_class", "diagnosis", "sample.title",
                     "general_class", "subclass")
    found <- intersect(candidates, names(dt))
    if (length(found) == 0) return(NULL)
    return(found[1])
}

# --- Helper: Split Reference vs Target ---
# Separates reference samples from projected (new) samples.
.split_ref_target <- function(meta, sample_ids = NULL) {
    if ("Batch" %in% names(meta)) {
        ref <- meta[is.na(Batch) | Batch != "New_Projected"]
        new <- meta[Batch == "New_Projected"]
     } else if (!is.null(sample_ids) && length(sample_ids) > 0) {
        ref <- meta[!Sample_ID_Full %in% sample_ids]
        new <- meta[Sample_ID_Full %in% sample_ids]
    } else {
        # No Batch, no sample_ids -> all samples are reference
        ref <- meta
        new <- meta[0]
    }
    list(ref = ref, new = new)
}


#' Subcluster Reference Classes via HDBSCAN
#'
#' For each reference class in the metadata, identifies spatial subclusters on the
#' UMAP embedding using HDBSCAN. Classes with few samples are treated as a single
#' cluster. Noise points (HDBSCAN cluster 0) are reassigned to the nearest
#' subcluster centroid.
#'
#' @param metadata_input Character path or data.frame/data.table. The combined
#'   metadata containing UMAP1, UMAP2, Sample_ID_Full, and a class column.
#' @param sample_ids Character vector. Optional IDs of target (new) samples to
#'   exclude from the reference set. If NULL, uses the Batch column.
#' @param min_cluster_size Integer. Minimum number of points to form a subcluster
#'   (passed to \code{dbscan::hdbscan} as \code{minPts}). Default is 3. Must be >= 2.
#' @param min_samples_for_split Integer. Classes with fewer samples than this will
#'   not be split (treated as a single cluster). Default is 4. Must be > min_cluster_size.
#'
#' @return A \code{data.table} with one row per reference sample, containing:
#'   \describe{
#'     \item{Sample_ID_Full}{Character. Sample identifier.}
#'     \item{UMAP1, UMAP2}{Numeric. UMAP coordinates.}
#'     \item{Class_Label}{Character. Original class label.}
#'     \item{Subcluster_ID}{Integer. Subcluster index within the class (1, 2, ...).}
#'     \item{Subcluster_Label}{Character. Human-readable label, e.g. "Ependymoma [A]"
#'       or "Ependymoma" if only one subcluster.}
#'     \item{Subcluster_N}{Integer. Number of samples in this subcluster.}
#'     \item{Centroid_UMAP1, Centroid_UMAP2}{Numeric. Centroid of the subcluster.}
#'   }
#'
#' @details
#' The algorithm proceeds as follows for each class:
#' \enumerate{
#'   \item If the class has fewer than \code{min_samples_for_split} samples,
#'     assign all to subcluster 1 (no splitting).
#'   \item Otherwise, run \code{dbscan::hdbscan()} on UMAP coordinates.
#'   \item If HDBSCAN finds only noise (all cluster = 0) or a single cluster,
#'     assign all to subcluster 1.
#'   \item Reassign noise points (cluster = 0) to the subcluster whose centroid
#'     is closest in Euclidean distance.
#'   \item Label subclusters alphabetically by size (largest = A).
#' }
#'
#' @export
#' @family subcluster
subcluster_reference_classes <- function(metadata_input,
                                          sample_ids = NULL,
                                          min_cluster_size = 3,
                                          min_samples_for_split = 4) {

    # --- Input Validation ---
    if (!requireNamespace("dbscan", quietly = TRUE)) {
        stop("Package 'dbscan' is required for HDBSCAN subclustering. ",
             "Install it with: install.packages('dbscan')")
    }

    if (min_cluster_size < 2) {
        stop("min_cluster_size must be >= 2. HDBSCAN requires at least 2 points per cluster.")
    }
    if (min_samples_for_split <= min_cluster_size) {
        stop("min_samples_for_split must be > min_cluster_size. ",
             "Otherwise classes cannot be meaningfully split.")
    }

    # Load data
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found: ", metadata_input)
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }

    req_cols <- c("UMAP1", "UMAP2", "Sample_ID_Full")
    missing <- setdiff(req_cols, names(meta))
    if (length(missing) > 0) {
        stop("Missing required columns in metadata: ", paste(missing, collapse = ", "))
    }

    # Separate reference from target
    split_data <- .split_ref_target(meta, sample_ids)
    ref <- split_data$ref

    if (nrow(ref) == 0) stop("Reference set is empty after splitting.")

    # Identify class column
    class_col <- .find_class_column(ref)
    if (is.null(class_col)) {
        stop("No class column found in metadata. Expected one of: ",
             "methylation_class, diagnosis, sample.title, general_class, subclass")
    }
    ref$Class_Label <- as.character(ref[[class_col]])

    # Remove samples with NA class
    n_before <- nrow(ref)
    ref <- ref[!is.na(Class_Label) & Class_Label != ""]
    n_removed <- n_before - nrow(ref)
    if (n_removed > 0) {
        message("Removed ", n_removed, " reference samples with missing class labels.")
    }

    # --- Subclustering Loop ---
    classes <- unique(ref$Class_Label)
    message("Subclustering ", length(classes), " reference classes (minPts=", min_cluster_size, ")...")

    results_list <- vector("list", length(classes))

    for (i in seq_along(classes)) {
        cls <- classes[i]
        idx <- which(ref$Class_Label == cls)
        n_cls <- length(idx)
        coords <- as.matrix(ref[idx, .(UMAP1, UMAP2)])

        if (n_cls < min_samples_for_split) {
            # Too few samples to split — single cluster
            cluster_ids <- rep(1L, n_cls)
            n_subclusters <- 1L
        } else {
            # Run HDBSCAN
            hdb_res <- dbscan::hdbscan(coords, minPts = min_cluster_size)
            raw_clusters <- hdb_res$cluster

            unique_real <- unique(raw_clusters[raw_clusters > 0])
            n_subclusters <- length(unique_real)

            if (n_subclusters == 0) {
                # All noise — treat as single cluster
                cluster_ids <- rep(1L, n_cls)
                n_subclusters <- 1L
            } else if (n_subclusters == 1) {
                # Single cluster found — assign all (including noise) to it
                cluster_ids <- rep(1L, n_cls)
            } else {
                # Multiple subclusters found — reassign noise to nearest centroid
                cluster_ids <- raw_clusters

                # Compute centroids of real clusters
                centroids <- do.call(rbind, lapply(unique_real, function(cid) {
                    colMeans(coords[cluster_ids == cid, , drop = FALSE])
                }))
                rownames(centroids) <- as.character(unique_real)

                # Reassign noise points
                noise_idx <- which(cluster_ids == 0)
                if (length(noise_idx) > 0) {
                    for (ni in noise_idx) {
                        dists_to_centroids <- sqrt(
                            (centroids[, 1] - coords[ni, 1])^2 +
                            (centroids[, 2] - coords[ni, 2])^2
                        )
                        cluster_ids[ni] <- unique_real[which.min(dists_to_centroids)]
                    }
                }

                # Re-index clusters sequentially (1, 2, 3, ...)
                # Order by cluster size (largest first → A, B, C)
                tab <- sort(table(cluster_ids), decreasing = TRUE)
                old_ids <- as.integer(names(tab))
                new_ids <- seq_along(old_ids)
                remap <- stats::setNames(new_ids, old_ids)
                cluster_ids <- as.integer(remap[as.character(cluster_ids)])
                n_subclusters <- length(new_ids)
            }
        }

        # Build labels
        if (n_subclusters == 1) {
            sub_labels <- rep(cls, n_cls)
        } else {
            # Alphabetical suffix: A = largest, B = second, ...
            sub_labels <- paste0(cls, " [", LETTERS[cluster_ids], "]")
        }

        # Compute centroids and sizes for each subcluster
        sub_info <- data.table::data.table(
            Sample_ID_Full = ref$Sample_ID_Full[idx],
            UMAP1 = coords[, 1],
            UMAP2 = coords[, 2],
            Class_Label = cls,
            Subcluster_ID = cluster_ids,
            Subcluster_Label = sub_labels
        )

        # Add centroid and size per subcluster
        centroids_dt <- sub_info[, .(
            Centroid_UMAP1 = mean(UMAP1),
            Centroid_UMAP2 = mean(UMAP2),
            Subcluster_N = .N
        ), by = Subcluster_ID]

        sub_info <- merge(sub_info, centroids_dt, by = "Subcluster_ID", all.x = TRUE)

        results_list[[i]] <- sub_info
    }

    final_dt <- data.table::rbindlist(results_list)

    # Summary
    n_total_subs <- length(unique(paste0(final_dt$Class_Label, "_", final_dt$Subcluster_ID)))
    n_multi <- sum(sapply(classes, function(cls) {
        length(unique(final_dt[Class_Label == cls]$Subcluster_ID)) > 1
    }))

    message("\nSubclustering complete:")
    message("  Total classes: ", length(classes))
    message("  Classes with multiple subclusters: ", n_multi)
    message("  Total subclusters: ", n_total_subs)

    return(final_dt)
}


#' Compute Subcluster-Based Distances
#'
#' For each target (projected) sample, computes the mean Euclidean distance to every
#' reference subcluster (all members, not KNN). Returns a ranked table of closest
#' subclusters per target sample.
#'
#' This is the subcluster-aware alternative to \code{get_top_3_distances()}.
#'
#' @param metadata_input Character path or data.frame/data.table. The combined metadata.
#' @param subcluster_table data.table. Output of \code{subcluster_reference_classes()}.
#'   If NULL, subclustering is computed on the fly with default parameters.
#' @param sample_ids Character vector. Optional IDs of target samples. If NULL, uses
#'   samples with Batch == "New_Projected".
#' @param n_top Integer. Number of top closest subclusters to report per sample.
#'   Default is 5.
#' @param min_cluster_size Integer. Passed to \code{subcluster_reference_classes()} if
#'   \code{subcluster_table} is NULL. Default is 3.
#' @param collapse_to_class Logical. If TRUE, after computing per-subcluster distances,
#'   only the closest subcluster per class is kept (so you get one entry per class).
#'   Default is TRUE.
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{Sample_ID}{Character. Target sample identifier.}
#'     \item{Rank}{Integer. Rank of proximity (1 = closest).}
#'     \item{Subcluster_Label}{Character. Label of the subcluster.}
#'     \item{Class_Label}{Character. Parent class name.}
#'     \item{Mean_Distance}{Numeric. Mean Euclidean distance to all subcluster members.}
#'     \item{Subcluster_N}{Integer. Number of samples in the subcluster.}
#'     \item{Centroid_UMAP1, Centroid_UMAP2}{Numeric. Subcluster centroid.}
#'   }
#'
#' @export
#' @family subcluster
get_subcluster_distances <- function(metadata_input,
                                     subcluster_table = NULL,
                                     sample_ids = NULL,
                                     n_top = 5,
                                     min_cluster_size = 3,
                                     collapse_to_class = TRUE) {

    # Load metadata
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found: ", metadata_input)
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }

    # Get target samples
    split_data <- .split_ref_target(meta, sample_ids)
    new_samples <- split_data$new

    if (!is.null(sample_ids)) {
        # Also pick up samples by grep (consistent with existing pipeline)
        indices <- unique(unlist(lapply(sample_ids, function(x) grep(x, meta$Sample_ID_Full))))
        if (length(indices) > 0) {
            extra <- meta[indices]
            new_samples <- data.table::rbindlist(list(new_samples, extra), use.names = TRUE, fill = TRUE)
            new_samples <- unique(new_samples, by = "Sample_ID_Full")
        }
    }

    if (nrow(new_samples) == 0) {
        warning("No target samples found.")
        return(NULL)
    }

    # Compute subclusters if not provided
    if (is.null(subcluster_table)) {
        message("No precomputed subcluster table provided. Computing now...")
        subcluster_table <- subcluster_reference_classes(
            metadata_input = metadata_input,
            sample_ids = if (!is.null(sample_ids)) sample_ids else new_samples$Sample_ID_Full,
            min_cluster_size = min_cluster_size
        )
    }

    # Validate subcluster table
    req_cols <- c("Sample_ID_Full", "UMAP1", "UMAP2", "Class_Label",
                  "Subcluster_ID", "Subcluster_Label", "Subcluster_N",
                  "Centroid_UMAP1", "Centroid_UMAP2")
    missing <- setdiff(req_cols, names(subcluster_table))
    if (length(missing) > 0) {
        stop("Subcluster table is missing columns: ", paste(missing, collapse = ", "))
    }

    # Precompute unique subclusters for efficiency
    unique_subs <- unique(subcluster_table[, .(Class_Label, Subcluster_ID, Subcluster_Label,
                                                Subcluster_N, Centroid_UMAP1, Centroid_UMAP2)])

    message("Computing distances for ", nrow(new_samples), " target samples to ",
            nrow(unique_subs), " subclusters...")

    results_list <- vector("list", nrow(new_samples))

    for (i in seq_len(nrow(new_samples))) {
        target <- new_samples[i, ]
        t_u1 <- target$UMAP1
        t_u2 <- target$UMAP2
        sid <- target$Sample_ID_Full

        # Calculate distance to every reference sample in the subcluster table
        subcluster_table[, dist_to_target := sqrt((UMAP1 - t_u1)^2 + (UMAP2 - t_u2)^2)]

        # Aggregate: mean distance per subcluster (ALL members, not KNN)
        sub_dists <- subcluster_table[, .(
            Mean_Distance = mean(dist_to_target, na.rm = TRUE),
            SD_Distance   = stats::sd(dist_to_target, na.rm = TRUE),
            Subcluster_N  = .N
        ), by = .(Class_Label, Subcluster_ID, Subcluster_Label,
                  Centroid_UMAP1, Centroid_UMAP2)]

        # Collapse: keep only the closest subcluster per class
        if (collapse_to_class) {
            sub_dists <- sub_dists[sub_dists[, .I[which.min(Mean_Distance)], by = Class_Label]$V1]
        }

        # Rank and take top N
        sub_dists <- sub_dists[order(Mean_Distance)]
        sub_dists <- sub_dists[seq_len(min(n_top, nrow(sub_dists)))]
        sub_dists[, Rank := seq_len(.N)]
        sub_dists[, Sample_ID := sid]

        results_list[[i]] <- sub_dists
    }

    # Cleanup temporary column
    subcluster_table[, dist_to_target := NULL]

    final_dt <- data.table::rbindlist(results_list)

    # Reorder columns for readability
    col_order <- c("Sample_ID", "Rank", "Subcluster_Label", "Class_Label",
                   "Mean_Distance", "SD_Distance", "Subcluster_N",
                   "Subcluster_ID", "Centroid_UMAP1", "Centroid_UMAP2")
    col_order <- intersect(col_order, names(final_dt))
    final_dt <- final_dt[, ..col_order]
    final_dt[, Mean_Distance := round(Mean_Distance, 4)]
    final_dt[, SD_Distance := round(SD_Distance, 4)]

    message("Subcluster distance computation complete.")
    return(final_dt)
}


#' Get Top 3 Subcluster Distances (Wide Format)
#'
#' Convenience wrapper around \code{get_subcluster_distances()} that returns
#' results in wide format, analogous to \code{get_top_3_distances()}.
#'
#' @inheritParams get_subcluster_distances
#' @return A \code{data.table} with one row per target sample and columns:
#'   Sample_ID, Group1_Name, Group1_Dist, Group1_N, Group2_Name, ...
#'
#' @export
#' @family subcluster
get_top_3_subcluster_distances <- function(metadata_input,
                                            subcluster_table = NULL,
                                            sample_ids = NULL,
                                            min_cluster_size = 3) {

    long_dt <- get_subcluster_distances(
        metadata_input = metadata_input,
        subcluster_table = subcluster_table,
        sample_ids = sample_ids,
        n_top = 3,
        min_cluster_size = min_cluster_size,
        collapse_to_class = TRUE
    )

    if (is.null(long_dt) || nrow(long_dt) == 0) return(NULL)

    # Pivot to wide format
    results_list <- lapply(unique(long_dt$Sample_ID), function(sid) {
        sub <- long_dt[Sample_ID == sid]
        # Pad to 3 rows if fewer
        while (nrow(sub) < 3) {
            sub <- rbind(sub, data.table::data.table(
                Sample_ID = sid, Rank = nrow(sub) + 1,
                Subcluster_Label = NA_character_, Class_Label = NA_character_,
                Mean_Distance = NA_real_, SD_Distance = NA_real_,
                Subcluster_N = NA_integer_
            ), fill = TRUE)
        }
        data.table::data.table(
            Sample_ID    = sid,
            Group1_Name  = sub$Subcluster_Label[1],
            Group1_Dist  = sub$Mean_Distance[1],
            Group1_N     = sub$Subcluster_N[1],
            Group2_Name  = sub$Subcluster_Label[2],
            Group2_Dist  = sub$Mean_Distance[2],
            Group2_N     = sub$Subcluster_N[2],
            Group3_Name  = sub$Subcluster_Label[3],
            Group3_Dist  = sub$Mean_Distance[3],
            Group3_N     = sub$Subcluster_N[3]
        )
    })

    data.table::rbindlist(results_list)
}


#' Plot UMAP Zoom with Subcluster Distances (Single Sample)
#'
#' Generates a focused UMAP plot for a specific target sample, using
#' subcluster-based distance computation instead of raw KNN. The bar plot
#' shows mean distance to each subcluster (closest subcluster per class),
#' and the UMAP plot highlights the local neighborhood.
#'
#' Drop-in replacement for \code{plot_umap_zoom()} — same signature plus
#' subcluster-specific parameters.
#'
#' @param metadata Character or data.frame. Path to metadata or the object itself.
#' @param sample_id Character. ID of the target sample.
#' @param subcluster_table data.table. Precomputed subclusters from
#'   \code{subcluster_reference_classes()}. If NULL, computed on the fly.
#' @param cases_files Character vector. Optional list of other projected IDs to hide.
#' @param fixed.colors Named character vector. Custom color palette.
#' @param circle_radius Numeric. Circle radius on the UMAP. Default is 1.5.
#' @param distance_multiplier Numeric. Multiplier for zoom area. Default is 2.
#' @param n_top_groups Integer. Top groups in bar plot. Default is 10.
#' @param min_cluster_size Integer. For on-the-fly subclustering. Default is 3.
#' @param collapse_to_class Logical. Show only closest subcluster per class in
#'   bar plot. Default is TRUE.
#' @param show_subcluster_boundaries Logical. Draw convex hulls around subclusters
#'   in the UMAP. Default is FALSE.
#'
#' @return A patchwork ggplot object.
#' @export
#' @family subcluster
plot_umap_zoom_subclusters <- function(metadata,
                                        sample_id,
                                        subcluster_table = NULL,
                                        cases_files = NULL,
                                        fixed.colors = NULL,
                                        circle_radius = 1.5,
                                        distance_multiplier = 2,
                                        n_top_groups = 10,
                                        min_cluster_size = 3,
                                        collapse_to_class = TRUE,
                                        show_subcluster_boundaries = FALSE) {

    # --- 1. Data Preparation ---
    if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata)

    target.idx_all <- grep(sample_id, metadata$Sample_ID_Full)
    if (length(target.idx_all) == 0) {
        warning(paste("ERROR: Sample '", sample_id, "' not found!"))
        return(NULL)
    }
    target.idx_orig <- target.idx_all[1]

    # Filter out other projected cases
    if (!is.null(cases_files) && length(cases_files) > 0) {
        all_case_indices <- unique(unlist(lapply(cases_files, function(x) grep(x, metadata$Sample_ID_Full))))
        indices_to_remove <- setdiff(all_case_indices, target.idx_orig)
        if (length(indices_to_remove) > 0) {
            meta_proc <- metadata[-indices_to_remove, ]
        } else {
            meta_proc <- metadata
        }
    } else {
        meta_proc <- metadata
    }

    target.idx.new <- grep(sample_id, meta_proc$Sample_ID_Full)[1]
    target.row <- meta_proc[target.idx.new, ]
    target.umap1 <- target.row$UMAP1
    target.umap2 <- target.row$UMAP2

    # Handle Class Label naming
    if (!"sample.title" %in% names(meta_proc)) {
        if ("general_class" %in% names(meta_proc)) meta_proc$sample.title <- meta_proc$general_class
        else if ("methylation_class" %in% names(meta_proc)) meta_proc$sample.title <- meta_proc$methylation_class
        else meta_proc$sample.title <- "Other"
    }
    meta_proc$sample.title <- as.character(meta_proc$sample.title)
    meta_proc$sample.title[target.idx.new] <- "Target Sample"

    # --- 2. Compute Subcluster Distances ---
    # Build subcluster table if not provided
    if (is.null(subcluster_table)) {
        subcluster_table <- subcluster_reference_classes(
            metadata_input = meta_proc,
            sample_ids = target.row$Sample_ID_Full,
            min_cluster_size = min_cluster_size
        )
    }

    # Compute distances from target to subclusters
    sub_dt <- data.table::as.data.table(subcluster_table)
    sub_dt[, dist_to_target := sqrt((UMAP1 - target.umap1)^2 + (UMAP2 - target.umap2)^2)]

    distance_stats <- sub_dt[, .(
        mean_distance = mean(dist_to_target, na.rm = TRUE),
        sd_distance   = stats::sd(dist_to_target, na.rm = TRUE),
        n_samples     = .N
    ), by = .(Class_Label, Subcluster_ID, Subcluster_Label)]

    # Collapse to closest subcluster per class if requested
    if (collapse_to_class) {
        distance_stats <- distance_stats[distance_stats[, .I[which.min(mean_distance)], by = Class_Label]$V1]
        # Use Subcluster_Label as the display name
        distance_stats[, display_label := Subcluster_Label]
    } else {
        distance_stats[, display_label := Subcluster_Label]
    }

    distance_stats <- distance_stats[order(mean_distance)]
    top_groups <- distance_stats[seq_len(min(n_top_groups, nrow(distance_stats)))]
    top_groups[, is_first := ifelse(seq_len(.N) == 1, "first", "other")]

    # Clean up
    sub_dt[, dist_to_target := NULL]

    # --- 3. Define Zoom Area ---
    meta_proc$dist_to_target <- sqrt((meta_proc$UMAP1 - target.umap1)^2 +
                                     (meta_proc$UMAP2 - target.umap2)^2)
    max_distance_zoom <- circle_radius * distance_multiplier
    meta_zoom <- meta_proc[meta_proc$dist_to_target <= max_distance_zoom, ]

    circle.df <- data.frame(
        x = target.umap1 + circle_radius * cos(seq(0, 2 * pi, length.out = 100)),
        y = target.umap2 + circle_radius * sin(seq(0, 2 * pi, length.out = 100))
    )

    # --- 4. Color Management ---
    classes_umap <- unique(meta_zoom$sample.title[!is.na(meta_zoom$sample.title)])
    classes_dist <- unique(top_groups$display_label[!is.na(top_groups$display_label)])
    all_classes <- unique(c(classes_umap, classes_dist))

    hardcoded_colors <- c(
        "Target Sample"  = "#D55E00",
        "ETMR"           = "#FF0000",
        "MB, G4"         = "#38b6ff",
        "MB, G3"         = "#005383",
        "MB, SHH INF"    = "#048f6a",
        "MB, SHH CHL AD" = "#36f2c0",
        "MB, WNT"        = "#47006d",
        "DLGNT"          = "#F032E6"
    )

    use_colors <- if (is.null(fixed.colors)) hardcoded_colors else fixed.colors

    missing_classes <- setdiff(all_classes, names(use_colors))
    if (length(missing_classes) > 0) {
        generated_colors <- scales::hue_pal()(length(missing_classes))
        names(generated_colors) <- missing_classes
        current_palette <- c(use_colors, generated_colors)
    } else {
        current_palette <- use_colors
    }

    border_palette <- current_palette
    if ("Target Sample" %in% names(border_palette)) border_palette["Target Sample"] <- "black"

    # --- 5. Bar Plot (Subcluster Distances) ---
    # Prepare jitter data: reference samples belonging to the top groups' subclusters
    jitter_data <- sub_dt[paste0(Class_Label, "_", Subcluster_ID) %in%
                          top_groups[, paste0(Class_Label, "_", Subcluster_ID)]]
    jitter_data[, dist_to_target := sqrt((UMAP1 - target.umap1)^2 + (UMAP2 - target.umap2)^2)]

    # Map jitter data to display labels
    label_map <- top_groups[, .(Class_Label, Subcluster_ID, display_label)]
    jitter_data <- merge(jitter_data, label_map, by = c("Class_Label", "Subcluster_ID"), all.x = TRUE)

    p_bar_dist <- ggplot2::ggplot(
        top_groups,
        ggplot2::aes(x = stats::reorder(display_label, -mean_distance), y = mean_distance)
    ) +
        ggplot2::geom_bar(
            ggplot2::aes(fill = is_first), stat = "identity",
            alpha = 0.9, width = 0.8, show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = c("first" = "#1F4E79", "other" = "grey80")) +
        ggplot2::geom_jitter(
            data = jitter_data,
            ggplot2::aes(x = display_label, y = dist_to_target),
            color = "grey40", size = 1.5, width = 0.15, alpha = 0.5, show.legend = FALSE
        ) +
        ggplot2::geom_text(ggplot2::aes(
            y = mean_distance + 0.05 * max(top_groups$mean_distance),
            label = paste0("n=", n_samples)
        ), hjust = 0, size = 5, fontface = "bold") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 16) +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor   = ggplot2::element_blank(),
            axis.title.y       = ggplot2::element_blank(),
            axis.text.y        = ggplot2::element_text(size = 12, face = "bold"),
            axis.text.x        = ggplot2::element_text(size = 12),
            plot.title          = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
            plot.margin         = ggplot2::margin(b = 10)
        ) +
        ggplot2::labs(
            y = "Mean Distance (Subcluster)",
            x = "",
            title = "Closest Groups (Subcluster)"
        )

    # --- 6. UMAP Plot ---
    meta_zoom$label_display <- NA
    target_point_data <- meta_zoom[meta_zoom$sample.title == "Target Sample", ]

    p_umap <- ggplot2::ggplot(meta_zoom, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(
            ggplot2::aes(fill = sample.title, color = sample.title),
            shape = 21, size = 4, alpha = 0.7, stroke = 0.8
        ) +
        ggplot2::scale_fill_manual(values = current_palette, na.value = "grey85", drop = TRUE) +
        ggplot2::scale_color_manual(values = border_palette, na.value = "grey85", drop = TRUE) +
        ggplot2::geom_path(
            data = circle.df, ggplot2::aes(x = x, y = y),
            color = "red", linewidth = 1, linetype = "dashed", inherit.aes = FALSE
        ) +
        ggplot2::geom_point(
            data = target_point_data,
            ggplot2::aes(x = UMAP1, y = UMAP2),
            fill = current_palette["Target Sample"], color = "black",
            shape = 21, size = 3.8, stroke = 0.8, inherit.aes = FALSE
        ) +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw(base_size = 16) +
        ggplot2::theme(
            legend.text  = ggplot2::element_text(size = 14),
            legend.title = ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank(),
            plot.title   = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)
        ) +
        ggplot2::labs(x = NULL, y = NULL, title = "Local Neighborhood")

    # Optional: draw convex hulls for subclusters visible in the zoom area
    if (show_subcluster_boundaries) {
        sub_in_zoom <- sub_dt[UMAP1 >= min(meta_zoom$UMAP1) & UMAP1 <= max(meta_zoom$UMAP1) &
                              UMAP2 >= min(meta_zoom$UMAP2) & UMAP2 <= max(meta_zoom$UMAP2)]
        if (nrow(sub_in_zoom) > 0) {
            # Compute convex hulls per subcluster
            hulls <- sub_in_zoom[, {
                if (.N >= 3) {
                    hull_idx <- grDevices::chull(UMAP1, UMAP2)
                    hull_idx <- c(hull_idx, hull_idx[1])  # close polygon
                    .(UMAP1 = UMAP1[hull_idx], UMAP2 = UMAP2[hull_idx])
                } else {
                    NULL
                }
            }, by = .(Subcluster_Label)]

            if (nrow(hulls) > 0) {
                p_umap <- p_umap +
                    ggplot2::geom_polygon(
                        data = hulls,
                        ggplot2::aes(x = UMAP1, y = UMAP2, group = Subcluster_Label),
                        fill = NA, color = "grey30", linetype = "dotted",
                        linewidth = 0.6, alpha = 0.8, inherit.aes = FALSE
                    )
            }
        }
    }

    # --- 7. Combine ---
    n_legend_items <- length(unique(meta_zoom$sample.title[!is.na(meta_zoom$sample.title)]))
    legend_rows <- ifelse(n_legend_items > 6, 2, 1)

    combined_plot <- patchwork::wrap_plots(p_umap, p_bar_dist, ncol = 2, widths = c(1.5, 1)) +
        patchwork::plot_layout(guides = 'collect') &
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = legend_rows, override.aes = list(size = 5, alpha = 1)),
            fill  = ggplot2::guide_legend(nrow = legend_rows, override.aes = list(size = 5, alpha = 1))
        ) &
        ggplot2::theme(
            legend.position  = "bottom",
            legend.direction = "horizontal",
            legend.box       = "horizontal",
            legend.margin    = ggplot2::margin(t = 10),
            legend.text      = ggplot2::element_text(size = 14),
            legend.key.size  = ggplot2::unit(0.8, "cm")
        )

    return(combined_plot)
}


#' Plot UMAP Zoom Subclusters Combined (Batch Wrapper)
#'
#' Iterates over multiple projected samples, generates subcluster-based zoom plots
#' for each, and saves the result to a PDF. Analogous to \code{plot_umap_zoom_combined()}.
#'
#' Precomputes subclusters once and reuses them across all samples for efficiency.
#'
#' @param metadata.in Character or data.frame. The combined metadata.
#' @param output.path Character. Output directory.
#' @param name Character. Base name for the output file.
#' @param sample_id Character vector. Specific IDs to plot. If NULL, uses "New_Projected".
#' @param ncol Integer. Grid columns. Default is 2.
#' @param width Numeric or "Auto". Plot width.
#' @param height Numeric or "Auto". Plot height.
#' @param n_top_groups Integer. Top groups in bar plot. Default is 10.
#' @param min_cluster_size Integer. HDBSCAN minPts. Default is 3.
#' @param collapse_to_class Logical. Show only closest subcluster per class. Default TRUE.
#' @param show_subcluster_boundaries Logical. Draw convex hulls. Default FALSE.
#'
#' @return NULL (saves PDF to disk).
#' @export
#' @family subcluster
plot_umap_zoom_subclusters_combined <- function(metadata.in, output.path, name,
                                                 sample_id = NULL,
                                                 ncol = 2,
                                                 width = "Auto",
                                                 height = "Auto",
                                                 n_top_groups = 10,
                                                 min_cluster_size = 3,
                                                 collapse_to_class = TRUE,
                                                 show_subcluster_boundaries = FALSE) {

    # Load metadata
    if (is.character(metadata.in)) {
        if (!file.exists(metadata.in)) stop("Metadata file not found.")
        metadata <- data.table::fread(metadata.in)
    } else {
        metadata <- data.table::as.data.table(metadata.in)
    }

    if (!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)
    output.file <- file.path(output.path, paste0(name, ".subcluster.results.plots.pdf"))

    # Identify cases
    if (is.null(sample_id)) {
        if ("Batch" %in% names(metadata)) {
            cases_to_plot <- metadata$Sample_ID_Full[metadata$Batch == "New_Projected"]
        } else {
            stop("Cannot identify new samples (missing 'Batch' column).")
        }
    } else {
        cases_to_plot <- sample_id
    }

    if (length(cases_to_plot) == 0) {
        warning("No cases to plot.")
        return(invisible(NULL))
    }

    # PRECOMPUTE subclusters ONCE (key optimization)
    message("Precomputing subclusters for all reference classes...")
    subcluster_table <- subcluster_reference_classes(
        metadata_input = metadata,
        sample_ids = cases_to_plot,
        min_cluster_size = min_cluster_size
    )

    # Generate plots
    plot_list <- vector("list", length(cases_to_plot))
    for (k in seq_along(cases_to_plot)) {
        sid <- cases_to_plot[k]
        message(paste("Processing Subcluster Plot:", sid))

        p <- tryCatch({
            plot_umap_zoom_subclusters(
                metadata = metadata,
                sample_id = sid,
                subcluster_table = subcluster_table,
                cases_files = cases_to_plot,
                circle_radius = 0.3,
                distance_multiplier = 1,
                n_top_groups = n_top_groups,
                min_cluster_size = min_cluster_size,
                collapse_to_class = collapse_to_class,
                show_subcluster_boundaries = show_subcluster_boundaries,
                fixed.colors = NULL
            )
        }, error = function(e) {
            warning("Error plotting ", sid, ": ", e$message)
            NULL
        })

        if (!is.null(p)) {
            p_header <- ggplot2::ggplot() +
                ggplot2::annotate("text", x = 0.5, y = 0.5,
                                  label = paste0("Case ", k, ": ", sid),
                                  size = 6, fontface = "bold") +
                ggplot2::theme_void() +
                ggplot2::theme(plot.margin = ggplot2::margin(b = 5, t = 10))

            plot_list[[k]] <- patchwork::wrap_plots(p_header, p, ncol = 1, heights = c(1, 15))
        } else {
            plot_list[[k]] <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(ggplot2::aes(0, 0, label = paste("No Data:", sid)))
        }
    }

    plot_list <- plot_list[!sapply(plot_list, is.null)]

    if (height == "Auto") height <- max(9, 5 * length(cases_to_plot) / 2)
    if (width == "Auto") width <- max(16, min(ncol, length(cases_to_plot)) * 10)

    if (length(plot_list) > 0) {
        final_grid <- patchwork::wrap_plots(plot_list, ncol = ncol)
        ggplot2::ggsave(output.file, final_grid, width = width, height = height, limitsize = FALSE)
        message("Subcluster plots saved to: ", output.file)
    } else {
        warning("No plots generated.")
    }

    invisible(NULL)
}


#' Visualize Subclusters on Reference Atlas
#'
#' Diagnostic plot showing all subclusters identified across the reference UMAP,
#' useful for inspecting the quality of subclustering before running analyses.
#'
#' @param subcluster_table data.table. Output of \code{subcluster_reference_classes()}.
#' @param highlight_class Character. Optional: if specified, only shows subclusters
#'   for this class (others dimmed).
#' @param output_path Character. Optional PDF output path.
#'
#' @return A ggplot object.
#' @export
#' @family subcluster
plot_subcluster_atlas <- function(subcluster_table,
                                   highlight_class = NULL,
                                   output_path = NULL) {

    dt <- data.table::as.data.table(subcluster_table)

    if (!is.null(highlight_class)) {
        dt[, PlotGroup := ifelse(Class_Label == highlight_class, Subcluster_Label, "Other")]
        dt[, PlotAlpha := ifelse(Class_Label == highlight_class, 0.8, 0.15)]
        dt[, PlotSize := ifelse(Class_Label == highlight_class, 2, 0.5)]

        # Color only highlighted class
        highlight_labels <- unique(dt[Class_Label == highlight_class]$Subcluster_Label)
        n_subs <- length(highlight_labels)
        sub_colors <- c(scales::brewer_pal(palette = "Set1")(min(n_subs, 9)))
        if (n_subs > 9) sub_colors <- c(sub_colors, scales::hue_pal()(n_subs - 9))
        names(sub_colors) <- highlight_labels[seq_along(sub_colors)]
        sub_colors["Other"] <- "grey85"

        title_text <- paste("Subclusters of:", highlight_class)
    } else {
        # Show all classes, color by subcluster label
        dt[, PlotGroup := Subcluster_Label]
        dt[, PlotAlpha := 0.6]
        dt[, PlotSize := 1]
        sub_colors <- NULL
        title_text <- "Subcluster Atlas (All Classes)"
    }

    # Order so highlighted points render on top
    dt <- dt[order(PlotAlpha)]

    p <- ggplot2::ggplot(dt, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(ggplot2::aes(color = PlotGroup),
                            size = dt$PlotSize, alpha = dt$PlotAlpha) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = "Subcluster") +
        ggplot2::coord_fixed()

    if (!is.null(sub_colors)) {
        p <- p + ggplot2::scale_color_manual(values = sub_colors, na.value = "grey85")
    }

    # Hide legend if too many groups
    n_groups <- length(unique(dt$PlotGroup))
    if (n_groups > 30) {
        p <- p + ggplot2::theme(legend.position = "none")
    }

    if (!is.null(output_path)) {
        ggplot2::ggsave(output_path, p, width = 14, height = 10, limitsize = FALSE)
        message("Subcluster atlas saved to: ", output_path)
    }

    return(p)
}