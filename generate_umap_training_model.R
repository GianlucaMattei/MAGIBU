#' Generate UMAP Training Model from Processed Methylation Data
#'
#' Reads processed methylation TSV files, performs variance-based feature selection,
#' PCA dimensionality reduction, and UMAP embedding. Saves a model bundle that
#' supports projection of both array and Nanopore samples via genomic positions.
#'
#' @param tsv.path Character. Directory containing processed .tsv files.
#'   Each file must have columns: \code{Probe_ID}, \code{Beta}, \code{chr_hg38}, \code{pos_hg38}.
#' @param save.path Character. Directory where model files will be saved.
#' @param metadata_file Character (Optional). Path to external CSV/TSV with clinical metadata.
#' @param name Character. Prefix for output files. Default: \code{"umap_model"}.
#' @param n_top_vars Integer. Number of most variable probes to select. Default: 50000.
#' @param pca_components Integer. Number of PCs to retain. Default: 50.
#' @param n_neighbors Integer. UMAP neighbors. Default: 30.
#' @param min_dist Numeric. UMAP minimum distance. Default: 0.3.
#' @param n_cores Integer. Cores for parallel loading. Default: 8.
#'
#' @return A list with \code{umap_model} and \code{metadata}.
#'
#' @details
#' The saved \code{.rds} model contains:
#' \describe{
#'   \item{pca_model}{Full \code{prcomp} object (supports \code{predict()}).}
#'   \item{selected_cpgs}{Probe names (e.g. \code{"cg18222083"}), ordered as in the PCA matrix.}
#'   \item{probe_to_position}{\code{data.table} mapping \code{Probe_ID} to \code{pos_id} (\code{"chr:pos"}).}
#'   \item{common_positions}{Named character vector: \code{cg* -> "chr:pos"}, ordered as \code{selected_cpgs}.
#'     This is the canonical reference for matching new samples by genomic position.}
#'   \item{col_means_training}{Named numeric vector of per-probe beta means (for imputation).}
#' }
#'
#' @importFrom foreach %dopar%
#' @export
generate_umap_training_model <- function(tsv.path, 
                                         save.path, 
                                         metadata_file = NULL,
                                         name = "umap_model",
                                         n_top_vars = 50000, 
                                         pca_components = 50, 
                                         n_neighbors = 30, 
                                         min_dist = 0.3,
                                         n_cores = 8) {
    
    if (!requireNamespace("foreach", quietly = TRUE) || !requireNamespace("doParallel", quietly = TRUE)) {
        stop("Packages 'foreach' and 'doParallel' are required.")
    }
    
    if (file.exists(save.path) && !dir.exists(save.path)) {
        stop(paste("Error: save.path '", save.path, "' exists and is a file, not a directory."))
    }
    if (!dir.exists(save.path)) dir.create(save.path, recursive = TRUE)
    
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(n_cores)
      RhpcBLASctl::omp_set_num_threads(n_cores)
    }

    # =========================================================================
    # 1. Load Data (Parallel)
    # =========================================================================
    tsv_files <- list.files(tsv.path, pattern = "\\.tsv$", full.names = TRUE)
    
    if (length(tsv_files) == 0) stop("No TSV files found in the specified path.")
    message(paste("Loading", length(tsv_files), "samples from", tsv.path, "..."))
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    `%dopar%` <- foreach::`%dopar%`
    
    data_list <- foreach::foreach(f = tsv_files, .packages = "data.table") %dopar% {
        
        header <- names(data.table::fread(f, nrows = 0))
        
        # --- Beta column detection ---
        target_col <- NULL
        if ("Beta" %in% header) {
            target_col <- "Beta"
        } else if ("beta_value" %in% header) {
            target_col <- "beta_value"
        } else {
            beta_candidates <- grep("beta", header, ignore.case = TRUE, value = TRUE)
            if (length(beta_candidates) > 0) target_col <- beta_candidates[1]
        }
        if (is.null(target_col)) {
            if (length(header) >= 2) target_col <- header[2]
            else stop(paste("Cannot identify Beta column in:", basename(f)))
        }
        
        # --- ID column detection ---
        id_col <- if ("Probe_ID" %in% header) "Probe_ID" 
                  else if ("cpg_id" %in% header) "cpg_id" 
                  else header[1]
        
        # --- Read with positions if available ---
        has_positions <- all(c("chr_hg38", "pos_hg38") %in% header)
        cols_to_read <- if (has_positions) {
            unique(c(id_col, target_col, "chr_hg38", "pos_hg38"))
        } else {
            c(id_col, target_col)
        }
        
        dt <- data.table::fread(f, select = cols_to_read)
        
        sample_name <- tools::file_path_sans_ext(basename(f))
        data.table::setnames(dt, old = id_col, new = "Probe_ID")
        data.table::setnames(dt, old = target_col, new = sample_name)
        
        return(dt)
    }
    
    message("Merging datasets...")
    
    # --- Extract position annotation from first file that has it ---
    position_lookup <- NULL
    for (dt in data_list) {
        if (all(c("chr_hg38", "pos_hg38") %in% names(dt))) {
            position_lookup <- unique(dt[!is.na(chr_hg38) & !is.na(pos_hg38), 
                                         .(Probe_ID, chr_hg38, pos_hg38)])
            position_lookup[, pos_id := paste0(chr_hg38, ":", pos_hg38)]
            message(paste("Position annotation extracted:", nrow(position_lookup), "probes with coordinates."))
            break
        }
    }
    
    if (is.null(position_lookup)) {
        warning("No position columns (chr_hg38, pos_hg38) found in training TSVs. ",
                "Nanopore projection will NOT work without external annotation.")
    }
    
    # --- Strip position columns before merge (keep only Probe_ID + beta) ---
    data_list_clean <- lapply(data_list, function(dt) {
        cols_to_drop <- intersect(c("chr_hg38", "pos_hg38"), names(dt))
        if (length(cols_to_drop) > 0) dt[, (cols_to_drop) := NULL]
        return(dt)
    })
    
    # Inner join: only probes shared across ALL samples survive
    merged_dt <- Reduce(function(x, y) merge(x, y, by = "Probe_ID", all = FALSE), data_list_clean)
    
    cpg_ids <- merged_dt$Probe_ID
    beta_matrix <- as.matrix(merged_dt[, -1, with = FALSE])
    rownames(beta_matrix) <- cpg_ids
    
    rm(data_list, data_list_clean, merged_dt)
    gc()
    
    message(paste("Matrix dimensions:", nrow(beta_matrix), "CpGs x", ncol(beta_matrix), "Samples"))
    
    # =========================================================================
    # 2. Feature Selection (Variance)
    # =========================================================================
    message("Calculating probe variance...")
    
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        probe_vars <- matrixStats::rowVars(beta_matrix)
    } else {
        probe_vars <- apply(beta_matrix, 1, stats::var)
    }
    
    top_indices <- order(probe_vars, decreasing = TRUE)[1:min(n_top_vars, length(probe_vars))]
    selected_matrix <- beta_matrix[top_indices, ]
    selected_cpgs <- rownames(selected_matrix)
    
    message(paste("Selected top", length(selected_cpgs), "variable CpGs."))
    
    # =========================================================================
    # 3. PCA
    # =========================================================================
    message("Running PCA...")
    
    t_matrix <- t(selected_matrix)
    pca_res <- stats::prcomp(t_matrix, center = TRUE, scale. = FALSE) 
    pca_data <- pca_res$x[, 1:pca_components]
    
    # =========================================================================
    # 4. Build Model Components
    # =========================================================================
    
    # Position lookup for selected CpGs
    probe_position_selected <- NULL
    common_positions <- NULL
    
    if (!is.null(position_lookup)) {
        probe_position_selected <- position_lookup[Probe_ID %in% selected_cpgs]
        
        n_with_pos <- nrow(probe_position_selected)
        n_without <- length(selected_cpgs) - n_with_pos
        if (n_without > 0) {
            message(paste("Note:", n_without, "selected probes lack genomic coordinates."))
        }
        
        # Named vector ordered like selected_cpgs: cg* -> "chr:pos"
        pos_map <- stats::setNames(probe_position_selected$pos_id, probe_position_selected$Probe_ID)
        common_positions <- pos_map[selected_cpgs]  # NA where probe has no position
        
        message(paste("Common positions:", sum(!is.na(common_positions)), "with coords,", 
                      sum(is.na(common_positions)), "without."))
    }
    
    # Column means for imputation
    col_means <- colMeans(t_matrix, na.rm = TRUE)
    
    model_components <- list(
        pca_model          = pca_res,
        pca_components     = pca_components,                 # N of PCs used for UMAP
        selected_cpgs      = selected_cpgs,
        probe_to_position  = probe_position_selected,
        common_positions   = common_positions,
        col_means_training = col_means
    )
    
    # =========================================================================
    # 5. UMAP
    # =========================================================================
    message("Running UMAP...")
    
    if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required.")
    
    set.seed(42)
    
    umap_model <- uwot::umap(
        pca_data, 
        n_neighbors = n_neighbors, 
        min_dist = min_dist,
        metric = "euclidean",
        ret_model = TRUE
    )
    
    umap_coords <- as.data.frame(umap_model$embedding)
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$Sample_ID_Full <- rownames(pca_data)
    
    # =========================================================================
    # 6. Metadata Integration
    # =========================================================================
    if (!is.null(metadata_file) && file.exists(metadata_file)) {
        message("Merging with external metadata...")
        ext_meta <- data.table::fread(metadata_file)
        final_meta <- merge(umap_coords, ext_meta, 
                           by.x = "Sample_ID_Full", by.y = colnames(ext_meta)[1], all.x = TRUE)
    } else {
        final_meta <- umap_coords
    }
    
    # =========================================================================
    # 7. Save
    # =========================================================================
    message("Saving model files...")
    
    data.table::fwrite(final_meta, 
                       file.path(save.path, paste0(name, ".umap_metadata_result.tsv")), sep = "\t")
    saveRDS(model_components, 
            file.path(save.path, paste0(name, ".umap_model.components.rds")))
    uwot::save_uwot(umap_model, 
                     file.path(save.path, paste0(name, ".umap_model.uwot")))
    
    message("=== Model generation complete ===")
    message("Model contains: pca_model, selected_cpgs, probe_to_position, common_positions, col_means_training")
    
    return(list(umap_model = umap_model, metadata = final_meta))
}