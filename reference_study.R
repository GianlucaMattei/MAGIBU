#' Analyze Reference Class Neighbors (Distance Based)
#' 
#' For each sample of a specific reference class (e.g., DLGNT), this function calculates 
#' the mean Local KNN distance to all other classes. It aggregates results to show 
#' which classes are topologically closest on average.
#' 
#' @param metadata_input Character or data.frame. Path to metadata file or the metadata object itself.
#' @param target_class_name Character. Exact name of the class to analyze (e.g., "DLGNT").
#' @param knn_filter Integer. Number of neighbors for local distance calculation (default 5).
#' 
#' @return A tibble with global average distances and standard deviations for each neighbor class.
#' @importFrom dplyr %>%
#' @export
analyze_reference_class_neighbors <- function(metadata_input, target_class_name, knn_filter = 5) {
    
    # 1. Data Preparation
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }
    
    # Identify correct class column
    class_cols <- c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")
    class_col <- intersect(class_cols, names(meta))[1]
    if (is.na(class_col)) stop("Unable to find a valid class column in metadata.")
    
    # Standardize class column name
    meta$Class_Label <- meta[[class_col]]
    
    # 2. Select Target Samples
    target_samples <- meta[meta$Class_Label == target_class_name, ]
    
    if (nrow(target_samples) == 0) {
        stop(paste("No samples found for class:", target_class_name))
    }
    
    # Calculate total samples per class for final report context
    class_counts <- meta %>% 
        dplyr::group_by(Class_Label) %>% 
        dplyr::summarise(Total_N = dplyr::n())
    
    message(paste("Analyzing distances for", nrow(target_samples), "samples of class", target_class_name, "..."))
    
    # 3. Iterate over target samples
    results_list <- lapply(1:nrow(target_samples), function(i) {
        current_sample <- target_samples[i,]
        current_id <- current_sample$Sample_ID_Full
        
        # Reference Set = All samples except self (Leave-One-Out)
        ref_set <- meta[meta$Sample_ID_Full != current_id, ]
        
        # Calculate Euclidean Distance to all other points
        dists <- sqrt((ref_set$UMAP1 - current_sample$UMAP1)^2 + 
                      (ref_set$UMAP2 - current_sample$UMAP2)^2)
        
        # Temporary data.table for sorting
        temp_df <- data.table::data.table(
            Neighbor_Class = ref_set$Class_Label,
            Dist = dists
        )
        
        # Handle infinite/large knn_filter
        limit_n <- if(is.infinite(knn_filter)) nrow(temp_df) else knn_filter

        # Calculate local mean distance (KNN) for each neighbor class
        class_stats <- temp_df %>%
            dplyr::group_by(Neighbor_Class) %>%
            dplyr::summarise(
                Local_Mean_Dist = mean(sort(Dist)[1:min(dplyr::n(), limit_n)], na.rm = TRUE),
                .groups = 'drop'
            )
        
        return(class_stats)
    })
    
    # 4. Global Aggregation
    all_results <- dplyr::bind_rows(results_list)
    
    # Calculate mean of mean distances across all target samples
    final_stats <- all_results %>%
        dplyr::group_by(Neighbor_Class) %>%
        dplyr::summarise(
            Global_Avg_Distance = mean(Local_Mean_Dist, na.rm = TRUE),
            SD_Distance = stats::sd(Local_Mean_Dist, na.rm = TRUE)
        ) %>%
        # Add total count context
        dplyr::left_join(class_counts, by = c("Neighbor_Class" = "Class_Label")) %>%
        dplyr::arrange(Global_Avg_Distance)
    
    message("\n--- Mean Distance Analysis (Local KNN-", knn_filter, ") for ", target_class_name, " ---")
    
    return(final_stats)
}

#' Analyze Reference Class Neighbors (Detailed per Sample)
#' 
#' Similar to `analyze_reference_class_neighbors`, but returns a wide-format table
#' where each row corresponds to a single target sample, showing its closest classes 
#' and relative distances. Useful for assessing intra-class heterogeneity.
#' 
#' @param metadata_input Character or data.frame. Path to metadata file or the metadata object itself.
#' @param target_class_name Character. Exact name of the class to analyze.
#' @param knn_filter Integer. Number of neighbors for local distance calculation (default 5).
#' @param n_top_results Integer. Number of top closest classes to report per sample (default 3).
#' 
#' @return A tibble in wide format (Sample_ID, Neighbor_Class_1, Local_Mean_Dist_1, ...).
#' @importFrom dplyr %>%
#' @export
analyze_reference_class_neighbors_detailed <- function(metadata_input, target_class_name, knn_filter = 5, n_top_results = 3) {
    
    # 1. Data Preparation
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }
    
    # Check class column
    class_cols <- c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")
    class_col <- intersect(class_cols, names(meta))[1]
    if (is.na(class_col)) stop("Unable to find a valid class column in metadata.")
    meta$Class_Label <- meta[[class_col]]
    
    # 2. Select Target Samples
    target_samples <- meta[meta$Class_Label == target_class_name, ]
    if (nrow(target_samples) == 0) stop(paste("No samples found for class:", target_class_name))
    
    message(paste("Detailed analysis for", nrow(target_samples), "samples of class", target_class_name, "..."))
    
    # 3. Iterate calculations
    results_list <- lapply(1:nrow(target_samples), function(i) {
        current_sample <- target_samples[i,]
        current_id <- current_sample$Sample_ID_Full
        
        # Reference set = All except self
        ref_set <- meta[meta$Sample_ID_Full != current_id, ]
        
        dists <- sqrt((ref_set$UMAP1 - current_sample$UMAP1)^2 + 
                      (ref_set$UMAP2 - current_sample$UMAP2)^2)
        
        temp_df <- data.table::data.table(
            Neighbor_Class = ref_set$Class_Label,
            Dist = dists
        )
        
        limit_n <- if(is.infinite(knn_filter)) nrow(temp_df) else knn_filter
        
        # Calculate local mean distance per class
        class_stats <- temp_df %>%
            dplyr::group_by(Neighbor_Class) %>%
            dplyr::summarise(
                Local_Mean_Dist = round(mean(sort(Dist)[1:min(dplyr::n(), limit_n)], na.rm = TRUE), 4),
                .groups = 'drop'
            ) %>%
            dplyr::arrange(Local_Mean_Dist) %>%
            dplyr::slice_head(n = n_top_results)
        
        class_stats$Sample_ID <- current_id
        class_stats$Rank <- 1:nrow(class_stats)
        
        return(class_stats)
    })
    
    # 4. Reshape from Long to Wide
    long_df <- dplyr::bind_rows(results_list) %>%
        dplyr::select(Sample_ID, Rank, Neighbor_Class, Local_Mean_Dist)
    
    wide_df <- long_df %>%
        tidyr::pivot_wider(
            id_cols = Sample_ID,
            names_from = Rank,
            values_from = c(Neighbor_Class, Local_Mean_Dist),
            names_glue = "{.value}_{Rank}"
        )
    
    # Reorder columns to have pairs (Class_1, Dist_1, Class_2, Dist_2...)
    col_order <- c("Sample_ID")
    for(k in 1:n_top_results) {
        col_order <- c(col_order, paste0("Neighbor_Class_", k), paste0("Local_Mean_Dist_", k))
    }
    
    # Filter only existing columns (safety check)
    col_order <- intersect(col_order, names(wide_df))
    final_df <- wide_df %>% dplyr::select(dplyr::all_of(col_order))
    
    message("\n--- Analysis Complete. Returning detailed table ---")
    
    return(final_df)
}

#' Plot Reference UMAP Atlas
#' 
#' Visualizes the entire reference dataset (Atlas) colored by methylation class.
#' Useful for inspecting the general structure of the model without projected samples.
#' 
#' @param metadata_input Character or data.frame. Path to metadata file or the metadata object itself.
#' @param output_path Character (Optional). File path to save the plot as PDF.
#' 
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme labs coord_fixed ggsave
#' @export
plot_reference_umap <- function(metadata_input, output_path = NULL) {
    
    # 1. Data Preparation
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("File metadata non trovato.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }
    
    # 2. Filter Reference Samples
    # If Batch column exists, keep only training (NA or not New_Projected)
    if ("Batch" %in% names(meta)) {
        plot_data <- meta %>% dplyr::filter(is.na(Batch) | Batch != "New_Projected")
    } else {
        # If no Batch column, assume all are reference
        plot_data <- meta
    }
    
    # 3. Identify Class Column
    class_cols <- c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")
    class_col <- intersect(class_cols, names(plot_data))[1]
    
    if (is.na(class_col)) {
        warning("No class column found. Using single color.")
        plot_data$Class_Label <- "Reference"
    } else {
        plot_data$Class_Label <- plot_data[[class_col]]
    }
    
    message(paste("Generating reference map with", nrow(plot_data), "samples."))
    
    # 4. Generate Plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP1, y = UMAP2, color = Class_Label)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.6) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            legend.position = "right",
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            legend.text = ggplot2::element_text(size = 6),
            legend.key.size = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(
            title = "Reference UMAP Atlas",
            subtitle = paste("Total Samples:", nrow(plot_data), "| Class:", class_col),
            x = "UMAP 1", y = "UMAP 2",
            color = "Tumor Class"
        ) +
        ggplot2::coord_fixed()
        
    # Handle Large Legends
    n_classes <- length(unique(plot_data$Class_Label))
    if (n_classes > 40) {
        message("Note: Too many classes (", n_classes, "). Hiding legend for clarity.")
        p <- p + ggplot2::theme(legend.position = "none")
    }
    
    # 5. Output
    if (!is.null(output_path)) {
        ggplot2::ggsave(output_path, p, width = 14, height = 10, limitsize = FALSE)
        message("Reference plot saved to: ", output_path)
    }
    
    return(p)
}