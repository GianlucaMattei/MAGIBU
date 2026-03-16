#' Plot UMAP Zoom (Single Sample)
#'
#' Generates a focused UMAP plot for a specific target sample, highlighting its local neighborhood.
#' It includes a main scatter plot restricted to a specific radius and bar plots summarizing the
#' distance to the closest methylation classes.
#'
#' @param metadata Character or data.frame. Path to metadata file or the object itself.
#' @param sample_id Character. The ID of the target sample to highlight.
#' @param cases_files Character vector. Optional list of other case IDs to filter or handle specifically.
#' @param fixed.colors Named character vector. Custom color palette for classes.
#' @param circle_radius Numeric. Radius for the statistical circle calculation. Default is 1.5.
#' @param distance_multiplier Numeric. Multiplier for the zoom area radius. Default is 2.
#' @param n_top_groups Integer. Number of top groups to show in the bar plot. Default is 10.
#' @param knn_filter Integer. Number of neighbors for local distance calculation. Default is 5.
#' @param color_palette Named character vector. Optional specific color palette.
#'
#' @return A patchwork ggplot object containing the UMAP zoom and distance bar plots.
#' @importFrom dplyr %>%
#' @export
plot_umap_zoom <- function(metadata, sample_id, cases_files = NULL, fixed.colors = NULL, 
                           circle_radius = 1.5, distance_multiplier = 2, n_top_groups = 10, 
                           knn_filter = 5, color_palette = NULL) {
    
    # --- 1. Data Preparation ---
    if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata)

    target.idx_all <- grep(sample_id, metadata$Sample_ID_Full)
    if (length(target.idx_all) == 0) {
        warning(paste("ERROR: Sample '", sample_id, "' not found!"))
        return(NULL)
    }
    target.idx_orig <- target.idx_all[1]
    
    # Filter out other cases from the view if necessary
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
    
    # Calculate Euclidean distance to target
    meta_proc$dist_to_target <- sqrt((meta_proc$UMAP1 - target.umap1)^2 + 
                                     (meta_proc$UMAP2 - target.umap2)^2)
    
    meta_others <- meta_proc %>% 
        dplyr::filter(Sample_ID_Full != target.row$Sample_ID_Full & !is.na(sample.title))
    
    # --- Distance Calculation (Smart KNN) ---
    limit_n <- if (is.infinite(knn_filter)) nrow(meta_others) else knn_filter
    
    distance_stats <- meta_others %>%
        dplyr::group_by(sample.title) %>%
        dplyr::summarise(
            n_samples = dplyr::n(),
            # Mean distance of the k nearest neighbors
            mean_distance = mean(sort(dist_to_target)[1:min(dplyr::n(), limit_n)], na.rm = TRUE),
            sd_distance = stats::sd(dist_to_target, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        dplyr::arrange(mean_distance)
    
    top_groups <- distance_stats %>% dplyr::slice_head(n = n_top_groups)
    top_groups <- top_groups %>% dplyr::mutate(is_first = ifelse(dplyr::row_number() == 1, "first", "other"))
    
    # Define Zoom Area
    max_distance_zoom <- circle_radius * distance_multiplier
    meta_zoom <- meta_proc[meta_proc$dist_to_target <= max_distance_zoom, ]
    
    # Create Circle Data for plotting
    circle.df <- data.frame(
        x = target.umap1 + circle_radius * cos(seq(0, 2*pi, length.out = 100)),
        y = target.umap2 + circle_radius * sin(seq(0, 2*pi, length.out = 100))
    )

    # --- Color Management ---
    classes_umap <- unique(meta_zoom$sample.title[!is.na(meta_zoom$sample.title)])
    classes_dist <- unique(top_groups$sample.title[!is.na(top_groups$sample.title)])
    all_classes <- unique(c(classes_umap, classes_dist))
    
    hardcoded_colors <- c(
        "Target Sample"     = "#D55E00",
        "ETMR"              = "#FF0000",
        "MB, G4"            = "#38b6ff",
        "MB, G3"            = "#005383",
        "MB, SHH INF"       = "#048f6a",
        "MB, SHH CHL AD"    = "#36f2c0",
        "MB, WNT"           = "#47006d",
        "DLGNT"             = "#F032E6" # Vivid Magenta
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

    # --- Plotting ---
    y_lab_text <- if (is.infinite(knn_filter)) "Mean Euclidean Distance (Global)" else paste0("Avg Distance (Top ", knn_filter, " NN)")
    title_text <- if (is.infinite(knn_filter)) "Closest Groups (Global)" else "Closest Groups (Local)"

    # 1. Bar Plot (Distances)
    p_bar_dist <- ggplot2::ggplot(top_groups, 
                         ggplot2::aes(x = stats::reorder(sample.title, -mean_distance), y = mean_distance)) +
        ggplot2::geom_bar(ggplot2::aes(fill = is_first), stat = "identity", alpha = 0.9, width = 0.8, show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = c("first" = "#1F4E79", "other" = "grey80")) +
        ggplot2::geom_jitter(
            data = meta_others %>% dplyr::filter(sample.title %in% top_groups$sample.title),
            ggplot2::aes(x = sample.title, y = dist_to_target, color = sample.title), 
            size = 2, width = 0.15, show.legend = FALSE 
        ) +
        ggplot2::scale_color_manual(values = current_palette, drop = TRUE) + 
        # Add error bars only if global mean (infinite filter)
        {if (is.infinite(knn_filter)) ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_distance - sd_distance/sqrt(n_samples), 
                                                                          ymax = mean_distance + sd_distance/sqrt(n_samples)), 
                                                             width = 0.2, color = "red", alpha = 0.8)} +
        ggplot2::geom_text(ggplot2::aes(
            y = mean_distance + 0.05 * max(top_groups$mean_distance),
            label = paste0("n=", n_samples)
        ), hjust = 0, size = 5, fontface = "bold") + 
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 16) + 
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(), 
                       axis.title.y = ggplot2::element_blank(), 
                       axis.text.y = ggplot2::element_text(size = 14, face = "bold"), 
                       axis.text.x = ggplot2::element_text(size = 12), 
                       plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5), 
                       plot.margin = ggplot2::margin(b = 10)) +
        ggplot2::labs(y = y_lab_text, x = "", title = title_text)

    # 2. UMAP Plot
    meta_zoom$label_display <- NA
    target_point_data <- meta_zoom[meta_zoom$sample.title == "Target Sample", ]
    label_data <- meta_zoom[!is.na(meta_zoom$label_display), ]

    p_umap <- ggplot2::ggplot(meta_zoom, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point(ggplot2::aes(fill = sample.title, color = sample.title), shape = 21, size = 4, alpha = 0.7, stroke = 0.8) +
        ggplot2::scale_fill_manual(values = current_palette, na.value = "grey85", drop = TRUE) +
        ggplot2::scale_color_manual(values = border_palette, na.value = "grey85", drop = TRUE) +
        ggplot2::geom_path(data = circle.df, ggplot2::aes(x = x, y = y), color = "red", linewidth = 1, linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::geom_point(data = target_point_data, ggplot2::aes(x = UMAP1, y = UMAP2), fill = current_palette["Target Sample"], color = "black", shape = 21, size = 3.8, stroke = 0.8, inherit.aes = FALSE) +
        ggplot2::coord_fixed() + 
        ggplot2::theme_bw(base_size = 16) + 
        ggplot2::theme(legend.text = ggplot2::element_text(size = 14), 
                       legend.title = ggplot2::element_blank(), 
                       panel.grid = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)) +
        ggplot2::labs(x = NULL, y = NULL, title = "Local Neighborhood")
    
    if (nrow(label_data) > 0) {
        p_umap <- p_umap + ggrepel::geom_text_repel(data = label_data, ggplot2::aes(label = label_display), size = 5, force = 10, max.overlaps = 50, min.segment.length = 0, fontface = "bold")
    }

    n_legend_items <- length(unique(meta_zoom$sample.title[!is.na(meta_zoom$sample.title)]))
    legend_rows <- ifelse(n_legend_items > 6, 2, 1)

    combined_plot <- patchwork::wrap_plots(p_umap, p_bar_dist, ncol = 2, widths = c(1.5, 1)) + 
                     patchwork::plot_layout(guides = 'collect') & 
                     ggplot2::guides(
                         color = ggplot2::guide_legend(nrow = legend_rows, override.aes = list(size = 5, alpha = 1)),
                         fill = ggplot2::guide_legend(nrow = legend_rows, override.aes = list(size = 5, alpha = 1))
                     ) & 
                     ggplot2::theme(legend.position = "bottom", 
                           legend.direction = "horizontal",
                           legend.box = "horizontal",
                           legend.margin = ggplot2::margin(t = 10),
                           legend.text = ggplot2::element_text(size = 14), 
                           legend.key.size = ggplot2::unit(0.8, "cm"))
    
    return(combined_plot)
}

#' Plot UMAP Zoom Combined (Wrapper)
#'
#' Wrapper function that iterates over multiple samples (e.g. batch "New_Projected"),
#' generates zoom plots for each, and optionally saves them to a PDF.
#'
#' @param metadata.in Character or data.frame. Metadata object or file path.
#' @param output.path Character. Directory to save the output PDF.
#' @param name Character. Base name for the output file.
#' @param sample_id Character vector. Optional specific IDs to plot.
#' @param ncol Integer. Number of columns in the grid layout.
#' @param width Numeric/Character. Plot width (or "Auto").
#' @param height Numeric/Character. Plot height (or "Auto").
#' @param n_top_groups Integer. Number of top groups for bar plots.
#'
#' @return NULL (Saves file to disk).
#' @export
plot_umap_zoom_combined <- function(metadata.in, output.path, name, sample_id = NULL, 
                                    ncol = 2, width = "Auto", height = "Auto", n_top_groups = 10) {
    
    if (is.character(metadata.in)) { 
        if (!file.exists(metadata.in)) stop("Metadata file not found.")
        metadata <- data.table::fread(metadata.in) 
    } else { 
        metadata <- data.table::as.data.table(metadata.in) 
    }
    
    if (!dir.exists(output.path)) dir.create(output.path, recursive = TRUE)
    output.file <- file.path(output.path, paste0(name, ".results.plots.pdf"))

    if (is.null(sample_id)) {
        if ("Batch" %in% names(metadata)) { 
            cases_to_plot <- metadata$Sample_ID_Full[metadata$Batch == "New_Projected"] 
        } else { 
            stop("Cannot identify new samples (missing 'Batch' column).") 
        }
    } else { 
        cases_to_plot <- sample_id 
    }

    plot_list <- list()
    for (k in seq_along(cases_to_plot)) {
        sid <- cases_to_plot[k]
        message(paste("Processing Plot:", sid))
        p <- plot_umap_zoom(metadata, sid, cases_files = cases_to_plot, circle_radius = 0.3, distance_multiplier = 1, n_top_groups = n_top_groups, fixed.colors = NULL)
        if (!is.null(p)) {
            p_header <- ggplot2::ggplot() + 
                ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("Case ", k, ": ", sid), size = 6, fontface = "bold") + 
                ggplot2::theme_void() + 
                ggplot2::theme(plot.margin = ggplot2::margin(b = 5, t = 10))
            
            final_p <- patchwork::wrap_plots(p_header, p, ncol = 1, heights = c(1, 15))
            plot_list[[k]] <- final_p
        } else { 
            plot_list[[k]] <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::geom_text(ggplot2::aes(0, 0, label = paste("No Data:", sid))) 
        }
    }
    
    plot_list <- plot_list[!sapply(plot_list, is.null)]
    if (height == "Auto") { height <- max(9, 5 * length(cases_to_plot) / 2) }
    if (width == "Auto") { width <- max(16, min(ncol, length(cases_to_plot)) * 10) }
    
    if (length(plot_list) > 0) {
        final_grid <- patchwork::wrap_plots(plot_list, ncol = ncol)
        ggplot2::ggsave(output.file, final_grid, width = width, height = height, limitsize = FALSE)
        message(paste("Plots saved to:", output.file))
    } else { 
        warning("No plots generated.") 
    }
}

#' Get Top 3 Closest Groups
#'
#' Calculates the 3 closest tumor groups for given samples based on local neighborhood density (KNN).
#' Returns a summarized dataframe.
#'
#' @param metadata_input Character or data.frame. Metadata.
#' @param sample_ids Character vector. IDs of samples to analyze.
#' @param knn_filter Integer. Number of neighbors for local mean calculation.
#'
#' @return A data.table containing the top 3 groups and their mean distances.
#' @export
get_top_3_distances <- function(metadata_input, sample_ids = NULL, knn_filter = 5) {
    
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- data.table::as.data.table(metadata_input)
    }
    
    # Split Reference vs Target
    if ("Batch" %in% names(meta)) {
        ref_samples <- meta[Batch != "New_Projected" | is.na(Batch)]
    } else {
        if (is.null(sample_ids)) stop("'Batch' column missing. Cannot identify reference.")
        ref_samples <- meta[!Sample_ID_Full %in% sample_ids]
    }
    
    if (nrow(ref_samples) == 0) stop("CRITICAL ERROR: Reference set empty!")
    
    # Select Target
    if (is.null(sample_ids)) {
         if ("Batch" %in% names(meta)) {
            new_samples <- meta[Batch == "New_Projected"]
         } else { stop("Define sample_ids or have Batch column.") }
    } else {
        indices <- unique(unlist(lapply(sample_ids, function(x) grep(x, meta$Sample_ID_Full))))
        new_samples <- meta[indices]
    }
    
    if (nrow(new_samples) == 0) {
        warning("No new samples found.")
        return(NULL)
    }
    
    # Identify Class Column
    class_col <- NULL
    for (col in c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")) {
        if (col %in% names(ref_samples)) {
            class_col <- col
            break
        }
    }
    if (is.null(class_col)) ref_samples$Class_Label <- "Unknown" else ref_samples$Class_Label <- ref_samples[[class_col]]
    
    # Loop over targets
    results_list <- list()
    msg_knn <- if(is.infinite(knn_filter)) "Global Average" else paste0("Local KNN-", knn_filter)
    message(paste("Calculating Top 3 Groups (", msg_knn, ") for", nrow(new_samples), "samples..."))
    
    for (i in 1:nrow(new_samples)) {
        target <- new_samples[i,]
        sid <- target$Sample_ID_Full
        
        dists <- sqrt((ref_samples$UMAP1 - target$UMAP1)^2 + 
                      (ref_samples$UMAP2 - target$UMAP2)^2)
        
        temp_df <- data.table::data.table(
            Ref_Class = ref_samples$Class_Label,
            Dist = dists
        )
        
        # Aggregate by class (Local KNN logic)
        group_stats <- temp_df %>%
            dplyr::group_by(Ref_Class) %>%
            dplyr::summarise(
                Mean_Dist = mean(sort(Dist)[1:min(dplyr::n(), if(is.infinite(knn_filter)) nrow(ref_samples) else knn_filter)], na.rm = TRUE),
                N_Samples = dplyr::n()
            ) %>%
            dplyr::arrange(Mean_Dist)
        
        top3 <- group_stats %>% dplyr::slice_head(n = 3)
        
        # Fill NA if < 3
        while (nrow(top3) < 3) {
            top3 <- dplyr::bind_rows(top3, data.frame(Ref_Class=NA, Mean_Dist=NA, N_Samples=NA))
        }
        
        res_row <- data.table::data.table(
            Sample_ID = sid,
            Group1_Name = top3$Ref_Class[1], Group1_Dist = round(top3$Mean_Dist[1], 4), Group1_N = top3$N_Samples[1],
            Group2_Name = top3$Ref_Class[2], Group2_Dist = round(top3$Mean_Dist[2], 4), Group2_N = top3$N_Samples[2],
            Group3_Name = top3$Ref_Class[3], Group3_Dist = round(top3$Mean_Dist[3], 4), Group3_N = top3$N_Samples[3]
        )
        results_list[[i]] <- res_row
    }
    
    final_df <- data.table::rbindlist(results_list)
    return(final_df)
}

#' Plot Global Projection
#'
#' Generates a global UMAP map showing the reference set colored by class and 
#' new samples highlighted with numbered labels.
#'
#' @param metadata_input Character or data.frame. Metadata.
#' @param output_path Character. Optional path to save PDF.
#' @param show_labels Logical. Legacy parameter.
#' @return A ggplot object.
#' @export
plot_global_projection <- function(metadata_input, output_path = NULL, show_labels = TRUE) {
    
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("File metadata not found.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- as.data.frame(metadata_input)
    }
    
    # Robust splitting
    bg_data <- meta %>% dplyr::filter(is.na(Batch) | Batch != "New_Projected")
    fg_data <- meta %>% dplyr::filter(Batch == "New_Projected")
    
    if (nrow(fg_data) == 0) stop("No new samples (Batch == 'New_Projected') found.")
    
    # Class Column
    class_cols <- c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")
    class_col <- intersect(class_cols, names(bg_data))[1]
    
    if (is.na(class_col)) {
        warning("No class column found for reference.")
        bg_data$Ref_Class <- "Reference"
    } else {
        bg_data$Ref_Class <- bg_data[[class_col]]
    }
    
    fg_data$Plot_Num <- 1:nrow(fg_data)
    
    message("Generating Global Map...")
    
    p <- ggplot2::ggplot() +
        # Layer 1: Background
        ggplot2::geom_point(data = bg_data, 
                   ggplot2::aes(x = UMAP1, y = UMAP2, color = Ref_Class), 
                   size = 0.8, alpha = 0.5) +
        
        # Layer 2: Foreground (Dark Grey dots)
        ggplot2::geom_point(data = fg_data,
                   ggplot2::aes(x = UMAP1, y = UMAP2),
                   color = "darkgrey", size = 1) + 
        
        # Layer 3: Numbers (Black, above dot)
        ggplot2::geom_text(data = fg_data,
                  ggplot2::aes(x = UMAP1, y = UMAP2, label = Plot_Num),
                  color = "black", fontface = "bold", size = 2.5, vjust = -1) +
        
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            legend.position = "right",
            legend.title = ggplot2::element_text(size = 10, face = "bold"),
            legend.text = ggplot2::element_text(size = 6),
            legend.key.size = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(
            title = "Global Projection Overview",
            subtitle = paste0("Reference colored by ", class_col, ". New samples numbered 1-", nrow(fg_data)),
            x = "UMAP 1", y = "UMAP 2", color = "Tumor Class"
        ) +
        ggplot2::coord_fixed()
    
    if (length(unique(bg_data$Ref_Class)) > 40) {
        message("Too many reference classes. Hiding legend.")
        p <- p + ggplot2::theme(legend.position = "none")
    }
    
    if (!is.null(output_path)) {
        ggplot2::ggsave(output_path, p, width = 16, height = 12, limitsize = FALSE)
        message("Plot saved to: ", output_path)
        
        legend_df <- fg_data[, c("Plot_Num", "Sample_ID_Full")]
        names(legend_df) <- c("Map_Number", "Sample_ID")
        legend_file <- sub("\\.pdf$", "_legend_samples.csv", output_path)
        utils::write.csv(legend_df, legend_file, row.names = FALSE)
        message("Legend CSV saved to: ", legend_file)
    }
    
    return(p)
}

#' Plot UMAP Highlight Specific
#'
#' Creates a custom plot highlighting a specific target sample and a specific reference class 
#' while dimming all other samples. Supports zooming.
#'
#' @param metadata_input Metadata.
#' @param target_id Character. ID of the sample to highlight.
#' @param highlight_class Character. Name of the reference class to highlight.
#' @param output_path Character. Optional output PDF path.
#' @param zoom Logical. If TRUE, crops the plot to the relevant area.
#' @param target_point_size Numeric. Size of the target point.
#' @param highlight_point_size Numeric. Size of the highlight class points.
#' @return A ggplot object.
#' @export
plot_umap_highlight_specific <- function(metadata_input, target_id, highlight_class, output_path = NULL, 
                                         zoom = FALSE, target_point_size = 0.5, highlight_point_size = 0.5) {
    
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("File metadata not found.")
        meta <- data.table::fread(metadata_input)
    } else {
        meta <- as.data.frame(metadata_input)
    }
    
    target.idx <- grep(target_id, meta$Sample_ID_Full)
    if (length(target.idx) == 0) stop("Target sample not found.")
    target_data <- meta[target.idx, ]
    
    class_cols <- c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")
    class_col <- intersect(class_cols, names(meta))[1]
    if (is.na(class_col)) stop("No class column found.")
    meta$Ref_Class <- meta[[class_col]]
    
    meta$PlotGroup <- "Background"
    meta$PlotGroup[meta$Ref_Class == highlight_class] <- "Highlight_Class"
    meta$PlotGroup[target.idx] <- "Target"
    
    meta$PlotGroup <- factor(meta$PlotGroup, levels = c("Background", "Highlight_Class", "Target"))
    meta <- meta[order(meta$PlotGroup), ]
    
    my_colors <- c(
        "Background" = "grey85",
        "Highlight_Class" = "#F032E6", # Vivid Magenta
        "Target" = "#D55E00"           # Case 5 Orange
    )
    
    my_labels <- c("Background" = "Other Reference", "Highlight_Class" = highlight_class, "Target" = target_id)
    
    message("Generating highlight plot...")
    
    p <- ggplot2::ggplot(meta, ggplot2::aes(x = UMAP1, y = UMAP2, color = PlotGroup, size = PlotGroup, shape = PlotGroup)) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_color_manual(values = my_colors, labels = my_labels, name = "Legend") +
        ggplot2::scale_size_manual(values = c("Background" = 0.5, "Highlight_Class" = highlight_point_size, "Target" = target_point_size), guide = "none") +
        ggplot2::scale_shape_manual(values = c("Background" = 16, "Highlight_Class" = 8, "Target" = 16), guide = "none") +
        
        # Extra circle for target
        ggplot2::geom_point(data = target_data, ggplot2::aes(x = UMAP1, y = UMAP2), 
                   color = "black", shape = 21, fill = NA, size = target_point_size, stroke = 0.3, inherit.aes = FALSE) +
        
        ggrepel::geom_label_repel(data = target_data, ggplot2::aes(label = Sample_ID_Full),
                         color = "black", size = 3, fontface = "bold", 
                         box.padding = 0.5, point.padding = 0.5,
                         fill = "white", alpha = 0.8,
                         show.legend = FALSE, inherit.aes = FALSE, x = target_data$UMAP1, y = target_data$UMAP2) +
        
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::labs(title = paste("Target Placement:", target_id),
             subtitle = paste("Comparison with", highlight_class),
             x = "UMAP 1", y = "UMAP 2")
    
    if (zoom) {
        interest_points <- meta %>% dplyr::filter(PlotGroup %in% c("Target", "Highlight_Class"))
        if (nrow(interest_points) > 0) {
            min_x <- min(interest_points$UMAP1, na.rm = TRUE); max_x <- max(interest_points$UMAP1, na.rm = TRUE)
            min_y <- min(interest_points$UMAP2, na.rm = TRUE); max_y <- max(interest_points$UMAP2, na.rm = TRUE)
            
            span_x <- max_x - min_x; if(span_x == 0) span_x <- 1
            span_y <- max_y - min_y; if(span_y == 0) span_y <- 1
            
            final_xlim <- c(min_x - (span_x * 0.1), max_x + (span_x * 0.1))
            final_ylim <- c(min_y - (span_y * 0.1), max_y + (span_y * 0.1))
            
            p <- p + ggplot2::coord_fixed(xlim = final_xlim, ylim = final_ylim)
            message("Zoom applied (+10% margin).")
        } else {
            p <- p + ggplot2::coord_fixed()
        }
    } else {
        p <- p + ggplot2::coord_fixed()
    }
    
    if (!is.null(output_path)) {
        ggplot2::ggsave(output_path, p, width = 10, height = 8)
        message("Highlight plot saved: ", output_path)
    }
    
    return(p)
}

#' Analyze Reference Class Neighbors (Global Mean)
#'
#' Calculates the mean distance of a reference class to all other classes (Global Average).
#'
#' @param metadata_input Metadata.
#' @param target_class_name Class to analyze.
#' @param knn_filter Neighbors count.
#' @return Summary data.table.
#' @export
analyze_reference_class_neighbors <- function(metadata_input, target_class_name, knn_filter = 5) {
    if (is.character(metadata_input)) meta <- data.table::fread(metadata_input) else meta <- data.table::as.data.table(metadata_input)
    
    class_col <- intersect(c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass"), names(meta))[1]
    meta$Class_Label <- meta[[class_col]]
    
    target_samples <- meta[Class_Label == target_class_name]
    if (nrow(target_samples) == 0) stop("Target class not found.")
    
    class_counts <- meta %>% dplyr::group_by(Class_Label) %>% dplyr::summarise(Total_N = dplyr::n())
    message(paste("Analyzing", nrow(target_samples), "samples..."))
    
    results_list <- lapply(1:nrow(target_samples), function(i) {
        curr <- target_samples[i,]
        ref_set <- meta[Sample_ID_Full != curr$Sample_ID_Full]
        dists <- sqrt((ref_set$UMAP1 - curr$UMAP1)^2 + (ref_set$UMAP2 - curr$UMAP2)^2)
        
        limit_n <- if(is.infinite(knn_filter)) nrow(ref_set) else knn_filter
        
        data.table::data.table(Neighbor_Class = ref_set$Class_Label, Dist = dists) %>%
            dplyr::group_by(Neighbor_Class) %>%
            dplyr::summarise(Local_Mean_Dist = mean(sort(Dist)[1:min(dplyr::n(), limit_n)], na.rm = TRUE), .groups = 'drop')
    })
    
    final_stats <- dplyr::bind_rows(results_list) %>%
        dplyr::group_by(Neighbor_Class) %>%
        dplyr::summarise(Global_Avg_Distance = mean(Local_Mean_Dist, na.rm = TRUE), SD_Distance = stats::sd(Local_Mean_Dist, na.rm = TRUE)) %>%
        dplyr::left_join(class_counts, by = c("Neighbor_Class" = "Class_Label")) %>%
        dplyr::arrange(Global_Avg_Distance)
    
    return(final_stats)
}

#' Analyze Reference Class Neighbors (Detailed per Sample)
#'
#' Detailed analysis returning neighbors for each individual sample in the target class.
#'
#' @param metadata_input Metadata.
#' @param target_class_name Class to analyze.
#' @param knn_filter Neighbors count.
#' @param n_top_results Top N classes to return per sample.
#' @return Wide format data.table.
#' @export
analyze_reference_class_neighbors_detailed <- function(metadata_input, target_class_name, knn_filter = 5, n_top_results = 3) {
    if (is.character(metadata_input)) meta <- data.table::fread(metadata_input) else meta <- data.table::as.data.table(metadata_input)
    
    class_col <- intersect(c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass"), names(meta))[1]
    meta$Class_Label <- meta[[class_col]]
    
    target_samples <- meta[Class_Label == target_class_name]
    if (nrow(target_samples) == 0) stop("Target class not found.")
    
    message(paste("Detailed analysis for", nrow(target_samples), "samples..."))
    
    results_list <- lapply(1:nrow(target_samples), function(i) {
        curr <- target_samples[i,]
        ref_set <- meta[Sample_ID_Full != curr$Sample_ID_Full]
        dists <- sqrt((ref_set$UMAP1 - curr$UMAP1)^2 + (ref_set$UMAP2 - curr$UMAP2)^2)
        
        limit_n <- if(is.infinite(knn_filter)) nrow(ref_set) else knn_filter
        
        stats <- data.table::data.table(Neighbor_Class = ref_set$Class_Label, Dist = dists) %>%
            dplyr::group_by(Neighbor_Class) %>%
            dplyr::summarise(Local_Mean_Dist = round(mean(sort(Dist)[1:min(dplyr::n(), limit_n)], na.rm = TRUE), 4), .groups = 'drop') %>%
            dplyr::arrange(Local_Mean_Dist) %>%
            dplyr::slice_head(n = n_top_results)
        
        stats$Sample_ID <- curr$Sample_ID_Full
        stats$Rank <- 1:nrow(stats)
        return(stats)
    })
    
    wide_df <- dplyr::bind_rows(results_list) %>%
        tidyr::pivot_wider(id_cols = Sample_ID, names_from = Rank, values_from = c(Neighbor_Class, Local_Mean_Dist), names_glue = "{.value}_{Rank}")
    
    col_order <- c("Sample_ID", paste0(rep(c("Neighbor_Class_", "Local_Mean_Dist_"), n_top_results), rep(1:n_top_results, each=2)))
    return(wide_df %>% dplyr::select(dplyr::all_of(intersect(col_order, names(wide_df)))))
}

#' Plot Reference UMAP Atlas
#' 
#' Plots the reference map only.
#' @param metadata_input Metadata.
#' @param output_path PDF path.
#' @export
plot_reference_umap <- function(metadata_input, output_path = NULL) {
    if (is.character(metadata_input)) meta <- data.table::fread(metadata_input) else meta <- as.data.frame(metadata_input)
    
    if ("Batch" %in% names(meta)) plot_data <- meta %>% dplyr::filter(is.na(Batch) | Batch != "New_Projected") else plot_data <- meta
    
    class_col <- intersect(c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass"), names(plot_data))[1]
    if (is.na(class_col)) plot_data$Class_Label <- "Reference" else plot_data$Class_Label <- plot_data[[class_col]]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP1, y = UMAP2, color = Class_Label)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.6) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right") +
        ggplot2::labs(title = "Reference UMAP Atlas", color = "Tumor Class") +
        ggplot2::coord_fixed()
    
    if (length(unique(plot_data$Class_Label)) > 40) p <- p + ggplot2::theme(legend.position = "none")
    
    if (!is.null(output_path)) {
        ggplot2::ggsave(output_path, p, width = 14, height = 10, limitsize = FALSE)
        message("Reference plot saved.")
    }
    return(p)
}
