get_nearest_neighbors <- function(metadata_input, sample_ids = NULL) {
    
    # 1. Load data
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("Metadata file not found.")
        meta <- fread(metadata_input)
    } else {
        meta <- as.data.table(metadata_input)
    }
    
    # 2. Split dataset
    # FIX: robust handling of NA values in the Batch column
    # Reference: everything that is NOT "New_Projected", NA INCLUDED
    if ("Batch" %in% names(meta)) {
        ref_samples <- meta[Batch != "New_Projected" | is.na(Batch)]
    } else {
        if(is.null(sample_ids)) stop("Missing 'Batch' column: cannot identify the reference samples automatically.")
        ref_samples <- meta[!Sample_ID_Full %in% sample_ids]
    }
    
    if (nrow(ref_samples) == 0) {
        stop("CRITICAL ERROR: no reference sample found. Check the 'Batch' column of the metadata.")
    }
    
    # Target Samples
    if (is.null(sample_ids)) {
         if ("Batch" %in% names(meta)) {
            new_samples <- meta[Batch == "New_Projected"]
         } else { stop("Give sample_ids or a Batch column.") }
    } else {
        indices <- unique(unlist(lapply(sample_ids, function(x) grep(x, meta$Sample_ID_Full))))
        new_samples <- meta[indices]
    }
    
    if (nrow(new_samples) == 0) {
        warning("No new sample found for the analysis.")
        return(NULL)
    }
    
    # Class column (renamed to 'Class_Label')
    class_col <- NULL
    for(col in c("methylation_class", "diagnosis", "sample.title", "general_class", "subclass")) {
        if(col %in% names(ref_samples)) {
            class_col <- col
            break
        }
    }
    if(is.null(class_col)) {
        ref_samples$Class_Label <- "Unknown"
    } else {
        ref_samples$Class_Label <- ref_samples[[class_col]]
    }
    
    # 3. Distances
    results_list <- list()
    
    message("Computing the 3 nearest neighbours of ", nrow(new_samples), " samples...")
    message("Reference dataset: ", nrow(ref_samples), " samples.")
    
    for (i in 1:nrow(new_samples)) {
        target <- new_samples[i,]
        sid <- target$Sample_ID_Full
        
        # Euclidean distance
        dists <- sqrt((ref_samples$UMAP1 - target$UMAP1)^2 + 
                      (ref_samples$UMAP2 - target$UMAP2)^2)
        
        temp_df <- data.table(
            Ref_ID = ref_samples$Sample_ID_Full,
            Ref_Class = ref_samples$Class_Label,
            Dist = dists
        )
        
        top3 <- temp_df[order(Dist)][1:3]
        
        res_row <- data.table(
            Sample_ID = sid,
            NN1_ID = top3$Ref_ID[1], NN1_Class = top3$Ref_Class[1], NN1_Dist = round(top3$Dist[1], 4),
            NN2_ID = top3$Ref_ID[2], NN2_Class = top3$Ref_Class[2], NN2_Dist = round(top3$Dist[2], 4),
            NN3_ID = top3$Ref_ID[3], NN3_Class = top3$Ref_Class[3], NN3_Dist = round(top3$Dist[3], 4)
        )
        results_list[[i]] <- res_row
    }
    
    final_df <- rbindlist(results_list)
    return(final_df)
}