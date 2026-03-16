get_nearest_neighbors <- function(metadata_input, sample_ids = NULL) {
    
    # 1. Caricamento Dati
    if (is.character(metadata_input)) {
        if (!file.exists(metadata_input)) stop("File metadata non trovato.")
        meta <- fread(metadata_input)
    } else {
        meta <- as.data.table(metadata_input)
    }
    
    # 2. Separazione Dataset
    # FIX: Gestione robusta dei valori NA nella colonna Batch
    # Reference: Tutto ciò che NON è "New_Projected", INCLUSI i NA
    if ("Batch" %in% names(meta)) {
        ref_samples <- meta[Batch != "New_Projected" | is.na(Batch)]
    } else {
        if(is.null(sample_ids)) stop("Colonna 'Batch' mancante: impossibile distinguere reference automaticamente.")
        ref_samples <- meta[!Sample_ID_Full %in% sample_ids]
    }
    
    if (nrow(ref_samples) == 0) {
        stop("ERRORE CRITICO: Nessun campione di reference trovato. Controlla la colonna 'Batch' dei metadati.")
    }
    
    # Target Samples
    if (is.null(sample_ids)) {
         if ("Batch" %in% names(meta)) {
            new_samples <- meta[Batch == "New_Projected"]
         } else { stop("Definire sample_ids o avere colonna Batch.") }
    } else {
        indices <- unique(unlist(lapply(sample_ids, function(x) grep(x, meta$Sample_ID_Full))))
        new_samples <- meta[indices]
    }
    
    if (nrow(new_samples) == 0) {
        warning("Nessun nuovo campione trovato per l'analisi.")
        return(NULL)
    }
    
    # Gestione Colonna Classe (Uniformiamo il nome a 'Class_Label')
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
    
    # 3. Calcolo Distanze
    results_list <- list()
    
    message("Calcolo dei 3 vicini più prossimi per ", nrow(new_samples), " campioni...")
    message("Dataset Reference: ", nrow(ref_samples), " campioni.")
    
    for (i in 1:nrow(new_samples)) {
        target <- new_samples[i,]
        sid <- target$Sample_ID_Full
        
        # Calcola distanza euclidea
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