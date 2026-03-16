library(data.table)
library(umap)
library(ggplot2)
library(parallel)
library(matrixStats)
library(irlba)
library(uwot)
library(dplyr)

# ==============================================================================
# FUNZIONE: generate_umap_training_model
# DESCRIZIONE: 
#   Questa funzione elabora i dati di metilazione (TSV), costruisce un modello
#   PCA e UMAP basato sulla strategia originale, e salva i componenti necessari
#   per proiettare futuri campioni.
#
# INPUT:
#   - training: Percorso alla cartella contenente i file .tsv di training
#   - path.to.annotation.data: Percorso al file annotazioni (.tsv)
#   - save.path: Percorso dove salvare i risultati (model.components, umap, metadata)
# ==============================================================================
plot_umap_zoom <- function(
    metadata, 
    sample_id, 
    fixed.colors = NULL, 
    circle_radius = 1.5, 
    distance_multiplier = 2, 
    n_top_groups = 10, 
    color_palette = NULL
) {

    # --- 1. Controlli Preliminari ---
    # Verifica esistenza colonne
    req_cols <- c("UMAP1", "UMAP2", "Sample_ID_Full")
    if (!all(req_cols %in% names(metadata))) {
        stop("Il metadata deve contenere le colonne: UMAP1, UMAP2, Sample_ID_Full")
    }

    # --- 2. Identificazione Target ---
    # Cerca il campione usando pattern matching
    target.idx_all <- grep(sample_id, metadata$Sample_ID_Full)
    
    if (length(target.idx_all) == 0) {
        warning(paste("ERRORE: Campione '", sample_id, "' non trovato nel metadata!"))
        return(NULL)
    }
    
    # Se ce n'è più di uno, prendiamo il primo ma avvisiamo
    target.idx_orig <- target.idx_all[1]
    target_row <- metadata[target.idx_orig, ]
    full_id <- target_row$Sample_ID_Full
    
    if (length(target.idx_all) > 1) {
        message("Nota: Trovati più campioni corrispondenti a '", sample_id, "'. Uso il primo: ", full_id)
    }

    # --- 3. Pulizia Metadata (Nascondi altri proiettati) ---
    # Logica: Se ci sono altri campioni proiettati (New_Projected) che non sono
    # il nostro target, li nascondiamo per pulire il plot.
    meta_proc <- metadata
    if ("Batch" %in% names(metadata)) {
        # Trova indici di TUTTI i nuovi proiettati
        all_new_indices <- which(metadata$Batch == "New_Projected")
        # Rimuovi dalla lista da nascondere il NOSTRO target
        indices_to_hide <- setdiff(all_new_indices, target.idx_orig)
        
        if (length(indices_to_hide) > 0) {
            meta_proc <- metadata[-indices_to_hide, ]
        }
    }
    
    # Ricalcoliamo l'indice del target nel dataset pulito
    target.idx <- which(meta_proc$Sample_ID_Full == full_id)
    target_row <- meta_proc[target.idx, ]

    # --- 4. Calcolo Distanze e Zoom ---
    # Calcola distanza euclidea di tutti i punti dal target
    meta_proc$dist_to_target <- sqrt((meta_proc$UMAP1 - target_row$UMAP1)^2 + 
                                     (meta_proc$UMAP2 - target_row$UMAP2)^2)
    
    # Definiamo i vicini (dentro il cerchio)
    neighbors_idx <- which(meta_proc$dist_to_target <= circle_radius)
    
    # Definiamo l'area di zoom (un rettangolo intorno al cerchio)
    zoom_radius <- circle_radius * distance_multiplier
    xlims <- c(target_row$UMAP1 - zoom_radius, target_row$UMAP1 + zoom_radius)
    ylims <- c(target_row$UMAP2 - zoom_radius, target_row$UMAP2 + zoom_radius)

    # --- 5. Gestione Colori ---
    # Identifichiamo la colonna che contiene la classe (es. methylation_class, diagnosis, etc)
    # Cerchiamo in ordine di preferenza
    class_cols <- c("methylation_class", "diagnosis", "subclass", "general_class", "sample.title")
    class_col <- intersect(class_cols, names(meta_proc))[1]
    
    if (is.na(class_col)) {
        warning("Nessuna colonna di classe trovata (es. methylation_class). Uso colori default.")
        meta_proc$PlotClass <- "Unknown"
    } else {
        meta_proc$PlotClass <- meta_proc[[class_col]]
    }

    # Filtriamo i gruppi più frequenti nell'area visualizzata per la legenda
    data_in_zoom <- meta_proc[meta_proc$UMAP1 >= xlims[1] & meta_proc$UMAP1 <= xlims[2] &
                              meta_proc$UMAP2 >= ylims[1] & meta_proc$UMAP2 <= ylims[2], ]
    
    top_groups <- data_in_zoom %>%
        count(PlotClass) %>%
        arrange(desc(n)) %>%
        slice_head(n = n_top_groups) %>%
        pull(PlotClass)
    
    # Tutto ciò che non è top group diventa "Other"
    meta_proc$LegendLabel <- ifelse(meta_proc$PlotClass %in% top_groups, 
                                    meta_proc$PlotClass, "Other")
    
    # Assegnazione Colori
    if (is.null(fixed.colors)) {
        unique_groups <- unique(c(top_groups, "Other"))
        n_colors <- length(unique_groups)
        
        if (!is.null(color_palette)) {
            my_cols <- colorRampPalette(color_palette)(n_colors)
        } else {
            # Palette default qualitativa
             my_cols <- scales::hue_pal()(n_colors)
        }
        names(my_cols) <- unique_groups
        my_cols["Other"] <- "grey80" # Colore fisso per Other
    } else {
        my_cols <- fixed.colors
        if (!"Other" %in% names(my_cols)) my_cols["Other"] <- "grey80"
    }

    # --- 6. Creazione Plot ---
    p <- ggplot(meta_proc, aes(x = UMAP1, y = UMAP2)) +
        # A. Punti di sfondo (grigi o colorati sbiaditi)
        geom_point(aes(color = LegendLabel), alpha = 0.4, size = 1) +
        
        # B. Cerchio di evidenziazione
        annotate("path",
                 x = target_row$UMAP1 + circle_radius * cos(seq(0, 2*pi, length.out=100)),
                 y = target_row$UMAP2 + circle_radius * sin(seq(0, 2*pi, length.out=100)),
                 color = "black", linetype = "dashed", size = 0.8) +
        
        # C. Il punto Target (Stella o punto grosso)
        geom_point(data = target_row, color = "black", fill = "red", shape = 23, size = 4, stroke = 1.5) +
        
        # D. Etichetta del Target
        geom_label_repel(data = target_row, aes(label = Sample_ID_Full), 
                         box.padding = 1, point.padding = 0.5, 
                         fill = "white", alpha = 0.9, fontface = "bold") +
        
        # E. Etichette per i vicini (opzionale: solo se pochi)
        # geom_text_repel(data = meta_proc[neighbors_idx,], aes(label = PlotClass), max.overlaps = 10, size = 3) +

        scale_color_manual(values = my_cols, name = "Class") +
        coord_cartesian(xlim = xlims, ylim = ylims) +
        theme_minimal() +
        labs(title = paste("Zoom on:", full_id),
             subtitle = paste("Radius:", circle_radius),
             x = "UMAP 1", y = "UMAP 2") +
        theme(legend.position = "right")

    return(p)
}