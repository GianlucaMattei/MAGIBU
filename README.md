# MaGiBu
A modular R framework for stable methylation analysis. MAGIBU projects heterogeneous data (Array/WGBS/Nanopore) onto a pre-computed reference. Features dynamic mapping and quantitative classification based on nearest class representatives to resolve heterogeneity.



1) ottieni i tsv con i beta per ottenere il reference
   funzione: process_methylation_data
   (presumibilmente sono idat)
   
2) genera la umap di riferimenti (ti serve anche l'annotazione)
   funzione: generate_umap_training_model
   generate_umap_training_model(path.450k, path.annotation, "/data2/gmattei/Magibu/Models", "capper.umap", n.probes = 50000, n.pc = 50, n.neighbors = 15, thr = 60)
2a) opzionale, controlla la umap reference (controlla se i campioni che devono essere vicino effettivamente lo sono)
   funzione: analyze_reference_class_neighbors (media globale) analyze_reference_class_neighbors_detailed (dettaglio per ogni campione) plot_reference_umap (plotta la mappa ottenuta fino ad ora)
   
4) ottieni i tsv dei beta dei campioni che vuoi utilizzare
   funzione: process_methylation_data

   a) se sono idat indichi la posizione dei file idat
   b) se sono nanopore indicalo in array.type = "Nanopore" e poi indica le colonne come si chiamano in nanopore.settings = c(chr = "V1", pos = "V2", beta = "V5"
      vanno indicati i tsv generati mediante modkit e resorted di poremeth2
        REF=hg38.fasta # your fasta reference
        THN=10
        modkit extract full --mapped-only --threads ${THN} --cpg --reference ${REF} dorado.output.bam modkit.output.tsv
        sh ModkitResorter.sh modkit.output.tsv

 5) proietta i nuovi campioni sull'UMAP di riferimento
    funzione: project_new_samples
    project_new_samples(training.metadata = "/path/to/umap_metadata_result.tsv", model.components = "/path/to/umap_model.components.rds", uwot.model = "/path/to/umap_umap_model.uwot", path.new.tsv = "/path/to/new_samples", save.path = "/path/to/save_directory")

 6) visualizza i risultati. ci sono vari modi:
    a) plot_umap_zoom_combined restituisce un pannello per ogni campione. Il pannello e costituito da un plot che mostra nel dettaglio (zoom) la umap e dove si posiziona il campione rispetto al reference, sulla destra dei barplot con la media delle distanze dai vari gruppi
    b) get_top_3_distances restituisce un dataframe dove per ogni campione ci viene detto quali sono le tre classi piu vicine e quanto distano in media
    c) plot_global_projection restituisce un pdf con la mappa e i campioni nel suo insieme
    d) plot_umap_highlight_specific restituisce un pdf con la mappa e i campioni nel suo insieme zoomata per far vedere rispetto a una determinata classe dove si posizionano i campioni


 7) Se 

    
