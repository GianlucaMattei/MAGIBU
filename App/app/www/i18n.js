// Interface translation. The page is written in English; every visible text that matches an entry
// below is replaced in the chosen language (also text that appears later, e.g. status messages).
// Plots, tables content and output files stay in English. Columns: fr, es, it, zh; German (de) in the DE table below.
(function () {
  var D = {
    // navigation
    "Prediction": ["Prédiction", "Predicción", "Predizione", "预测"],
    "Global map": ["Carte globale", "Mapa global", "Mappa globale", "全局图"],
    "Class zoom": ["Zoom sur les classes", "Zoom por clases", "Zoom classi", "类别放大"],
    "Build model": ["Créer un modèle", "Crear modelo", "Crea modello", "构建模型"],
    "About": ["À propos", "Acerca de", "Informazioni", "关于"],
    // prediction
    "1. Model": ["1. Modèle", "1. Modelo", "1. Modello", "1. 模型"],
    "Model folder": ["Dossier du modèle", "Carpeta del modelo", "Cartella del modello", "模型文件夹"],
    "2. Samples": ["2. Échantillons", "2. Muestras", "2. Campioni", "2. 样本"],
    "Single file": ["Fichier unique", "Archivo único", "File singolo", "单个文件"],
    "All files in a folder": ["Tous les fichiers d'un dossier", "Todos los archivos de una carpeta", "Tutti i file di una cartella", "文件夹中的所有文件"],
    "Selected files": ["Fichiers sélectionnés", "Archivos seleccionados", "File selezionati", "所选文件"],
    "Paths (one per line)": ["Chemins (un par ligne)", "Rutas (una por línea)", "Percorsi (uno per riga)", "路径（每行一个）"],
    "idat (_Grn or _Red, also .gz) or tsv": ["idat (_Grn ou _Red, aussi .gz) ou tsv", "idat (_Grn o _Red, también .gz) o tsv", "idat (_Grn o _Red, anche .gz) oppure tsv", "idat（_Grn 或 _Red，也可 .gz）或 tsv"],
    "Array type": ["Type de puce", "Tipo de array", "Tipo di array", "芯片类型"],
    "Automatic": ["Automatique", "Automático", "Automatico", "自动"],
    "3. Results": ["3. Résultats", "3. Resultados", "3. Risultati", "3. 结果"],
    "Output folder": ["Dossier de sortie", "Carpeta de salida", "Cartella di uscita", "输出文件夹"],
    "Analysis name": ["Nom de l'analyse", "Nombre del análisis", "Nome dell'analisi", "分析名称"],
    "Cores to use": ["Cœurs à utiliser", "Núcleos a usar", "Core da usare", "使用的核心数"],
    "Run prediction": ["Lancer la prédiction", "Iniciar predicción", "Avvia predizione", "运行预测"],
    "Stop": ["Arrêter", "Detener", "Interrompi", "停止"],
    "Reopen a previous analysis (folder)": ["Rouvrir une analyse précédente (dossier)", "Reabrir un análisis anterior (carpeta)", "Riapri un'analisi precedente (cartella)", "重新打开以前的分析（文件夹）"],
    "Load": ["Charger", "Cargar", "Carica", "加载"],
    "Browse...": ["Parcourir...", "Examinar...", "Sfoglia...", "浏览..."],
    "Summary": ["Résumé", "Resumen", "Riepilogo", "摘要"],
    "Sample names": ["Noms des échantillons", "Nombres de las muestras", "Nomi dei campioni", "样本名称"],
    "Distances by class": ["Distances par classe", "Distancias por clase", "Distanze per classe", "按类别的距离"],
    "Classic plot": ["Graphique classique", "Gráfico clásico", "Plot classico", "经典图"],
    "Log": ["Journal", "Registro", "Registro", "日志"],
    "Download distances (Excel)": ["Télécharger les distances (Excel)", "Descargar distancias (Excel)", "Scarica distanze (Excel)", "下载距离（Excel）"],
    "Open results folder": ["Ouvrir le dossier des résultats", "Abrir carpeta de resultados", "Apri cartella risultati", "打开结果文件夹"],
    "Give each sample a name (for example Case 1). Names appear in tables, plots and in the Excel file; the file name stays as sample ID.": [
      "Donnez un nom à chaque échantillon (par exemple Case 1). Les noms apparaissent dans les tableaux, les graphiques et le fichier Excel ; le nom du fichier reste l'identifiant de l'échantillon.",
      "Dé un nombre a cada muestra (por ejemplo Case 1). Los nombres aparecen en tablas, gráficos y en el archivo Excel; el nombre del archivo sigue siendo el ID de la muestra.",
      "Dai un nome a ogni campione (per esempio Case 1). I nomi compaiono in tabelle, grafici e nel file Excel; il nome del file resta l'ID del campione.",
      "为每个样本命名（例如 Case 1）。名称会显示在表格、图和 Excel 文件中；文件名仍作为样本 ID。"],
    "Number as Case 1, 2, 3...": ["Numéroter Case 1, 2, 3...", "Numerar como Case 1, 2, 3...", "Numera come Case 1, 2, 3...", "编号为 Case 1, 2, 3..."],
    "Clear names": ["Effacer les noms", "Borrar nombres", "Cancella nomi", "清除名称"],
    "Save names": ["Enregistrer les noms", "Guardar nombres", "Salva nomi", "保存名称"],
    "Sample": ["Échantillon", "Muestra", "Campione", "样本"],
    "Download PDF (all samples)": ["Télécharger le PDF (tous les échantillons)", "Descargar PDF (todas las muestras)", "Scarica PDF (tutti i campioni)", "下载 PDF（所有样本）"],
    "Sample ID": ["ID de l'échantillon", "ID de la muestra", "ID campione", "样本 ID"],
    "Name": ["Nom", "Nombre", "Nome", "名称"],
    "Closest class": ["Classe la plus proche", "Clase más cercana", "Classe più vicina", "最近的类别"],
    "Run or load an analysis first.": ["Lancez ou chargez d'abord une analyse.", "Primero ejecute o cargue un análisis.", "Prima esegui o carica un'analisi.", "请先运行或加载一个分析。"],
    // global map
    "Colour by": ["Couleur selon", "Color según", "Colora per", "着色依据"],
    "Class": ["Classe", "Clase", "Classe", "类别"],
    "Family (text before the comma)": ["Famille (texte avant la virgule)", "Familia (texto antes de la coma)", "Famiglia (testo prima della virgola)", "家族（逗号前的文字）"],
    "Show sample names": ["Afficher les noms des échantillons", "Mostrar nombres de las muestras", "Mostra i nomi dei campioni", "显示样本名称"],
    "Zoom on samples": ["Zoom sur les échantillons", "Zoom en las muestras", "Zoom sui campioni", "放大到样本"],
    "Whole map": ["Carte entière", "Mapa completo", "Mappa intera", "整张图"],
    "Mouse wheel: zoom in and out.": ["Molette de la souris : zoom avant et arrière.", "Rueda del ratón: acercar y alejar.", "Rotellina del mouse: ingrandisci e riduci.", "鼠标滚轮：放大和缩小。"],
    "Drag: zoom on an area.": ["Glisser : zoom sur une zone.", "Arrastrar: zoom en un área.", "Trascina: zoom su un'area.", "拖动：放大某个区域。"],
    "Double click: back to the whole map.": ["Double clic : retour à la carte entière.", "Doble clic: volver al mapa completo.", "Doppio clic: torna alla mappa intera.", "双击：返回整张图。"],
    "Click a class in the legend to hide it, double click to show only that one.": [
      "Cliquez sur une classe dans la légende pour la masquer, double-cliquez pour n'afficher qu'elle.",
      "Haga clic en una clase de la leyenda para ocultarla, doble clic para mostrar solo esa.",
      "Clic su una classe nella legenda per nasconderla, doppio clic per mostrare solo quella.",
      "单击图例中的类别可隐藏它，双击则只显示该类别。"],
    "Pass over a point to see its class.": ["Survolez un point pour voir sa classe.", "Pase sobre un punto para ver su clase.", "Passa sopra un punto per vederne la classe.", "将鼠标悬停在点上可查看其类别。"],
    "Download the current view:": ["Télécharger la vue actuelle :", "Descargar la vista actual:", "Scarica la vista attuale:", "下载当前视图："],
    // class zoom
    "The view widens until it contains every reference sample of the chosen classes and the chosen samples.": [
      "La vue s'élargit jusqu'à contenir tous les échantillons de référence des classes choisies et les échantillons choisis.",
      "La vista se amplía hasta contener todas las muestras de referencia de las clases elegidas y las muestras elegidas.",
      "La vista si allarga fino a contenere tutti i campioni di riferimento delle classi scelte e i campioni scelti.",
      "视图会扩大，直到包含所选类别的全部参考样本以及所选样本。"],
    "Samples to show": ["Échantillons à afficher", "Muestras a mostrar", "Campioni da mostrare", "要显示的样本"],
    "Classes to highlight": ["Classes à mettre en évidence", "Clases a resaltar", "Classi da evidenziare", "要突出显示的类别"],
    "choose one or more classes": ["choisissez une ou plusieurs classes", "elija una o más clases", "scegli una o più classi", "选择一个或多个类别"],
    "Add the closest": ["Ajouter les plus proches", "Añadir las más cercanas", "Aggiungi le più vicine", "添加最近的类别"],
    "Whole family": ["Toute la famille", "Toda la familia", "Tutta la famiglia", "整个家族"],
    "Clear": ["Vider", "Vaciar", "Svuota", "清空"],
    "Show other reference samples in grey": ["Afficher en gris les autres échantillons de référence", "Mostrar en gris las demás muestras de referencia", "Mostra in grigio gli altri campioni di riferimento", "以灰色显示其他参考样本"],
    "Class labels on the plot": ["Étiquettes des classes sur le graphique", "Etiquetas de clases en el gráfico", "Etichette delle classi nel grafico", "图中显示类别标签"],
    "Sample labels on the plot": ["Étiquettes des échantillons sur le graphique", "Etiquetas de muestras en el gráfico", "Etichette dei campioni nel grafico", "图中显示样本标签"],
    "Sample names are set in Prediction > Sample names.": ["Les noms des échantillons se définissent dans Prédiction > Noms des échantillons.", "Los nombres de las muestras se definen en Predicción > Nombres de las muestras.", "I nomi dei campioni si impostano in Predizione > Nomi dei campioni.", "样本名称在“预测 > 样本名称”中设置。"],
    "Title": ["Titre", "Título", "Titolo", "标题"],
    // build model
    "Dataset folder (idat or tsv)": ["Dossier du jeu de données (idat ou tsv)", "Carpeta del conjunto de datos (idat o tsv)", "Cartella del dataset (idat o tsv)", "数据集文件夹（idat 或 tsv）"],
    "Annotation file (tsv, csv)": ["Fichier d'annotation (tsv, csv)", "Archivo de anotación (tsv, csv)", "File di annotazione (tsv, csv)", "注释文件（tsv、csv）"],
    "Sample column": ["Colonne échantillon", "Columna de muestra", "Colonna campione", "样本列"],
    "Class column": ["Colonne classe", "Columna de clase", "Colonna classe", "类别列"],
    "Folder where the model is saved": ["Dossier où enregistrer le modèle", "Carpeta donde guardar el modelo", "Cartella dove salvare il modello", "模型保存文件夹"],
    "Model name": ["Nom du modèle", "Nombre del modelo", "Nome del modello", "模型名称"],
    "Advanced parameters": ["Paramètres avancés", "Parámetros avanzados", "Parametri avanzati", "高级参数"],
    "Most variable CpGs": ["CpG les plus variables", "CpG más variables", "CpG più variabili", "变异最大的 CpG 数"],
    "PCA components": ["Composantes ACP", "Componentes PCA", "Componenti PCA", "PCA 主成分数"],
    "UMAP neighbours": ["Voisins UMAP", "Vecinos UMAP", "Vicini UMAP", "UMAP 邻居数"],
    "UMAP minimum distance": ["Distance minimale UMAP", "Distancia mínima UMAP", "Distanza minima UMAP", "UMAP 最小距离"],
    "Model map": ["Carte du modèle", "Mapa del modelo", "Mappa del modello", "模型图"],
    // about
    "For research use only.": ["Réservé à la recherche.", "Solo para uso en investigación.", "Solo per uso di ricerca.", "仅供研究使用。"],
    "Not a diagnostic device: the result does not replace the diagnosis.": [
      "Ce n'est pas un dispositif de diagnostic : le résultat ne remplace pas le diagnostic.",
      "No es un dispositivo de diagnóstico: el resultado no sustituye al diagnóstico.",
      "Non è un dispositivo diagnostico: il risultato non sostituisce la diagnosi.",
      "本软件不是诊断器械：结果不能替代诊断。"],
    "The prediction projects each sample onto the UMAP map of the chosen model and measures its distances from the reference samples.": [
      "La prédiction projette chaque échantillon sur la carte UMAP du modèle choisi et mesure ses distances aux échantillons de référence.",
      "La predicción proyecta cada muestra sobre el mapa UMAP del modelo elegido y mide sus distancias a las muestras de referencia.",
      "La predizione proietta ogni campione sulla mappa UMAP del modello scelto e misura le distanze dai campioni di riferimento.",
      "预测会将每个样本投影到所选模型的 UMAP 图上，并测量其与参考样本的距离。"],
    "The UMAP projection is stochastic: running the same sample again can move it slightly (about 0.2 UMAP units), as in the original Magibu pipeline.": [
      "La projection UMAP est stochastique : relancer le même échantillon peut le déplacer légèrement (environ 0,2 unité UMAP), comme dans la pipeline Magibu d'origine.",
      "La proyección UMAP es estocástica: volver a procesar la misma muestra puede desplazarla ligeramente (unas 0,2 unidades UMAP), como en el pipeline original de Magibu.",
      "La proiezione UMAP è stocastica: rilanciando lo stesso campione può spostarsi leggermente (circa 0,2 unità UMAP), come nella pipeline Magibu originale.",
      "UMAP 投影具有随机性：再次运行同一样本，其位置可能略有变化（约 0.2 个 UMAP 单位），与原始 Magibu 流程相同。"],
    "Third party software included (R and R packages, Electron): see the THIRD_PARTY_LICENSES file in the program folder.": [
      "Logiciels tiers inclus (R et paquets R, Electron) : voir le fichier THIRD_PARTY_LICENSES dans le dossier du programme.",
      "Software de terceros incluido (R y paquetes de R, Electron): ver el archivo THIRD_PARTY_LICENSES en la carpeta del programa.",
      "Software di terzi inclusi (R e pacchetti R, Electron): vedere il file THIRD_PARTY_LICENSES nella cartella del programma.",
      "包含的第三方软件（R 及 R 包、Electron）：请参阅程序文件夹中的 THIRD_PARTY_LICENSES 文件。"],
    // status messages
    "A job is already running.": ["Une tâche est déjà en cours.", "Ya hay un trabajo en curso.", "C'è già un lavoro in corso.", "已有任务正在运行。"],
    "Job running... (details in the Log tab)": ["Tâche en cours... (détails dans l'onglet Journal)", "Trabajo en curso... (detalles en la pestaña Registro)", "Lavoro in corso... (dettagli nella scheda Registro)", "任务运行中……（详情见“日志”选项卡）"],
    "Job stopped.": ["Tâche interrompue.", "Trabajo detenido.", "Lavoro interrotto.", "任务已停止。"],
    "Choose the model folder.": ["Choisissez le dossier du modèle.", "Elija la carpeta del modelo.", "Indica la cartella del modello.", "请选择模型文件夹。"],
    "Give at least one sample file or folder.": ["Indiquez au moins un fichier ou dossier d'échantillons.", "Indique al menos un archivo o carpeta de muestras.", "Indica almeno un file o una cartella di campioni.", "请至少指定一个样本文件或文件夹。"],
    "Choose the output folder.": ["Choisissez le dossier de sortie.", "Elija la carpeta de salida.", "Indica la cartella di uscita.", "请选择输出文件夹。"],
    "Saving names...": ["Enregistrement des noms...", "Guardando nombres...", "Salvataggio dei nomi...", "正在保存名称……"],
    "Names saved: Excel file and classic PDF updated.": ["Noms enregistrés : fichier Excel et PDF classique mis à jour.", "Nombres guardados: archivo Excel y PDF clásico actualizados.", "Nomi salvati: file Excel e PDF classico aggiornati.", "名称已保存：Excel 文件和经典 PDF 已更新。"],
    "Choose a valid model in the Prediction tab, or load an analysis.": ["Choisissez un modèle valide dans l'onglet Prédiction, ou chargez une analyse.", "Elija un modelo válido en la pestaña Predicción, o cargue un análisis.", "Scegli un modello valido nella scheda Predizione, o carica un'analisi.", "请在“预测”选项卡中选择有效的模型，或加载一个分析。"],
    "No analysis loaded: only the model map is shown. Samples appear after a prediction.": [
      "Aucune analyse chargée : seule la carte du modèle est affichée. Les échantillons apparaissent après une prédiction.",
      "Ningún análisis cargado: solo se muestra el mapa del modelo. Las muestras aparecen tras una predicción.",
      "Nessuna analisi caricata: si vede solo la mappa del modello. I campioni compaiono dopo una predizione.",
      "未加载分析：仅显示模型图。运行预测后会显示样本。"],
    "No samples: run or load an analysis first.": ["Aucun échantillon : lancez ou chargez d'abord une analyse.", "Sin muestras: primero ejecute o cargue un análisis.", "Nessun campione: prima esegui o carica un'analisi.", "没有样本：请先运行或加载一个分析。"],
    "Choose the dataset folder.": ["Choisissez le dossier du jeu de données.", "Elija la carpeta del conjunto de datos.", "Indica la cartella del dataset.", "请选择数据集文件夹。"],
    "Choose the annotation file.": ["Choisissez le fichier d'annotation.", "Elija el archivo de anotación.", "Indica il file di annotazione.", "请选择注释文件。"],
    "Choose the folder where the model is saved.": ["Choisissez le dossier où enregistrer le modèle.", "Elija la carpeta donde guardar el modelo.", "Indica la cartella dove salvare il modello.", "请选择模型保存文件夹。"],
    // engine errors (the part after ": " is a path or a list and stays as it is)
    "Model folder not found": ["Dossier du modèle introuvable", "Carpeta del modelo no encontrada", "Cartella del modello non trovata", "未找到模型文件夹"],
    "Model file missing": ["Fichier du modèle manquant", "Falta un archivo del modelo", "Manca un file del modello", "缺少模型文件"],
    "More than one model in the folder": ["Plusieurs modèles dans le dossier", "Más de un modelo en la carpeta", "Più di un modello nella cartella", "文件夹中有多个模型"],
    "Path not found": ["Chemin introuvable", "Ruta no encontrada", "Percorso non trovato", "未找到路径"],
    "Unrecognised file (expected *_Grn.idat / *_Red.idat, also .gz, or .tsv)": [
      "Fichier non reconnu (attendu *_Grn.idat / *_Red.idat, aussi .gz, ou .tsv)",
      "Archivo no reconocido (se esperaba *_Grn.idat / *_Red.idat, también .gz, o .tsv)",
      "File non riconosciuto (attesi *_Grn.idat / *_Red.idat, anche .gz, oppure .tsv)",
      "无法识别的文件（应为 *_Grn.idat / *_Red.idat，也可 .gz，或 .tsv）"],
    "Same sample given more than once": ["Même échantillon indiqué plusieurs fois", "Misma muestra indicada más de una vez", "Stesso campione indicato più volte", "同一样本被指定了多次"],
    "Missing file (both channels _Grn and _Red are needed)": ["Fichier manquant (les deux canaux _Grn et _Red sont nécessaires)", "Falta un archivo (se necesitan ambos canales _Grn y _Red)", "Manca un file (servono entrambi i canali _Grn e _Red)", "缺少文件（需要 _Grn 和 _Red 两个通道）"],
    "idat conversion failed for": ["Échec de la conversion idat pour", "Falló la conversión idat de", "Conversione idat fallita per", "idat 转换失败"],
    "No Magibu result (*.magibu_result.rds) in": ["Aucun résultat Magibu (*.magibu_result.rds) dans", "Ningún resultado Magibu (*.magibu_result.rds) en", "Nessun risultato Magibu (*.magibu_result.rds) in", "未找到 Magibu 结果（*.magibu_result.rds）"],
    "Annotation file not found": ["Fichier d'annotation introuvable", "Archivo de anotación no encontrado", "File di annotazione non trovato", "未找到注释文件"],
    "Columns not found in the annotation": ["Colonnes introuvables dans l'annotation", "Columnas no encontradas en la anotación", "Colonne non trovate nell'annotazione", "注释中未找到列"],
    "Dataset folder not found": ["Dossier du jeu de données introuvable", "Carpeta del conjunto de datos no encontrada", "Cartella del dataset non trovata", "未找到数据集文件夹"],
    "The destination folder already contains a model, choose a new or empty one": [
      "Le dossier de destination contient déjà un modèle, choisissez-en un nouveau ou vide",
      "La carpeta de destino ya contiene un modelo, elija una nueva o vacía",
      "La cartella di destinazione contiene già un modello, scegline una nuova o vuota",
      "目标文件夹中已有模型，请选择新的或空的文件夹"],
    "No file or folder given.": ["Aucun fichier ou dossier indiqué.", "No se indicó ningún archivo o carpeta.", "Nessun file o cartella indicato.", "未指定文件或文件夹。"],
    "No idat or tsv file found in the given paths.": ["Aucun fichier idat ou tsv trouvé dans les chemins indiqués.", "No se encontró ningún archivo idat o tsv en las rutas indicadas.", "Nessun file idat o tsv trovato nei percorsi indicati.", "在指定路径中未找到 idat 或 tsv 文件。"],
    "The model has no annotated classes: distances cannot be computed.": ["Le modèle n'a pas de classes annotées : impossible de calculer les distances.", "El modelo no tiene clases anotadas: no se pueden calcular las distancias.", "Il modello non ha classi annotate: impossibile calcolare le distanze.", "模型没有已注释的类别：无法计算距离。"],
    "The annotation needs at least two columns: sample identifier and class.": ["L'annotation doit avoir au moins deux colonnes : identifiant de l'échantillon et classe.", "La anotación necesita al menos dos columnas: identificador de muestra y clase.", "L'annotazione deve avere almeno due colonne: identificativo del campione e classe.", "注释至少需要两列：样本标识和类别。"],
    // tables
    "Search:": ["Rechercher :", "Buscar:", "Cerca:", "搜索："],
    "Previous": ["Précédent", "Anterior", "Precedente", "上一页"],
    "Next": ["Suivant", "Siguiente", "Successivo", "下一页"],
    "No data available in table": ["Aucune donnée disponible", "No hay datos disponibles", "Nessun dato disponibile", "表中无数据"],
    "No matching records found": ["Aucun résultat", "No se encontraron resultados", "Nessun risultato", "未找到匹配记录"]
  };

  // German: kept apart and appended as 5th column
  var DE = { "Prediction": "Vorhersage", "Global map": "Gesamtkarte", "Class zoom": "Klassen-Zoom", "Build model": "Modell erstellen", "About": "Info", "1. Model": "1. Modell", "Model folder": "Modellordner", "2. Samples": "2. Proben", "Single file": "Einzelne Datei", "All files in a folder": "Alle Dateien eines Ordners", "Selected files": "Ausgewählte Dateien", "Paths (one per line)": "Pfade (einer pro Zeile)", "idat (_Grn or _Red, also .gz) or tsv": "idat (_Grn oder _Red, auch .gz) oder tsv", "Array type": "Array-Typ", "Automatic": "Automatisch", "3. Results": "3. Ergebnisse", "Output folder": "Ausgabeordner", "Analysis name": "Name der Analyse", "Cores to use": "Zu verwendende Kerne", "Run prediction": "Vorhersage starten", "Stop": "Abbrechen", "Reopen a previous analysis (folder)": "Frühere Analyse öffnen (Ordner)", "Load": "Laden", "Browse...": "Durchsuchen...", "Summary": "Übersicht", "Sample names": "Probennamen", "Distances by class": "Abstände nach Klasse", "Classic plot": "Klassische Grafik", "Log": "Protokoll", "Download distances (Excel)": "Abstände herunterladen (Excel)", "Open results folder": "Ergebnisordner öffnen", "Give each sample a name (for example Case 1). Names appear in tables, plots and in the Excel file; the file name stays as sample ID.": "Geben Sie jeder Probe einen Namen (zum Beispiel Case 1). Die Namen erscheinen in Tabellen, Grafiken und in der Excel-Datei; der Dateiname bleibt die Proben-ID.", "Number as Case 1, 2, 3...": "Als Case 1, 2, 3... nummerieren", "Clear names": "Namen löschen", "Save names": "Namen speichern", "Sample": "Probe", "Download PDF (all samples)": "PDF herunterladen (alle Proben)", "Sample ID": "Proben-ID", "Name": "Name", "Closest class": "Nächste Klasse", "Run or load an analysis first.": "Zuerst eine Analyse starten oder laden.", "Colour by": "Färben nach", "Class": "Klasse", "Family (text before the comma)": "Familie (Text vor dem Komma)", "Show sample names": "Probennamen anzeigen", "Zoom on samples": "Auf Proben zoomen", "Whole map": "Ganze Karte", "Mouse wheel: zoom in and out.": "Mausrad: vergrößern und verkleinern.", "Drag: zoom on an area.": "Ziehen: auf einen Bereich zoomen.", "Double click: back to the whole map.": "Doppelklick: zurück zur ganzen Karte.", "Click a class in the legend to hide it, double click to show only that one.": "Klick auf eine Klasse in der Legende blendet sie aus, Doppelklick zeigt nur diese.", "Pass over a point to see its class.": "Mit der Maus über einen Punkt fahren, um seine Klasse zu sehen.", "Download the current view:": "Aktuelle Ansicht herunterladen:", "The view widens until it contains every reference sample of the chosen classes and the chosen samples.": "Die Ansicht wird erweitert, bis sie alle Referenzproben der gewählten Klassen und die gewählten Proben enthält.", "Samples to show": "Anzuzeigende Proben", "Classes to highlight": "Hervorzuhebende Klassen", "choose one or more classes": "eine oder mehrere Klassen wählen", "Add the closest": "Nächste hinzufügen", "Whole family": "Ganze Familie", "Clear": "Leeren", "Show other reference samples in grey": "Andere Referenzproben grau anzeigen", "Class labels on the plot": "Klassenbeschriftungen in der Grafik", "Sample labels on the plot": "Probenbeschriftungen in der Grafik", "Sample names are set in Prediction > Sample names.": "Probennamen werden unter Vorhersage > Probennamen festgelegt.", "Title": "Titel", "Dataset folder (idat or tsv)": "Datensatzordner (idat oder tsv)", "Annotation file (tsv, csv)": "Annotationsdatei (tsv, csv)", "Sample column": "Probenspalte", "Class column": "Klassenspalte", "Folder where the model is saved": "Ordner zum Speichern des Modells", "Model name": "Modellname", "Advanced parameters": "Erweiterte Parameter", "Most variable CpGs": "Variabelste CpGs", "PCA components": "PCA-Komponenten", "UMAP neighbours": "UMAP-Nachbarn", "UMAP minimum distance": "UMAP-Mindestabstand", "Model map": "Modellkarte", "For research use only.": "Nur für Forschungszwecke.", "Not a diagnostic device: the result does not replace the diagnosis.": "Kein Diagnoseprodukt: Das Ergebnis ersetzt nicht die Diagnose.", "The prediction projects each sample onto the UMAP map of the chosen model and measures its distances from the reference samples.": "Die Vorhersage projiziert jede Probe auf die UMAP-Karte des gewählten Modells und misst ihre Abstände zu den Referenzproben.", "The UMAP projection is stochastic: running the same sample again can move it slightly (about 0.2 UMAP units), as in the original Magibu pipeline.": "Die UMAP-Projektion ist stochastisch: Wird dieselbe Probe erneut verarbeitet, kann sie sich leicht verschieben (etwa 0,2 UMAP-Einheiten), wie in der ursprünglichen Magibu-Pipeline.", "Third party software included (R and R packages, Electron): see the THIRD_PARTY_LICENSES file in the program folder.": "Enthaltene Software Dritter (R und R-Pakete, Electron): siehe Datei THIRD_PARTY_LICENSES im Programmordner.", "A job is already running.": "Es läuft bereits ein Auftrag.", "Job running... (details in the Log tab)": "Auftrag läuft... (Details im Reiter Protokoll)", "Job stopped.": "Auftrag abgebrochen.", "Choose the model folder.": "Modellordner wählen.", "Give at least one sample file or folder.": "Mindestens eine Probendatei oder einen Ordner angeben.", "Choose the output folder.": "Ausgabeordner wählen.", "Saving names...": "Namen werden gespeichert...", "Names saved: Excel file and classic PDF updated.": "Namen gespeichert: Excel-Datei und klassisches PDF aktualisiert.", "Choose a valid model in the Prediction tab, or load an analysis.": "Im Reiter Vorhersage ein gültiges Modell wählen oder eine Analyse laden.", "No analysis loaded: only the model map is shown. Samples appear after a prediction.": "Keine Analyse geladen: Nur die Modellkarte wird angezeigt. Die Proben erscheinen nach einer Vorhersage.", "No samples: run or load an analysis first.": "Keine Proben: zuerst eine Analyse starten oder laden.", "Choose the dataset folder.": "Datensatzordner wählen.", "Choose the annotation file.": "Annotationsdatei wählen.", "Choose the folder where the model is saved.": "Ordner zum Speichern des Modells wählen.", "Model folder not found": "Modellordner nicht gefunden", "Model file missing": "Modelldatei fehlt", "More than one model in the folder": "Mehr als ein Modell im Ordner", "Path not found": "Pfad nicht gefunden", "Unrecognised file (expected *_Grn.idat / *_Red.idat, also .gz, or .tsv)": "Unbekannte Datei (erwartet *_Grn.idat / *_Red.idat, auch .gz, oder .tsv)", "Same sample given more than once": "Dieselbe Probe mehrfach angegeben", "Missing file (both channels _Grn and _Red are needed)": "Datei fehlt (beide Kanäle _Grn und _Red werden benötigt)", "idat conversion failed for": "idat-Konvertierung fehlgeschlagen für", "No Magibu result (*.magibu_result.rds) in": "Kein Magibu-Ergebnis (*.magibu_result.rds) in", "Annotation file not found": "Annotationsdatei nicht gefunden", "Columns not found in the annotation": "Spalten in der Annotation nicht gefunden", "Dataset folder not found": "Datensatzordner nicht gefunden", "The destination folder already contains a model, choose a new or empty one": "Der Zielordner enthält bereits ein Modell, bitte einen neuen oder leeren wählen", "No file or folder given.": "Keine Datei und kein Ordner angegeben.", "No idat or tsv file found in the given paths.": "In den angegebenen Pfaden wurde keine idat- oder tsv-Datei gefunden.", "The model has no annotated classes: distances cannot be computed.": "Das Modell hat keine annotierten Klassen: Abstände können nicht berechnet werden.", "The annotation needs at least two columns: sample identifier and class.": "Die Annotation braucht mindestens zwei Spalten: Probenkennung und Klasse.", "Search:": "Suchen:", "Previous": "Zurück", "Next": "Weiter", "No data available in table": "Keine Daten vorhanden", "No matching records found": "Keine passenden Einträge gefunden" };
  Object.keys(DE).forEach(function (k) { if (D[k]) D[k][4] = DE[k]; });


  // added in 0.3: nanopore input (fr, es, it, zh, de)
  var D2 = { "Illumina array type (idat only)": [ "Type de puce Illumina (idat uniquement)", "Tipo de array Illumina (solo idat)", "Tipo di array Illumina (solo idat)", "Illumina 芯片类型（仅 idat）", "Illumina-Array-Typ (nur idat)" ], "idat (_Grn or _Red, also .gz), nanopore BAM (with .bai) or tsv": [ "idat (_Grn ou _Red, aussi .gz), BAM nanopore (avec .bai) ou tsv", "idat (_Grn o _Red, también .gz), BAM nanopore (con .bai) o tsv", "idat (_Grn o _Red, anche .gz), BAM nanopore (con .bai) oppure tsv", "idat（_Grn 或 _Red，也可 .gz）、纳米孔 BAM（需 .bai）或 tsv", "idat (_Grn oder _Red, auch .gz), Nanopore-BAM (mit .bai) oder tsv" ], "Nanopore options (BAM)": [ "Options nanopore (BAM)", "Opciones nanopore (BAM)", "Opzioni nanopore (BAM)", "纳米孔选项（BAM）", "Nanopore-Optionen (BAM)" ], "Minimum coverage per CpG": [ "Couverture minimale par CpG", "Cobertura mínima por CpG", "Copertura minima per CpG", "每个 CpG 的最低覆盖度", "Mindestabdeckung pro CpG" ], "Methylation": [ "Méthylation", "Metilación", "Metilazione", "甲基化", "Methylierung" ], "5mC only": [ "5mC seule", "Solo 5mC", "Solo 5mC", "仅 5mC", "Nur 5mC" ], "Unrecognised file (expected *_Grn.idat / *_Red.idat, also .gz, nanopore BAM or .tsv)": [ "Fichier non reconnu (attendu *_Grn.idat / *_Red.idat, aussi .gz, BAM nanopore ou .tsv)", "Archivo no reconocido (se esperaba *_Grn.idat / *_Red.idat, también .gz, BAM nanopore o .tsv)", "File non riconosciuto (attesi *_Grn.idat / *_Red.idat, anche .gz, BAM nanopore oppure .tsv)", "无法识别的文件（应为 *_Grn.idat / *_Red.idat，也可 .gz，纳米孔 BAM 或 .tsv）", "Unbekannte Datei (erwartet *_Grn.idat / *_Red.idat, auch .gz, Nanopore-BAM oder .tsv)" ]
  , "BAM index (.bai) not found, create it with samtools index": [ "Index BAM (.bai) introuvable, créez-le avec samtools index", "Índice BAM (.bai) no encontrado, créelo con samtools index", "Indice del BAM (.bai) non trovato, crearlo con samtools index", "未找到 BAM 索引（.bai），请用 samtools index 创建", "BAM-Index (.bai) nicht gefunden, mit samtools index erstellen" ]
  , "Reading BAM files is not available on this system (Windows and Linux only)": [ "La lecture des BAM n'est pas disponible sur ce système (Windows et Linux uniquement)", "La lectura de BAM no está disponible en este sistema (solo Windows y Linux)", "La lettura dei BAM non è disponibile su questo sistema (solo Windows e Linux)", "此系统不支持读取 BAM（仅限 Windows 和 Linux）", "BAM-Dateien können auf diesem System nicht gelesen werden (nur Windows und Linux)" ]
  , "BAM reading failed": [ "Échec de la lecture du BAM", "Error al leer el BAM", "Lettura del BAM non riuscita", "读取 BAM 失败", "Lesen des BAM fehlgeschlagen" ], "No idat, BAM or tsv file found in the given paths.": [ "Aucun fichier idat, BAM ou tsv trouvé dans les chemins indiqués.", "No se encontró ningún archivo idat, BAM o tsv en las rutas indicadas.", "Nessun file idat, BAM o tsv trovato nei percorsi indicati.", "在指定路径中未找到 idat、BAM 或 tsv 文件。", "In den angegebenen Pfaden wurde keine idat-, BAM- oder tsv-Datei gefunden." ], "Unrecognised methylation file (expected Magibu tsv or entropy tsv)": [ "Fichier de méthylation non reconnu (attendu tsv Magibu ou tsv entropy)", "Archivo de metilación no reconocido (se esperaba tsv de Magibu o tsv entropy)", "File di metilazione non riconosciuto (attesi tsv Magibu o tsv entropy)", "无法识别的甲基化文件（应为 Magibu tsv 或 entropy tsv）", "Unbekannte Methylierungsdatei (erwartet Magibu-tsv oder entropy-tsv)" ], "Empty file": [ "Fichier vide", "Archivo vacío", "File vuoto", "空文件", "Leere Datei" ], "No CpG near the model positions": [ "Aucun CpG proche des positions du modèle", "Ningún CpG cerca de las posiciones del modelo", "Nessun CpG vicino alle posizioni del modello", "模型位置附近没有 CpG", "Kein CpG in der Nähe der Modellpositionen" ], "Building a model from modkit files is not supported, convert them to Magibu tsv first": [ "Créer un modèle à partir de fichiers modkit n'est pas possible, convertissez-les d'abord en tsv Magibu", "No se puede crear un modelo a partir de archivos modkit, conviértalos antes a tsv de Magibu", "Non si può creare un modello da file modkit, convertirli prima in tsv Magibu", "不支持用 modkit 文件构建模型，请先转换为 Magibu tsv", "Ein Modell aus modkit-Dateien wird nicht unterstützt, zuerst in Magibu-tsv umwandeln" ] };
  Object.keys(D2).forEach(function (k) { D[k] = D2[k]; });

  // texts with variable parts: {1}, {2}... are the captured groups
  var P = [
    [/^Showing (\d+) to (\d+) of (\d+) entries$/,
     ["Lignes {1} à {2} sur {3}", "Mostrando {1} a {2} de {3} filas", "Righe da {1} a {2} di {3}", "显示第 {1} 至 {2} 项，共 {3} 项", "Zeilen {1} bis {2} von {3}"]],
    [/^\(filtered from (\d+) total entries\)$/,
     ["(filtré parmi {1} lignes)", "(filtrado de {1} filas)", "(filtrate da {1} righe)", "（从 {1} 项中筛选）", "(gefiltert aus {1} Zeilen)"]],
    [/^Analysis '(.+)': (\d+) samples\. Folder: (.+?)( \([\d.]+ min\))?$/,
     ["Analyse '{1}' : {2} échantillons. Dossier : {3}{4}", "Análisis '{1}': {2} muestras. Carpeta: {3}{4}",
      "Analisi '{1}': {2} campioni. Cartella: {3}{4}", "分析“{1}”：{2} 个样本。文件夹：{3}{4}", "Analyse '{1}': {2} Proben. Ordner: {3}{4}"]],
    [/^Model built in (.+) \(([\d.]+) min\)\. It is already selected in the Prediction tab\.$/,
     ["Modèle créé dans {1} ({2} min). Il est déjà sélectionné dans l'onglet Prédiction.",
      "Modelo creado en {1} ({2} min). Ya está seleccionado en la pestaña Predicción.",
      "Modello creato in {1} ({2} min). È già selezionato nella scheda Predizione.",
      "模型已创建于 {1}（{2} 分钟）。已在“预测”选项卡中选中。",
      "Modell erstellt in {1} ({2} min). Es ist im Reiter Vorhersage bereits ausgewählt."]],
    [/^Few model CpGs measured \(below 50%\): (.+)\. Their position on the map is unreliable\.$/,
     ["Peu de CpG du modèle mesurés (moins de 50 %) : {1}. Leur position sur la carte n'est pas fiable.",
      "Pocos CpG del modelo medidos (menos del 50%): {1}. Su posición en el mapa no es fiable.",
      "Pochi CpG del modello misurati (meno del 50%): {1}. La loro posizione sulla mappa non è affidabile.",
      "测得的模型 CpG 较少（低于 50%）：{1}。其在图上的位置不可靠。",
      "Wenige Modell-CpGs gemessen (unter 50 %): {1}. Ihre Position auf der Karte ist unzuverlässig."]],
    [/^At least 15 samples are needed to build a model \(found (\d+)\)\.$/,
     ["Il faut au moins 15 échantillons pour créer un modèle (trouvés : {1}).",
      "Se necesitan al menos 15 muestras para crear un modelo (encontradas: {1}).",
      "Servono almeno 15 campioni per creare un modello (trovati {1}).",
      "构建模型至少需要 15 个样本（找到 {1} 个）。",
      "Für ein Modell werden mindestens 15 Proben benötigt (gefunden: {1})."]]
  ];
  var ERR = ["Erreur : ", "Error: ", "Errore: ", "错误：", "Fehler: "];

  var IDX = { fr: 0, es: 1, it: 2, zh: 3, de: 4 };
  var lang = "en";
  var TXT = new WeakMap(), ATT = new WeakMap();
  // never translated: plots, logs, code, typed text, table cells (class names, IDs)
  var SKIP = "svg, script, style, textarea, pre, code, .plotly, .mg-flags, .shiny-plot-output, .dataTable tbody td:not(.dataTables_empty)";

  function core(t) {
    var i = IDX[lang], e = D[t];
    if (e) return e[i] || t;
    if (t.indexOf("Error: ") === 0) { var r = t.slice(7); return ERR[i] + (core(r) || r); }
    for (var k = 0; k < P.length; k++) {
      var m = t.match(P[k][0]);
      if (m) return P[k][1][i].replace(/\{(\d)\}/g, function (_, g) { return m[+g] || ""; });
    }
    var c = t.indexOf(": ");
    if (c > 0 && D[t.slice(0, c)]) return D[t.slice(0, c)][i] + (lang === "zh" ? "：" : ": ") + t.slice(c + 2);
    return null;
  }
  function tr(s) {
    if (lang === "en" || !s) return s;
    var t = s.trim();
    if (!t) return s;
    var out = core(t);
    return out === null ? s : s.replace(t, function () { return out; });
  }
  function skip(el) { return !el || (el.closest && el.closest(SKIP)); }

  function doText(n) {
    if (skip(n.parentElement)) return;
    var cur = n.nodeValue, r = TXT.get(n);
    var src = (r && r.out === cur) ? r.src : cur;
    var out = tr(src);
    TXT.set(n, { src: src, out: out });
    if (out !== cur) n.nodeValue = out;
  }
  function doAttr(el, a) {
    if (!el.hasAttribute(a) || skip(el.parentElement)) return;
    var cur = el.getAttribute(a), m = ATT.get(el) || {}, r = m[a];
    var src = (r && r.out === cur) ? r.src : cur;
    var out = tr(src);
    m[a] = { src: src, out: out }; ATT.set(el, m);
    if (out !== cur) el.setAttribute(a, out);
  }
  function walk(root) {
    if (!root) return;
    if (root.nodeType === 3) { doText(root); return; }
    if (root.nodeType !== 1 || skip(root)) return;
    doAttr(root, "placeholder");
    root.querySelectorAll("[placeholder]").forEach(function (el) { doAttr(el, "placeholder"); });
    var w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT), n;
    while ((n = w.nextNode())) doText(n);
  }

  function setLang(l, save) {
    if (l !== "en" && !(l in IDX)) l = "en";
    lang = l;
    document.documentElement.lang = l === "zh" ? "zh-CN" : l;
    document.querySelectorAll(".mg-flag").forEach(function (b) { b.classList.toggle("active", b.getAttribute("data-lang") === l); });
    try { localStorage.setItem("magibu-lang", l); } catch (e) {}
    walk(document.body);
    if (save && window.Shiny && Shiny.setInputValue) Shiny.setInputValue("mg_lang", l, { priority: "event" });
  }

  new MutationObserver(function (list) {
    if (lang === "en") return;
    list.forEach(function (m) {
      if (m.type === "characterData") doText(m.target);
      else if (m.type === "attributes") doAttr(m.target, m.attributeName);
      else m.addedNodes.forEach(walk);
    });
  }).observe(document.documentElement, { subtree: true, childList: true, characterData: true, attributes: true, attributeFilter: ["placeholder"] });

  $(function () {
    var saved = null;
    try { saved = localStorage.getItem("magibu-lang"); } catch (e) {}
    setLang(saved || "en", false);
    if (window.Shiny) Shiny.addCustomMessageHandler("mg-lang", function (l) { setLang(l, false); });
  });
  $(document).on("click", ".mg-flag", function () { setLang($(this).data("lang"), true); });
  window.magibuLang = setLang;
})();
