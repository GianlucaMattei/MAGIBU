# ==============================================================================
# Magibu App - interface (Shiny). Started by Electron (executable) or with
#   Rscript -e "shiny::runApp('app', launch.browser = TRUE)"
# Interface text is written in English; www/i18n.js translates it in the browser
# (flags at the top right). Plots and output files stay in English.
# ==============================================================================
suppressPackageStartupMessages({ library(shiny); library(bslib) })

app_dir <- normalizePath(getwd(), winslash = "/")
source(file.path(app_dir, "R", "magibu_lib.R"))
magibu_init(app_dir)

models_root <- Sys.getenv("MAGIBU_MODELS", file.path(dirname(app_dir), "models"))
default_model <- {
  d <- if (dir.exists(models_root)) list.dirs(models_root, recursive = FALSE) else character(0)
  d <- d[vapply(d, function(x) length(list.files(x, pattern = "\\.uwot$")) == 1, logical(1))]
  if (length(d) > 0) normalizePath(d[1], winslash = "/") else ""
}
default_out <- file.path(path.expand("~"), "Magibu_results", format(Sys.Date(), "%Y-%m-%d"))
max_cores <- max(1L, min(24L, parallel::detectCores() - 1L, na.rm = TRUE))

# chosen language, remembered between sessions
lang_file <- file.path(tools::R_user_dir("Magibu", "config"), "language.txt")
langs <- c("en", "de", "fr", "es", "it", "zh")
read_lang <- function() { l <- tryCatch(readLines(lang_file, 1, warn = FALSE), error = function(e) "en"); if (l %in% langs) l else "en" }
# day / night mode, also remembered
mode_file <- file.path(dirname(lang_file), "mode.txt")
read_mode <- function() { m <- tryCatch(readLines(mode_file, 1, warn = FALSE), error = function(e) "light"); if (identical(m, "dark")) "dark" else "light" }
save_setting <- function(value, file) try({ dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE); writeLines(value, file) }, silent = TRUE)

# path field + Browse button (native dialog when running inside the executable)
path_input <- function(id, label, value = "", mode = "dir", rows = NULL, placeholder = NULL) {
  field <- if (is.null(rows)) textInput(id, label, value, width = "100%", placeholder = placeholder)
           else textAreaInput(id, label, value, width = "100%", rows = rows, placeholder = placeholder)
  div(class = "mg-path", field,
      tags$button(type = "button", class = "btn btn-outline-secondary btn-sm mg-browse",
                  `data-target` = id, `data-mode` = mode, "Browse..."))
}

array_choices <- c("Automatic" = "auto", "EPIC v2" = "EpicV2", "EPIC v1 (850k)" = "EpicV1", "450k" = "Infinium")

status_box <- function(st) {
  if (is.null(st)) return(NULL)
  cls <- switch(st$type, ok = "alert-success", error = "alert-danger", run = "alert-info", warn = "alert-warning", "alert-secondary")
  div(class = paste("alert", cls, "mg-status"), st$msg)
}

flag_bar <- div(class = "mg-flags",
  lapply(list(c("en", "English"), c("de", "Deutsch"), c("fr", "Français"), c("es", "Español"), c("it", "Italiano"), c("zh", "中文")),
         function(l) tags$button(type = "button", class = "mg-flag", `data-lang` = l[1], title = l[2],
                                 tags$img(src = paste0("flags/", l[1], ".svg"), alt = l[2]))))

# colour theme: "v2" (light, petrol blue, with night mode) or the old "flatly", chosen with MAGIBU_THEME
ui_theme <- Sys.getenv("MAGIBU_THEME", "v2")
app_theme <- if (ui_theme == "v2") {
  bs_theme(version = 5, bg = "#FFFFFF", fg = "#1F2937", primary = "#0E7490", secondary = "#64748B",
           success = "#15803D", danger = "#B91C1C", info = "#0E7490",
           base_font = font_collection("Segoe UI", "system-ui", "-apple-system", "Roboto", "Helvetica Neue", "Arial", "sans-serif"))
} else bs_theme(version = 5, bootswatch = "flatly")

ui <- page_navbar(
  title = "Magibu", id = "tab", theme = app_theme,
  header = tags$head(tags$link(rel = "stylesheet", href = "magibu.css"),
                     if (ui_theme == "v2") tags$link(rel = "stylesheet", href = "magibu_v2.css"),
                     tags$script(src = "i18n.js"), tags$script(src = "magibu.js")),

  # ---------------------------------------------------------------- PREDICTION
  nav_panel("Prediction", value = "predict", layout_sidebar(
    sidebar = sidebar(width = 400,
      h5("1. Model"),
      path_input("pred_model", "Model folder", default_model, "dir"),
      h5("2. Samples"),
      radioButtons("pred_mode", NULL, inline = FALSE,
                   c("Single file" = "file", "All files in a folder" = "dir", "Selected files" = "files")),
      path_input("pred_paths", "Paths (one per line)", "", "pred", rows = 4,
                 placeholder = "idat (_Grn or _Red, also .gz), nanopore BAM (with .bai) or tsv"),
      selectInput("pred_array", "Illumina array type (idat only)", array_choices, width = "100%"),
      accordion(open = FALSE, accordion_panel("Nanopore options (BAM)",
        numericInput("pred_mincov", "Minimum coverage per CpG", 1, min = 1, step = 1, width = "100%"),
        radioButtons("pred_5hmc", "Methylation", c("5mC only" = "m", "5mC + 5hmC" = "mh"), inline = TRUE))),
      h5("3. Results"),
      path_input("pred_out", "Output folder", default_out, "dir"),
      textInput("pred_name", "Analysis name", "magibu", width = "100%"),
      numericInput("pred_cores", "Cores to use", magibu_default_cores(), min = 1, max = max_cores, width = "100%"),
      actionButton("pred_run", "Run prediction", class = "btn-primary btn-lg w-100"),
      actionButton("job_stop", "Stop", class = "btn-outline-danger btn-sm w-100"),
      hr(),
      path_input("prev_path", "Reopen a previous analysis (folder)", "", "dir"),
      actionButton("prev_load", "Load", class = "btn-outline-secondary btn-sm w-100")
    ),
    uiOutput("pred_status"),
    navset_card_tab(id = "pred_tabs",
      nav_panel("Summary", value = "summary",
        div(class = "mg-toolbar", downloadButton("dl_xlsx", "Download distances (Excel)", class = "btn-sm"),
            tags$button(type = "button", class = "btn btn-outline-secondary btn-sm mg-open", "Open results folder")),
        DT::DTOutput("tbl_summary")),
      nav_panel("Sample names", value = "names",
        p(class = "text-muted small", "Give each sample a name (for example Case 1). Names appear in tables, plots and in the Excel file; the file name stays as sample ID."),
        div(class = "mg-toolbar",
            actionButton("nm_case", "Number as Case 1, 2, 3...", class = "btn-outline-primary btn-sm"),
            actionButton("nm_reset", "Clear names", class = "btn-outline-secondary btn-sm"),
            actionButton("nm_save", "Save names", class = "btn-primary btn-sm")),
        uiOutput("names_table")),
      nav_panel("Distances by class", value = "dist",
        selectInput("dist_sample", "Sample", choices = NULL, width = "400px"),
        DT::DTOutput("tbl_dist")),
      nav_panel("Classic plot", value = "classic",
        div(class = "mg-toolbar", selectInput("classic_sample", "Sample", choices = NULL, width = "400px"),
            downloadButton("dl_pdf", "Download PDF (all samples)", class = "btn-sm")),
        plotOutput("plot_classic", height = "640px")),
      nav_panel("Log", value = "log", verbatimTextOutput("log_pred"))
    )
  )),

  # ---------------------------------------------------------------- GLOBAL MAP
  nav_panel("Global map", value = "map", layout_sidebar(
    sidebar = sidebar(width = 320,
      radioButtons("gm_color", "Colour by", c("Class" = "class", "Family (text before the comma)" = "family")),
      checkboxInput("gm_labels", "Show sample names", TRUE),
      div(class = "mg-toolbar",
          actionButton("gm_samples", "Zoom on samples", class = "btn-outline-primary btn-sm"),
          actionButton("gm_all", "Whole map", class = "btn-outline-secondary btn-sm")),
      div(class = "mg-help",
          p("Mouse wheel: zoom in and out."),
          p("Drag: zoom on an area."),
          p("Double click: back to the whole map."),
          p("Click a class in the legend to hide it, double click to show only that one."),
          p("Pass over a point to see its class.")),
      p(class = "text-muted small", "Download the current view:"),
      div(class = "mg-toolbar", downloadButton("dl_map_pdf", "PDF", class = "btn-sm"),
          downloadButton("dl_map_png", "PNG", class = "btn-sm"))
    ),
    uiOutput("map_status"),
    card(full_screen = TRUE, plotly::plotlyOutput("gmap", height = "780px"))
  )),

  # ---------------------------------------------------------------- CLASS ZOOM
  nav_panel("Class zoom", value = "zoom", layout_sidebar(
    sidebar = sidebar(width = 400,
      p(class = "text-muted small", "The view widens until it contains every reference sample of the chosen classes and the chosen samples."),
      selectizeInput("zoom_samples", "Samples to show", choices = NULL, multiple = TRUE, width = "100%",
                     options = list(plugins = list("remove_button"))),
      selectizeInput("zoom_classes", "Classes to highlight", choices = NULL, multiple = TRUE, width = "100%",
                     options = list(plugins = list("remove_button"), placeholder = "choose one or more classes")),
      div(class = "mg-toolbar",
          actionButton("zoom_near", "Add the closest", class = "btn-outline-primary btn-sm"),
          actionButton("zoom_family", "Whole family", class = "btn-outline-primary btn-sm"),
          actionButton("zoom_clear", "Clear", class = "btn-outline-secondary btn-sm")),
      checkboxInput("zoom_other", "Show other reference samples in grey", TRUE),
      checkboxInput("zoom_lab_cls", "Class labels on the plot", FALSE),
      checkboxInput("zoom_lab_smp", "Sample labels on the plot", TRUE),
      p(class = "text-muted small", "Sample names are set in Prediction > Sample names."),
      textInput("zoom_title", "Title", "", width = "100%"),
      div(class = "mg-toolbar", downloadButton("dl_zoom_pdf", "PDF", class = "btn-sm"),
          downloadButton("dl_zoom_png", "PNG", class = "btn-sm"))
    ),
    uiOutput("zoom_status"),
    card(full_screen = TRUE, plotOutput("plot_zoom", height = "760px"))
  )),

  # ---------------------------------------------------------------- BUILD MODEL
  nav_panel("Build model", value = "train", layout_sidebar(
    sidebar = sidebar(width = 400,
      path_input("tr_data", "Dataset folder (idat or tsv)", "", "dir"),
      path_input("tr_annot", "Annotation file (tsv, csv)", "", "file"),
      layout_columns(col_widths = c(6, 6),
        selectInput("tr_idcol", "Sample column", choices = NULL),
        selectInput("tr_clscol", "Class column", choices = NULL)),
      selectInput("tr_array", "Array type", array_choices, width = "100%"),
      path_input("tr_out", "Folder where the model is saved", "", "dir"),
      textInput("tr_name", "Model name", "model", width = "100%"),
      accordion(open = FALSE, accordion_panel("Advanced parameters",
        numericInput("tr_nvars", "Most variable CpGs", 50000, min = 1000, step = 1000),
        numericInput("tr_pcs", "PCA components", 50, min = 2),
        numericInput("tr_nn", "UMAP neighbours", 30, min = 2),
        numericInput("tr_mindist", "UMAP minimum distance", 0.3, min = 0, max = 1, step = 0.05),
        numericInput("tr_cores", "Cores to use", magibu_default_cores(), min = 1, max = max_cores))),
      actionButton("tr_run", "Build model", class = "btn-primary btn-lg w-100"),
      actionButton("job_stop2", "Stop", class = "btn-outline-danger btn-sm w-100")
    ),
    uiOutput("tr_status"),
    navset_card_tab(
      nav_panel("Model map", plotOutput("plot_ref", height = "700px")),
      nav_panel("Log", verbatimTextOutput("log_tr")))
  )),

  nav_panel("About", value = "info", card(
    h4(paste("Magibu App", MAGIBU_VERSION)),
    p(strong("For research use only."), " Not a diagnostic device: the result does not replace the diagnosis."),
    p("The prediction projects each sample onto the UMAP map of the chosen model and measures its distances from the reference samples."),
    p("The UMAP projection is stochastic: running the same sample again can move it slightly (about 0.2 UMAP units), as in the original Magibu pipeline."),
    p("Third party software included (R and R packages, Electron): see the THIRD_PARTY_LICENSES file in the program folder."))),

  nav_spacer(),
  nav_item(div(class = "mg-topbar", if (ui_theme == "v2") input_dark_mode(id = "dark_mode", mode = "light"), flag_bar))
)

server <- function(input, output, session) {
  rv <- reactiveValues(result = NULL, job = NULL, log = character(0), st_pred = NULL, st_tr = NULL,
                       ref_model = NULL, view = NULL)

  # ---- language: saved choice sent to the page; a flag click is saved for next time
  session$sendCustomMessage("mg-lang", read_lang())
  observeEvent(input$mg_lang, if (input$mg_lang %in% langs) save_setting(input$mg_lang, lang_file), ignoreInit = TRUE)
  # night mode: saved choice applied at start, every switch saved
  if (ui_theme == "v2" && read_mode() == "dark") toggle_dark_mode("dark")
  observeEvent(input$dark_mode, save_setting(input$dark_mode, mode_file), ignoreInit = TRUE)

  # ---- tests / reopening: ?tab=zoom, MAGIBU_AUTOLOAD=<results folder>
  observe({
    q <- parseQueryString(isolate(session$clientData$url_search))
    if (!is.null(q$tab)) nav_select("tab", q$tab)
    if (!is.null(q$sub)) nav_select("pred_tabs", q$sub)
    auto <- Sys.getenv("MAGIBU_AUTOLOAD", "")
    if (nzchar(auto)) load_result(auto)
  }) |> bindEvent(TRUE, once = TRUE)

  load_result <- function(path) {
    res <- tryCatch(readRDS(magibu_find_result(path)), error = function(e) e)
    if (inherits(res, "error")) { rv$st_pred <- list(type = "error", msg = conditionMessage(res)); return(invisible()) }
    rv$result <- res
    rv$st_pred <- list(type = "ok", msg = paste0("Analysis '", res$name, "': ", length(res$samples), " samples. Folder: ", res$out_dir))
  }
  observeEvent(input$prev_load, load_result(trimws(input$prev_path)))

  # ---- background jobs (separate R: the interface stays responsive)
  start_job <- function(kind, fun, args) {
    if (!is.null(rv$job)) { showNotification("A job is already running.", type = "warning"); return(invisible()) }
    logf <- tempfile(fileext = ".log")
    proc <- callr::r_bg(function(app_dir, fun, args) {
      source(file.path(app_dir, "R", "magibu_lib.R")); magibu_init(app_dir); do.call(fun, args)
    }, args = list(app_dir, fun, args), stdout = logf, stderr = "2>&1", supervise = TRUE)
    rv$log <- character(0)
    rv$job <- list(proc = proc, logf = logf, kind = kind, t0 = Sys.time())
    st <- list(type = "run", msg = "Job running... (details in the Log tab)")
    if (kind == "predict") rv$st_pred <- st else rv$st_tr <- st
  }

  observe({
    j <- rv$job; req(j)
    invalidateLater(800)
    if (file.exists(j$logf)) rv$log <- utils::tail(readLines(j$logf, warn = FALSE), 300)
    if (j$proc$is_alive()) return()
    res <- tryCatch(j$proc$get_result(), error = function(e) e)
    rv$job <- NULL
    mins <- round(as.numeric(difftime(Sys.time(), j$t0, units = "mins")), 1)
    if (inherits(res, "error")) {
      msg <- if (!is.null(res$parent)) conditionMessage(res$parent) else conditionMessage(res)
      st <- list(type = "error", msg = paste("Error:", msg))
      if (j$kind == "predict") rv$st_pred <- st else rv$st_tr <- st
    } else if (j$kind == "predict") {
      load_result(res)
      rv$st_pred$msg <- paste0(rv$st_pred$msg, " (", mins, " min)")
    } else {
      rv$ref_model <- res
      updateTextInput(session, "pred_model", value = res)
      rv$st_tr <- list(type = "ok", msg = paste0("Model built in ", res, " (", mins, " min). It is already selected in the Prediction tab."))
    }
  })

  stop_job <- function() {
    j <- rv$job
    if (is.null(j)) return()
    j$proc$kill_tree()
    rv$job <- NULL
    st <- list(type = "error", msg = "Job stopped.")
    if (j$kind == "predict") rv$st_pred <- st else rv$st_tr <- st
  }
  observeEvent(input$job_stop, stop_job()); observeEvent(input$job_stop2, stop_job())
  session$onSessionEnded(function() { j <- isolate(rv$job); if (!is.null(j)) try(j$proc$kill_tree(), silent = TRUE) })

  need <- function(cond, msg, slot) {
    if (!cond) { rv[[slot]] <- list(type = "error", msg = msg); FALSE } else TRUE
  }

  # ---------------------------------------------------------------- PREDICTION
  observeEvent(input$pred_run, {
    paths <- trimws(input$pred_paths)
    ok <- need(nzchar(trimws(input$pred_model)), "Choose the model folder.", "st_pred") &&
      need(nzchar(paths), "Give at least one sample file or folder.", "st_pred") &&
      need(nzchar(trimws(input$pred_out)), "Choose the output folder.", "st_pred")
    if (!ok) return()
    chk <- tryCatch({ magibu_model_files(trimws(input$pred_model)); magibu_collect_inputs(paths) }, error = function(e) e)
    if (inherits(chk, "error")) { rv$st_pred <- list(type = "error", msg = conditionMessage(chk)); return() }
    start_job("predict", "magibu_run_prediction",
              list(paths = paths, model_dir = trimws(input$pred_model), out_dir = trimws(input$pred_out),
                   name = input$pred_name, array_type = input$pred_array, n_cores = input$pred_cores,
                   min_cov = input$pred_mincov, include_5hmc = identical(input$pred_5hmc, "mh")))
  })

  output$pred_status <- renderUI(tagList(status_box(rv$st_pred), status_box(low_cov())))
  # samples with few model CpGs measured (e.g. short nanopore runs): position on the map unreliable
  low_cov <- reactive({
    sm <- rv$result$distances$Summary
    if (is.null(sm) || !"Model_CpG_coverage_pct" %in% names(sm)) return(NULL)
    low <- sm[!is.na(Model_CpG_coverage_pct) & Model_CpG_coverage_pct < 50]
    if (nrow(low) == 0) return(NULL)
    list(type = "warn", msg = paste0("Few model CpGs measured (below 50%): ",
                                     paste0(low$Sample, " ", low$Model_CpG_coverage_pct, "%", collapse = ", "),
                                     ". Their position on the map is unreliable."))
  })
  output$log_pred <- renderText(paste(rv$log, collapse = "\n"))
  output$log_tr <- renderText(paste(rv$log, collapse = "\n"))

  # id -> display name
  snames <- reactive({ req(rv$result); magibu_sample_names(rv$result) })
  named_choices <- function(nm) stats::setNames(names(nm), ifelse(nm == names(nm), nm, paste0(nm, "  (", names(nm), ")")))

  observeEvent(snames(), {
    ch <- named_choices(snames())
    keep <- function(x) if (!is.null(x) && x %in% ch) x else ch[[1]]
    updateSelectInput(session, "dist_sample", choices = ch, selected = keep(isolate(input$dist_sample)))
    updateSelectInput(session, "classic_sample", choices = ch, selected = keep(isolate(input$classic_sample)))
    sel <- isolate(input$zoom_samples); sel <- if (length(sel) > 0) intersect(sel, ch) else unname(ch)
    updateSelectizeInput(session, "zoom_samples", choices = ch, selected = sel, server = FALSE)
  })
  observeEvent(rv$result, session$sendCustomMessage("mg-outdir", rv$result$out_dir))

  dt_opts <- list(scrollX = TRUE, pageLength = 15, dom = "frtip")
  output$tbl_summary <- DT::renderDT({
    req(rv$result)
    DT::datatable(.mg_with_names(rv$result$distances$Summary, snames()), rownames = FALSE, options = dt_opts, selection = "none")
  })
  output$tbl_dist <- DT::renderDT({
    req(rv$result, input$dist_sample)
    d <- rv$result$distances$Distances
    DT::datatable(d[d$Sample == input$dist_sample, -1], rownames = FALSE, options = dt_opts, selection = "none")
  })
  output$plot_classic <- renderPlot({
    req(rv$result, input$classic_sample)
    suppressWarnings(plot_umap_zoom(as.data.frame(rv$result$metadata), input$classic_sample, cases_files = rv$result$samples,
                                    circle_radius = 0.3, distance_multiplier = 1, n_top_groups = 10))
  }, res = 72)
  output$dl_xlsx <- downloadHandler(
    filename = function() basename(rv$result$files$xlsx),
    content = function(file) file.copy(rv$result$files$xlsx, file))
  output$dl_pdf <- downloadHandler(
    filename = function() basename(rv$result$files$pdf),
    content = function(file) file.copy(rv$result$files$pdf, file))

  # ---- sample names: one field per sample (rebuilt only when another analysis is loaded)
  output$names_table <- renderUI({
    res <- rv$result
    if (is.null(res)) return(p(class = "text-muted", "Run or load an analysis first."))
    nm <- isolate(snames())
    tags$table(class = "table table-sm mg-names",
      tags$thead(tags$tr(tags$th("Sample ID"), tags$th("Name"), tags$th("Closest class"))),
      tags$tbody(lapply(seq_along(res$samples), function(i) {
        s <- res$samples[i]
        tags$tr(tags$td(tags$code(s)),
                tags$td(textInput(paste0("nm_", i), NULL, if (nm[[s]] == s) "" else nm[[s]], width = "100%", placeholder = s)),
                tags$td(res$distances$Summary$Group1_Class[res$distances$Summary$Sample == s]))
      })))
  })
  observeEvent(input$nm_case, {
    req(rv$result)
    for (i in seq_along(rv$result$samples)) updateTextInput(session, paste0("nm_", i), value = paste("Case", i))
  })
  observeEvent(input$nm_reset, {
    req(rv$result)
    for (i in seq_along(rv$result$samples)) updateTextInput(session, paste0("nm_", i), value = "")
  })
  observeEvent(input$nm_save, {
    res <- rv$result; req(res)
    lab <- stats::setNames(vapply(seq_along(res$samples), function(i) { v <- input[[paste0("nm_", i)]]; if (is.null(v)) "" else v }, ""),
                           res$samples)
    out <- withProgress(message = "Saving names...", value = 0.3,
                        tryCatch(magibu_apply_names(res, lab), error = function(e) e))
    if (inherits(out, "error")) { rv$st_pred <- list(type = "error", msg = paste("Error:", conditionMessage(out))); return() }
    rv$result <- out
    rv$st_pred <- list(type = "ok", msg = "Names saved: Excel file and classic PDF updated.")
  })

  # ---------------------------------------------------------------- GLOBAL MAP
  zoom_meta <- reactive({
    if (!is.null(rv$result)) return(rv$result$metadata)
    tryCatch(magibu_read_model_meta(trimws(input$pred_model)), error = function(e) NULL)
  })
  map_samples <- reactive(if (is.null(rv$result)) character(0) else rv$result$samples)
  map_names <- reactive(if (is.null(rv$result)) NULL else snames())

  output$map_status <- renderUI({
    if (is.null(zoom_meta())) status_box(list(type = "info", msg = "Choose a valid model in the Prediction tab, or load an analysis."))
    else if (is.null(rv$result)) status_box(list(type = "info", msg = "No analysis loaded: only the model map is shown. Samples appear after a prediction."))
  })

  output$gmap <- plotly::renderPlotly({
    m <- zoom_meta(); req(m)
    ref <- magibu_global_data(m, input$gm_color)
    pal <- magibu_class_colors(ref$group)
    if ("Unlabelled" %in% ref$group) pal["Unlabelled"] <- "#C8C8C8"
    ref[, hover := paste0("<b>", ifelse(is.na(sample.title), "Unlabelled", sample.title), "</b><br>", Sample_ID_Full)]
    ref[, group := factor(group, levels = names(pal))]
    p <- plotly::plot_ly(source = "gmap") |>
      plotly::add_markers(data = ref, x = ~UMAP1, y = ~UMAP2, color = ~group, colors = pal,
                          marker = list(size = 5, opacity = 0.85), text = ~hover, hoverinfo = "text")
    s <- map_samples()
    if (length(s) > 0) {
      mm <- magibu_standard_meta(m)
      fg <- mm[magibu_is_new(mm) & Sample_ID_Full %in% s]
      nm <- map_names()
      fg[, lab := unname(nm[Sample_ID_Full])]
      fg[, hover := paste0("<b>", lab, "</b><br>", Sample_ID_Full)]
      p <- p |> plotly::add_trace(data = fg, x = ~UMAP1, y = ~UMAP2, type = "scatter",
                                  mode = if (isTRUE(input$gm_labels)) "markers+text" else "markers",
                                  text = ~lab, textposition = "top right", hovertext = ~hover, hoverinfo = "text",
                                  textfont = list(size = 12, color = "black"), name = "Samples", inherit = FALSE,
                                  marker = list(symbol = "asterisk-open", size = 10, color = "black",
                                                line = list(color = "black", width = 1.8)))
    }
    rv$view <- NULL
    p |>
      plotly::layout(xaxis = list(title = "UMAP 1", zeroline = FALSE, showgrid = FALSE),
                     yaxis = list(title = "UMAP 2", zeroline = FALSE, showgrid = FALSE, scaleanchor = "x"),
                     legend = list(font = list(size = 10), itemsizing = "constant"),
                     dragmode = "zoom", plot_bgcolor = "white", margin = list(t = 20)) |>
      plotly::config(scrollZoom = TRUE, displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"),
                     toImageButtonOptions = list(format = "png", scale = 3, filename = "magibu_map")) |>
      plotly::event_register("plotly_relayout")
  })

  # current view (used by the downloads)
  observeEvent(plotly::event_data("plotly_relayout", source = "gmap"), {
    ev <- plotly::event_data("plotly_relayout", source = "gmap")
    if (isTRUE(ev[["xaxis.autorange"]])) { rv$view <- NULL; return() }
    x <- c(ev[["xaxis.range[0]"]], ev[["xaxis.range[1]"]]); y <- c(ev[["yaxis.range[0]"]], ev[["yaxis.range[1]"]])
    if (is.null(x) && !is.null(ev[["xaxis.range"]])) x <- unlist(ev[["xaxis.range"]])
    if (is.null(y) && !is.null(ev[["yaxis.range"]])) y <- unlist(ev[["yaxis.range"]])
    if (length(x) == 2 && length(y) == 2) rv$view <- list(x = sort(as.numeric(x)), y = sort(as.numeric(y)))
  })
  observeEvent(input$gm_samples, {
    m <- zoom_meta(); s <- map_samples()
    if (length(s) == 0) { showNotification("No samples: run or load an analysis first.", type = "warning"); return() }
    mm <- magibu_standard_meta(m); fg <- mm[magibu_is_new(mm) & Sample_ID_Full %in% s]
    pad <- max(1.5, 0.15 * max(diff(range(fg$UMAP1)), diff(range(fg$UMAP2))))
    x <- range(fg$UMAP1) + c(-pad, pad); y <- range(fg$UMAP2) + c(-pad, pad)
    plotly::plotlyProxy("gmap", session) |>
      plotly::plotlyProxyInvoke("relayout", list(xaxis.range = x, yaxis.range = y))
    rv$view <- list(x = x, y = y)
  })
  observeEvent(input$gm_all, {
    plotly::plotlyProxy("gmap", session) |>
      plotly::plotlyProxyInvoke("relayout", list(xaxis.autorange = TRUE, yaxis.autorange = TRUE))
    rv$view <- NULL
  })
  map_plot <- function() {
    v <- rv$view
    magibu_global_plot(zoom_meta(), map_samples(), map_names(), color_by = input$gm_color,
                       xlim = v$x, ylim = v$y, label_samples = isTRUE(input$gm_labels))
  }
  map_file <- function(ext) paste0(if (!is.null(rv$result)) rv$result$name else "magibu", ".map.", ext)
  output$dl_map_pdf <- downloadHandler(filename = function() map_file("pdf"),
    content = function(file) ggplot2::ggsave(file, map_plot(), width = 16, height = 10, device = "pdf"))
  output$dl_map_png <- downloadHandler(filename = function() map_file("png"),
    content = function(file) ggplot2::ggsave(file, map_plot(), width = 16, height = 10, dpi = 200, device = "png"))

  # ---------------------------------------------------------------- CLASS ZOOM
  observeEvent(zoom_meta(), {
    m <- zoom_meta()
    cls <- sort(unique(m$sample.title[!magibu_is_new(m) & !is.na(m$sample.title)]))
    grp <- split(cls, magibu_family(cls))
    grp <- lapply(grp, function(x) stats::setNames(as.list(x), x))
    keep <- intersect(isolate(input$zoom_classes), cls)
    if (length(keep) == 0 && !is.null(rv$result)) keep <- magibu_nearest_classes(rv$result$distances, rv$result$samples, 3)
    updateSelectizeInput(session, "zoom_classes", choices = grp, selected = keep, server = FALSE)
  })
  observeEvent(input$zoom_near, {
    req(rv$result)
    s <- if (length(input$zoom_samples) > 0) input$zoom_samples else rv$result$samples
    updateSelectizeInput(session, "zoom_classes", selected = union(input$zoom_classes, magibu_nearest_classes(rv$result$distances, s, 3)))
  })
  observeEvent(input$zoom_family, {
    m <- zoom_meta(); req(m, length(input$zoom_classes) > 0)
    cls <- unique(m$sample.title[!magibu_is_new(m) & !is.na(m$sample.title)])
    updateSelectizeInput(session, "zoom_classes", selected = union(input$zoom_classes, cls[magibu_family(cls) %in% magibu_family(input$zoom_classes)]))
  })
  observeEvent(input$zoom_clear, updateSelectizeInput(session, "zoom_classes", selected = character(0)))

  zoom_plot <- reactive({
    m <- zoom_meta(); req(m)
    magibu_zoom_plot(m, classes = input$zoom_classes, samples = input$zoom_samples,
                     labels = map_names(), show_other = input$zoom_other,
                     label_classes = input$zoom_lab_cls, label_samples = input$zoom_lab_smp,
                     title = if (nzchar(input$zoom_title)) input$zoom_title else NULL)
  })
  output$zoom_status <- renderUI({
    if (is.null(zoom_meta())) status_box(list(type = "info", msg = "Choose a valid model in the Prediction tab, or load an analysis."))
    else if (is.null(rv$result)) status_box(list(type = "info", msg = "No analysis loaded: only the model map is shown. Samples appear after a prediction."))
  })
  output$plot_zoom <- renderPlot(zoom_plot(), res = 80)
  zoom_file <- function(ext) paste0(if (!is.null(rv$result)) rv$result$name else "magibu", ".zoom.", ext)
  output$dl_zoom_pdf <- downloadHandler(filename = function() zoom_file("pdf"),
    content = function(file) ggplot2::ggsave(file, zoom_plot(), width = 14, height = 10, device = "pdf"))
  output$dl_zoom_png <- downloadHandler(filename = function() zoom_file("png"),
    content = function(file) ggplot2::ggsave(file, zoom_plot(), width = 14, height = 10, dpi = 200, device = "png"))

  # ---------------------------------------------------------------- BUILD MODEL
  observeEvent(input$tr_annot, {
    f <- trimws(input$tr_annot)
    cols <- if (nzchar(f) && file.exists(f) && !dir.exists(f)) tryCatch(names(data.table::fread(f, nrows = 0)), error = function(e) NULL)
    updateSelectInput(session, "tr_idcol", choices = cols, selected = cols[1])
    updateSelectInput(session, "tr_clscol", choices = cols, selected = cols[min(2, length(cols))])
  })
  observeEvent(input$tr_run, {
    v <- lapply(list(data = input$tr_data, annot = input$tr_annot, out = input$tr_out), trimws)
    ok <- need(nzchar(v$data) && dir.exists(v$data), "Choose the dataset folder.", "st_tr") &&
      need(nzchar(v$annot) && file.exists(v$annot), "Choose the annotation file.", "st_tr") &&
      need(nzchar(v$out), "Choose the folder where the model is saved.", "st_tr")
    if (!ok) return()
    start_job("train", "magibu_run_training",
              list(dataset_dir = v$data, annotation_file = v$annot, model_dir = v$out, name = input$tr_name,
                   id_col = input$tr_idcol, class_col = input$tr_clscol, array_type = input$tr_array,
                   n_top_vars = input$tr_nvars, pca_components = input$tr_pcs, n_neighbors = input$tr_nn,
                   min_dist = input$tr_mindist, n_cores = input$tr_cores))
  })
  output$tr_status <- renderUI(status_box(rv$st_tr))
  output$plot_ref <- renderPlot({
    req(rv$ref_model)
    magibu_reference_map(magibu_read_model_meta(rv$ref_model), basename(rv$ref_model))
  }, res = 80)
}

shinyApp(ui, server)
