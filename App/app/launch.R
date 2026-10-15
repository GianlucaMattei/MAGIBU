# Starts the app from the Electron shell: port and folder come from the environment (no code on the command line,
# so paths with spaces are safe).
app_dir <- Sys.getenv("MAGIBU_APP_DIR")
port <- as.integer(Sys.getenv("MAGIBU_PORT"))
shiny::runApp(app_dir, port = port, host = "127.0.0.1", launch.browser = FALSE)
