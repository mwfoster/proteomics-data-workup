app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

cached <- data.frame(SourceLabel = "cached", Value = 1, stringsAsFactors = FALSE)
stopifnot(identical(app_env$project_input_table(NULL, cached), cached))

uploaded <- data.frame(SourceLabel = "uploaded", Value = 2, stringsAsFactors = FALSE)
upload_path <- tempfile(fileext = ".csv")
write.csv(uploaded, upload_path, row.names = FALSE)
upload_info <- data.frame(name = "identifications.csv", datapath = upload_path, stringsAsFactors = FALSE)
selected <- app_env$project_input_table(upload_info, cached)
stopifnot(identical(as.character(selected$SourceLabel), "uploaded"))
stopifnot(identical(as.integer(selected$Value), 2L))

cat("Identification inputs prefer current uploads and fall back to restored project tables.\n")
