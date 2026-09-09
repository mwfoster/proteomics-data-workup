app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

metadata <- data.frame(
  Sample = c("sample_1", "sample_2"),
  RunOrder = c(2, 1),
  Condition = c("Control", "Case"),
  stringsAsFactors = FALSE
)

view <- app_env$prepare_metadata_editor_view(metadata, "Condition")
stopifnot(identical(colnames(view$display), "Condition"))
stopifnot(identical(view$sample_keys, c("sample_2", "sample_1")))
stopifnot(identical(as.character(view$display$Condition), c("Case", "Control")))

cat("Metadata editor retains hidden sample keys in displayed row order.\n")
