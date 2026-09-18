app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

overview <- data.frame(
  Condition = c("Control", "Treatment"),
  Replicate = c("1", "1"),
  stringsAsFactors = FALSE
)
metadata <- data.frame(
  Condition = c("Control", "Treatment"),
  Replicate = c("1", "1"),
  Subject = c("S01", "S02"),
  stringsAsFactors = FALSE
)

labels <- app_env$identification_metadata_labels(
  overview,
  metadata,
  selected_columns = c("Subject", "Condition"),
  fallback = c("Control_1", "Treatment_1")
)
stopifnot(identical(labels, c("S01_Control", "S02_Treatment")))

fallback <- app_env$identification_metadata_labels(
  overview,
  metadata,
  selected_columns = character(0),
  fallback = c("Control_1", "Treatment_1")
)
stopifnot(identical(fallback, c("Control_1", "Treatment_1")))

cat("Identification plots can build labels from selected metadata columns.\n")
