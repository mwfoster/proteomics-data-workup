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

overview_by_run <- data.frame(
  FileName = c("run-a.htrms", "run-b.htrms"),
  Condition = c("FA", "O3"),
  Replicate = c("1", "1"),
  stringsAsFactors = FALSE
)
metadata_by_run <- data.frame(
  Sample = c("run-a.htrms", "run-b.htrms"),
  SampleName = c("Subject01 FA", "Subject01 O3"),
  Condition = c("Apical FA", "Apical O3"),
  Replicate = c("8", "8"),
  stringsAsFactors = FALSE
)
run_labels <- app_env$identification_metadata_labels(
  overview_by_run,
  metadata_by_run,
  selected_columns = "SampleName",
  fallback = c("FA_1", "O3_1")
)
stopifnot(identical(run_labels, c("Subject01 FA", "Subject01 O3")))

cat("Identification plots can build labels from selected metadata columns.\n")
