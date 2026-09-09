app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

metadata <- data.frame(
  Sample = sprintf("run_%02d.raw", seq_len(49)),
  SampleName = sprintf("subject_%02d.raw", seq_len(49)),
  stringsAsFactors = FALSE
)

flags <- app_env$sample_exclusion_flags(metadata, "run_01")
stopifnot(length(flags) == nrow(metadata))
stopifnot(identical(which(flags), 1L))

cat("Sample exclusion flags retain one value per metadata row.\n")
