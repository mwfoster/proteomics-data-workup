app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

expression <- matrix(
  seq_len(12),
  nrow = 4,
  dimnames = list(c("S1", "S2", "S3", "S4"), c("P1", "P2", "P3"))
)
metadata <- data.frame(
  Sample = c("S4", "S2", "S1", "S3"),
  Condition = c("Control", "Case", "Case", "Control"),
  Batch = c("B2", "B1", "B1", "B2"),
  stringsAsFactors = FALSE
)

subset <- app_env$subset_pca_samples(expression, metadata, "Condition", "Case")
stopifnot(identical(rownames(subset$expression), c("S1", "S2")))
stopifnot(identical(as.character(subset$metadata$Sample), c("S1", "S2")))
stopifnot(identical(as.character(subset$metadata$Condition), c("Case", "Case")))
stopifnot(identical(subset$excluded_samples, 2L))

all_samples <- app_env$subset_pca_samples(expression, metadata, "", character(0))
stopifnot(identical(rownames(all_samples$expression), rownames(expression)))
stopifnot(identical(as.character(all_samples$metadata$Sample), rownames(expression)))
stopifnot(identical(all_samples$excluded_samples, 0L))

cat("PCA metadata subset keeps only selected levels and preserves matrix sample order.\n")
