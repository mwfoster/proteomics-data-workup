app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

input <- matrix(c(1, NaN, Inf, -Inf, 5, NA_real_), nrow = 2)
output <- app_env$normalize_pca_missing_values(input)

stopifnot(is.matrix(output))
stopifnot(identical(dim(output), dim(input)))
stopifnot(all(is.na(output[c(2, 3, 4, 6)])))
stopifnot(identical(output[c(1, 5)], c(1, 5)))

cat("PCA input converts all non-finite values to NA.\n")
