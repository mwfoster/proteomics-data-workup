compose <- paste(readLines("compose.yaml", warn = FALSE), collapse = "\n")

stopifnot(!grepl("caddy:", compose, fixed = TRUE))
stopifnot(grepl('"127.0.0.1:6875:3838"', compose, fixed = TRUE))
stopifnot(!grepl('"0.0.0.0:3838:3838"', compose, fixed = TRUE))
stopifnot(grepl("Rscript", compose, fixed = TRUE))
stopifnot(!grepl('"curl"', compose, fixed = TRUE))

cat("VM nginx deployment configuration tests passed.\n")
