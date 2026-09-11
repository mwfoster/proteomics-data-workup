license <- paste(readLines("LICENSE", warn = FALSE), collapse = "\n")
readme <- paste(readLines("README.md", warn = FALSE), collapse = "\n")

stopifnot(grepl("MIT License", license, fixed = TRUE))
stopifnot(grepl("Copyright (c) 2026 Matthew W. Foster", license, fixed = TRUE))
stopifnot(grepl("Permission is hereby granted, free of charge", license, fixed = TRUE))
stopifnot(grepl("## License", readme, fixed = TRUE))
stopifnot(grepl("MIT License", readme, fixed = TRUE))

cat("Release license tests passed.\n")
