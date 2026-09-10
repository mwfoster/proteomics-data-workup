source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

disabled <- proteomics_download_copy_destination(FALSE, "C:/missing")
stopifnot(!disabled$copy, identical(disabled$message, ""))

missing <- proteomics_download_copy_destination(TRUE, file.path(tempdir(), "folder-that-does-not-exist"))
stopifnot(!missing$copy, grepl("does not exist", missing$message, fixed = TRUE))

available <- proteomics_download_copy_destination(TRUE, tempdir())
stopifnot(available$copy, identical(available$path, normalizePath(tempdir(), winslash = "/", mustWork = TRUE)))

cat("Download-copy destination tests passed.\n")
