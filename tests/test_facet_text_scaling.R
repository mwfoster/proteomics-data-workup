source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

stopifnot(identical(facet_text_scale(1), 1))
stopifnot(identical(facet_text_scale(2), 0.9))
stopifnot(identical(facet_text_scale(3), 0.8))
stopifnot(identical(facet_text_scale(4), 0.7))
stopifnot(identical(facet_text_scale(8), 0.7))
stopifnot(identical(facet_text_scale(NA), 1))

cat("Facet text scaling tests passed.\n")
