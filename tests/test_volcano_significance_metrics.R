source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

p_value <- c(0.01, 0.2)
fdr <- c(0.02, 0.25)
without_fdr <- volcano_significance_metric_values(p_value, fdr, FALSE)
stopifnot(identical(names(without_fdr), "p-value"), identical(without_fdr[[1]], p_value))

with_fdr <- volcano_significance_metric_values(p_value, fdr, TRUE)
stopifnot(identical(names(with_fdr), c("BH FDR", "p-value")))

cat("Volcano significance metric tests passed.\n")
