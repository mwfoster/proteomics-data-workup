source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

columns <- c("Protein", "ProteinName", "ProteinDescription", "Status", "Log2FoldChange", "PValue", "BH_FDR")
stopifnot(identical(volcano_hits_display_columns(columns, TRUE), columns))
stopifnot(identical(volcano_hits_display_columns(columns, FALSE), setdiff(columns, "BH_FDR")))

cat("Volcano FDR visibility tests passed.\n")
