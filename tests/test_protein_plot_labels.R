source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

info <- data.frame(
  Protein = c("gene_key_1", "gene_key_2", "gene_key_3"),
  PG.ProteinNames = c("CXCL8_HUMAN", "", "CXCL8_HUMAN"),
  PG.ProteinGroups = c("P10145", "Q99999", "P22222"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
features <- c("gene_key_1", "gene_key_2", "gene_key_3", "unmatched")

name_labels <- protein_plot_labels(info, features, "protein_name")
stopifnot(identical(unname(name_labels), c("CXCL8_HUMAN", "Q99999", "CXCL8_HUMAN.1", "unmatched")))
stopifnot(identical(names(name_labels), features))

accession_labels <- protein_plot_labels(info, features, "accession")
stopifnot(identical(unname(accession_labels), c("P10145", "Q99999", "P22222", "unmatched")))

single_plot <- protein_plot_display(info, "gene_key_1", "protein_name", "", "Selected proteins")
stopifnot(identical(single_plot$title, "CXCL8_HUMAN"))
stopifnot(identical(unname(single_plot$labels), "CXCL8_HUMAN"))

multiple_plot <- protein_plot_display(info, c("gene_key_1", "gene_key_3"), "accession", "", "Selected proteins")
stopifnot(identical(multiple_plot$title, "Selected proteins"))
stopifnot(identical(unname(multiple_plot$labels), c("P10145", "P22222")))

custom_plot <- protein_plot_display(info, "gene_key_1", "accession", "My title", "Selected proteins")
stopifnot(identical(custom_plot$title, "My title"))

cat("Protein plot label tests passed.\n")
