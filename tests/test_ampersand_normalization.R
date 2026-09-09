source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

values <- c("SEL & BRT", "SEL &amp; BRT", "SEL &#38; BRT", "SEL &#x26; BRT")
stopifnot(identical(normalize_proteomics_text(values), rep("SEL & BRT", 4)))

metadata <- data.frame(
  Sample = c("S1", "S2"),
  Condition = c("SEL &amp; BRT", "Veh"),
  stringsAsFactors = FALSE
)
normalized <- normalize_proteomics_metadata(metadata)
stopifnot(identical(normalized$Condition, c("SEL & BRT", "Veh")))

comparison <- make_stats_comparison_id("Condition", "SEL &amp; BRT", "Veh")
parsed <- parse_stats_comparison_id(comparison)
stopifnot(identical(parsed$numerator, "SEL & BRT"))
stopifnot(identical(stats_comparison_label(comparison), "Condition: SEL & BRT vs Veh"))

cat("Encoded and literal ampersands normalize to the same metadata comparison value.\n")
