app_path <- if (file.exists("app.R")) "app.R" else "../app.R"
app_env <- new.env(parent = globalenv())
sys.source(app_path, envir = app_env)

html <- htmltools::renderTags(app_env$ui)$html
for (input_id in c("feature_protein_label", "script_box_protein_label", "correlation_protein_label")) {
  select_tag <- regmatches(html, regexpr(paste0('(?s)<select id="', input_id, '"[^>]*>.*?</select>'), html, perl = TRUE))
  stopifnot(length(select_tag) == 1L, nzchar(select_tag))
  stopifnot(grepl('value="protein_name"[^>]*selected', select_tag))
  stopifnot(grepl('value="accession"', select_tag, fixed = TRUE))
}

cat("Protein plot label controls passed.\n")
