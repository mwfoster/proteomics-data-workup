app_expressions <- parse(file = "app.R")

function_assignment <- function(name) {
  matches <- vapply(app_expressions, function(expression) {
    is.call(expression) &&
      identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name(name))
  }, logical(1))

  stopifnot(sum(matches) == 1L)
  app_expressions[[which(matches)]]
}

test_environment <- new.env(parent = globalenv())
eval(function_assignment("calculate_protein_comparison"), envir = test_environment)

report <- data.frame(
  Protein = "P1",
  numerator_1 = 2,
  numerator_2 = 8,
  denominator_1 = 1,
  denominator_2 = 2,
  check.names = FALSE
)

result <- test_environment$calculate_protein_comparison(
  report,
  numerator_cols = c("numerator_1", "numerator_2"),
  denominator_cols = c("denominator_1", "denominator_2"),
  paired = TRUE
)

# The paired log2 ratios are 1 and 2, so their arithmetic mean is 1.5.
# This fixture distinguishes the requested method from log2(mean(2, 4)).
stopifnot(isTRUE(all.equal(result$data$log2_fc[[1L]], 1.5, tolerance = 1e-12)))

eval(function_assignment("replicate_pair_plan"), envir = test_environment)
md <- data.frame(Condition = c("A", "A", "B", "B"),
                 Replicate = c(1, 2, 1, 2), Subject = c("X", "Y", "Y", "X"))
labels <- c("AX", "AY", "BY", "BX")
default_plan <- test_environment$replicate_pair_plan(md, labels, "A", "B")
subject_plan <- test_environment$replicate_pair_plan(md, labels, "A", "B", pair_col = "Subject")
stopifnot(default_plan$balanced, subject_plan$balanced)
stopifnot(identical(default_plan$denominator_labels, c("BY", "BX")))
stopifnot(identical(subject_plan$denominator_labels, c("BX", "BY")))
stopifnot(!test_environment$replicate_pair_plan(md, labels, "A", "B", pair_col = "Missing")$balanced)
md$Subject[4] <- "Y"
stopifnot(!test_environment$replicate_pair_plan(md, labels, "A", "B", pair_col = "Subject")$balanced)

message("Paired fold-change and metadata pairing tests passed.")
