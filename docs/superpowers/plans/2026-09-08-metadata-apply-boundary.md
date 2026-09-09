# Metadata Apply Boundary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Make metadata tab a draft workspace and update every downstream tab only after one successful **Apply metadata changes** action.

**Architecture:** Extract pure metadata-state helpers, retain the current metadata builder as `draft_metadata()`, and introduce one immutable applied snapshot consumed by the existing `built_metadata()` downstream interface. The Apply observer validates and atomically replaces metadata, Table S1 scope, and committed SPQC edits; Make metadata previews remain draft-aware and SPQC cell edits no longer replace the DataTable.

**Tech Stack:** R 4.4, Shiny reactive values, DT editable tables, existing DuckDB/RDS project helpers, base-R tests.

**Spec:** `docs/superpowers/specs/2026-09-08-metadata-apply-boundary-design.md`

## Global Constraints

- Do not change SPQC inference, exclusion matching, joins, abundance processing, imputation, batch correction, or statistical calculations.
- Do not add R package dependencies.
- A failed Apply must leave the preceding applied snapshot unchanged.
- Project restore must begin in a fully applied, clean state.
- Project save must treat applied metadata as authoritative and must not silently persist draft metadata.
- The app directory is not a Git repository; replace each commit step with a recorded verification checkpoint unless it becomes a repository before execution.

---

### Task 1: Pure applied-metadata state helpers

**Files:**
- Create: `metadata_state.R`
- Create: `tests/test_metadata_apply_state.R`
- Modify: `app.R:20-40` (source the helper)

**Interfaces:**
- Produces: `empty_proteomics_metadata_state() -> list`
- Produces: `merge_proteomics_metadata_edits(committed, draft) -> data.frame`
- Produces: `build_proteomics_applied_state(metadata, columns, committed_edits, draft_edits) -> list(metadata, columns, spqc_metadata_edits)`
- Consumes: `apply_proteomics_metadata_cell_edits()` from `project_duckdb.R`

- [ ] **Step 1: Write the failing helper tests**

```r
source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")
source(if (file.exists("metadata_state.R")) "metadata_state.R" else "../metadata_state.R")

metadata <- data.frame(
  Sample = c("SPQC_1", "Sample_1"),
  Condition = c("SPQC", "Control"),
  Batch = c("1", "1"),
  stringsAsFactors = FALSE
)
committed <- data.frame(Sample = "SPQC_1", Column = "Batch", Value = "1", stringsAsFactors = FALSE)
draft <- data.frame(Sample = "SPQC_1", Column = "Batch", Value = "2", stringsAsFactors = FALSE)

state <- build_proteomics_applied_state(metadata, c("Sample", "Condition"), committed, draft)
stopifnot(identical(state$metadata$Batch, c("2", "1")))
stopifnot(identical(state$columns, c("Sample", "Condition")))
stopifnot(nrow(state$spqc_metadata_edits) == 1L, state$spqc_metadata_edits$Value[[1L]] == "2")

before <- state
try(build_proteomics_applied_state(metadata, "MissingColumn", committed, draft), silent = TRUE)
stopifnot(identical(state, before))
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' tests\test_metadata_apply_state.R
```

Expected: FAIL because `metadata_state.R` or `build_proteomics_applied_state()` does not exist.

- [ ] **Step 3: Implement the helpers**

```r
empty_proteomics_metadata_state <- function() {
  list(metadata = NULL, columns = character(0), spqc_metadata_edits = empty_spqc_metadata_edits())
}

merge_proteomics_metadata_edits <- function(committed, draft) {
  if (!nrow(draft)) return(committed)
  draft_keys <- paste(draft$Sample, draft$Column, sep = "\r")
  if (nrow(committed)) {
    committed_keys <- paste(committed$Sample, committed$Column, sep = "\r")
    committed <- committed[!committed_keys %in% draft_keys, , drop = FALSE]
  }
  rbind(committed, draft)
}

build_proteomics_applied_state <- function(metadata, columns, committed_edits, draft_edits) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(metadata) || !"Sample" %in% colnames(metadata)) stop("Metadata must contain at least one Sample row.")
  if (anyDuplicated(as.character(metadata$Sample))) stop("Metadata Sample values must be unique before applying changes.")
  columns <- unique(as.character(columns))
  if (!length(columns) || any(!columns %in% colnames(metadata))) stop("Every selected Table S1 column must exist in metadata.")
  edits <- merge_proteomics_metadata_edits(committed_edits, draft_edits)
  list(
    metadata = apply_proteomics_metadata_cell_edits(metadata, edits),
    columns = columns,
    spqc_metadata_edits = edits
  )
}
```

Source `metadata_state.R` after `project_duckdb.R` so its edit-overlay dependency is available.

- [ ] **Step 4: Run the helper test and verify GREEN**

Expected: `Metadata apply state tests passed.` and exit code 0.

- [ ] **Step 5: Record checkpoint**

Record that `metadata_state.R`, its source line, and `tests/test_metadata_apply_state.R` are complete. Do not attempt a Git commit in the non-repository directory.

---

### Task 2: Split live draft metadata from the applied snapshot

**Files:**
- Modify: `app.R:1930-1945` (state declarations)
- Modify: `app.R:3520-3670` (metadata builder)
- Modify: `tests/test_metadata_scope_server.R`
- Create: `tests/test_metadata_draft_boundary_static.R`

**Interfaces:**
- Consumes: `empty_proteomics_metadata_state()` from Task 1
- Produces: `draft_metadata() -> data.frame`
- Produces: `applied_metadata_state <- reactiveVal(list)`
- Preserves: `built_metadata() -> data.frame` as the downstream applied-metadata accessor

- [ ] **Step 1: Write the failing boundary test**

Add static assertions:

```r
text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
stopifnot(grepl("applied_metadata_state <- reactiveVal(empty_proteomics_metadata_state())", text, fixed = TRUE))
stopifnot(grepl("draft_metadata <- reactive({", text, fixed = TRUE))
stopifnot(grepl("built_metadata <- reactive({\n  state <- applied_metadata_state()", text, fixed = TRUE))
stopifnot(grepl('validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))', text, fixed = TRUE))
```

Extend the live server test, where dependencies are available, to capture `built_metadata()` after an initial Apply, change `input$spqc_assignment_mode`, and assert that `built_metadata()` is unchanged before the next Apply.

- [ ] **Step 2: Run both tests and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' tests\test_metadata_draft_boundary_static.R
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' tests\test_metadata_scope_server.R
```

Expected: the static test fails because the draft/applied split is absent; the live test either fails for immediate propagation or reports its existing dependency skip.

- [ ] **Step 3: Introduce the applied state and rename the builder**

Add:

```r
applied_metadata_state <- reactiveVal(empty_proteomics_metadata_state())
metadata_apply_revision <- reactiveVal(0L)
metadata_apply_message <- reactiveVal("Metadata changes have not been applied.")
```

Rename the existing `built_metadata <- reactive({ ... })` body to `draft_metadata <- reactive({ ... })`. Then define the stable downstream interface:

```r
built_metadata <- reactive({
  metadata_apply_revision()
  state <- applied_metadata_state()
  validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))
  state$metadata
})
```

Do not change downstream consumers yet; preserving the `built_metadata()` name freezes them automatically.

- [ ] **Step 4: Initialize restored project state**

In both DuckDB/RDS restoration (`load_project_into_session`) and legacy ZIP restoration, build the applied state from restored metadata plus restored `metadata_export_columns`. If saved columns are empty, use `metadata_default_columns(colnames(metadata))`. Clear draft SPQC edits and set the apply message to `"All metadata changes applied."`.

- [ ] **Step 5: Run the boundary tests and verify GREEN**

Expected: the static boundary checks pass. The live test passes when its Shiny dependencies are available and otherwise retains its explicit skip message.

- [ ] **Step 6: Record checkpoint**

Record the draft/applied split and restore initialization as a completed unit.

---

### Task 3: Make one Apply action commit the whole metadata tab atomically

**Files:**
- Modify: `app.R:1025-1120` (metadata UI)
- Modify: `app.R:2780-2860` (replacement and Apply observers)
- Modify: `app.R:7200-7260` (remove SPQC-only Apply observer)
- Modify: `tests/test_spqc_edit_stability_static.R`
- Create: `tests/test_metadata_apply_transaction_static.R`

**Interfaces:**
- Consumes: `draft_metadata()` and `build_proteomics_applied_state()`
- Updates: `applied_metadata_state`, `spqc_metadata_edits`, `spqc_metadata_draft_edits`, `metadata_apply_revision`
- Produces UI input: `input$apply_metadata_changes`

- [ ] **Step 1: Write the failing transaction/UI tests**

Require these exact contracts:

```r
stopifnot(grepl('actionButton("apply_metadata_changes", "Apply metadata changes", class = "btn-primary")', text, fixed = TRUE))
stopifnot(!grepl('actionButton("apply_spqc_metadata_edits"', text, fixed = TRUE))
stopifnot(grepl('actionButton("apply_metadata_replacement", "Load modified metadata into draft")', text, fixed = TRUE))
stopifnot(grepl('candidate <- build_proteomics_applied_state(', text, fixed = TRUE))
stopifnot(grepl('applied_metadata_state(candidate)', text, fixed = TRUE))
stopifnot(grepl('metadata_apply_revision(isolate(metadata_apply_revision()) + 1L)', text, fixed = TRUE))
stopifnot(grepl('autosave_active_project("applied metadata changes", include_derived = FALSE)', text, fixed = TRUE))
```

Also assert that `applied_metadata_state(candidate)` occurs after the candidate construction and validation lines, ensuring the old state is not mutated before validation completes.

- [ ] **Step 2: Run the transaction test and verify RED**

Expected: FAIL because `apply_metadata_changes` and the atomic candidate assignment are absent.

- [ ] **Step 3: Consolidate the UI**

Move one prominent button and status near the **Build analysis metadata** heading:

```r
actionButton("apply_metadata_changes", "Apply metadata changes", class = "btn-primary"),
textOutput("metadata_apply_status")
```

Remove the SPQC-only Apply button. Rename the replacement action to **Load modified metadata into draft** and revise its note to say downstream tabs are unchanged until the main Apply action.

- [ ] **Step 4: Make replacement upload draft-only**

Change its observer to validate against `draft_metadata()`, then set `active_metadata_override(replacement)` and clear draft SPQC edits. Do not update `project_db_cache`, invalidate batch/statistics caches, or autosave in this observer.

- [ ] **Step 5: Implement the atomic Apply observer**

```r
observeEvent(input$apply_metadata_changes, {
  tryCatch({
    candidate <- build_proteomics_applied_state(
      metadata = draft_metadata(),
      columns = input$metadata_export_columns,
      committed_edits = spqc_metadata_edits(),
      draft_edits = spqc_metadata_draft_edits()
    )
    cache <- invalidate_proteomics_project_cache(project_db_cache(), "meta_file")
    cache$metadata <- candidate$metadata
    project_db_cache(cache)
    cached_batch_corrected_s3_result(NULL)
    spqc_metadata_edits(candidate$spqc_metadata_edits)
    spqc_metadata_draft_edits(empty_spqc_metadata_edits())
    applied_metadata_state(candidate)
    metadata_apply_revision(isolate(metadata_apply_revision()) + 1L)
    metadata_apply_message("All metadata changes applied.")
    autosave_active_project("applied metadata changes", include_derived = FALSE)
    showNotification("Metadata changes applied to all tabs.", type = "message")
  }, error = function(e) {
    metadata_apply_message(paste0("Metadata was not applied: ", conditionMessage(e)))
    showNotification(conditionMessage(e), type = "error", duration = NULL)
  })
})
```

Do not assign any applied-state reactive before `candidate` has been created successfully.

- [ ] **Step 6: Add applied/dirty status**

Render `metadata_apply_status` from the persistent apply message plus whether draft SPQC edits exist. Metadata controls may mark the message `"Unapplied metadata changes."`; this status change must not mutate `applied_metadata_state` or `metadata_apply_revision`.

- [ ] **Step 7: Run transaction tests and verify GREEN**

Expected: both static tests pass and the old SPQC-only Apply ID is absent.

- [ ] **Step 8: Record checkpoint**

Record the UI and atomic Apply transaction as a completed unit.

---

### Task 4: Keep Make metadata previews live and stabilize SPQC editing

**Files:**
- Modify: `app.R:6970-7050` (metadata notes and exclusion preview)
- Modify: `app.R:7170-7260` (metadata and SPQC DataTables)
- Modify: `tests/test_spqc_edit_stability_static.R`
- Modify: `tests/test_metadata_draft_boundary_static.R`

**Interfaces:**
- Consumes: `draft_metadata()` for Make metadata previews
- Consumes: `spqc_metadata_draft_edits()` as an overlay
- Must not call: `replaceData()` from the SPQC cell-edit observer

- [ ] **Step 1: Strengthen the failing SPQC regression test**

```r
edit_start <- regexpr("observeEvent(input$spqc_metadata_preview_cell_edit", text, fixed = TRUE)[1]
apply_start <- regexpr("observeEvent(input$apply_metadata_changes", text, fixed = TRUE)[1]
edit_block <- substr(text, edit_start, apply_start - 1L)
stopifnot(grepl("spqc_metadata_draft_edits(draft_edits)", edit_block, fixed = TRUE))
stopifnot(!grepl("replaceData(", edit_block, fixed = TRUE))
stopifnot(!grepl("spqc_metadata_edits(", edit_block, fixed = TRUE))
```

Require `draft_metadata()` in `metadata_build_note`, `sample_exclusion_preview`, `metadata_preview`, and `spqc_metadata_rows`.

- [ ] **Step 2: Run the test and verify RED**

Expected: FAIL because the edit block still calls `replaceData()` and previews still use `built_metadata()`.

- [ ] **Step 3: Remove proxy replacement from SPQC edits**

Retain the bounds checks and Sample/Column/Value draft record, but remove `spqc_metadata_proxy`, `updated_table`, and `replaceData()`. DT already leaves the browser-edited value visible; recording the edit must not redraw the widget.

- [ ] **Step 4: Point tab-local previews to the draft**

Use `draft_metadata()` plus the draft SPQC overlay for:

- metadata build note;
- exclusion note and preview;
- Table S1 metadata preview on the Make metadata tab;
- detected SPQC rows and editable SPQC table.

Keep authoritative Table S1/project downloads on applied `built_metadata()` unless their label explicitly says draft.

- [ ] **Step 5: Verify the focused tests GREEN**

Expected: SPQC stability and draft-boundary tests pass.

- [ ] **Step 6: Manual local UI check**

Run the app, edit two SPQC cells on a non-first DataTable page, and verify:

1. the table remains visible and on the same page;
2. downstream group choices do not change before Apply;
3. both edits remain visible;
4. one Apply updates downstream group choices and clears the unapplied status.

- [ ] **Step 7: Record checkpoint**

Record the exact manual steps and outcome because the static test cannot prove browser focus/scroll behavior.

---

### Task 5: Apply Table S1 scope downstream and preserve project semantics

**Files:**
- Modify: `app.R:2670-2710` (save payload)
- Modify: `app.R:3740-3980` (metadata scopes/selectors)
- Modify: `tests/test_metadata_scope.R`
- Modify: `tests/test_project_restored_tables.R`
- Modify: `tests/test_project_settings_static.R`

**Interfaces:**
- Consumes: `applied_metadata_state()$columns`
- Preserves: `downstream_metadata_columns() -> character`
- Persists: applied metadata, applied column scope, committed SPQC edits

- [ ] **Step 1: Write failing scope and persistence assertions**

```r
stopifnot(grepl("state <- applied_metadata_state()", text, fixed = TRUE))
stopifnot(grepl("downstream_metadata_columns <- reactive({\n  state <- applied_metadata_state()", text, fixed = TRUE))
stopifnot(!grepl("table_s1_metadata_choices(available, input$metadata_export_columns)", text, fixed = TRUE))
```

Extend project tests so restored metadata and restored `metadata_export_columns` initialize the same applied state, while changing the selector alone does not alter `downstream_metadata_columns()`.

- [ ] **Step 2: Run tests and verify RED**

Expected: FAIL because downstream scope still reads `input$metadata_export_columns` directly.

- [ ] **Step 3: Freeze downstream Table S1 scope**

```r
downstream_metadata_columns <- reactive({
  state <- applied_metadata_state()
  validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))
  table_s1_metadata_choices(colnames(state$metadata), state$columns)
})
```

The Make metadata selector continues to display the draft selection, but no downstream selector may read it directly.

- [ ] **Step 4: Persist only applied state**

Build save payload metadata and `settings$metadata_export_columns` from `applied_metadata_state()`. Keep `spqc_metadata_draft_edits` out of the authoritative project payload. If unapplied changes exist during manual save, show a warning that the saved project's active analysis metadata is the last applied version.

- [ ] **Step 5: Verify project restore**

When loading a project, set the applied snapshot before incrementing `project_restore_token`, then update the draft controls from restored settings. Confirm no stale draft edit survives the restore.

- [ ] **Step 6: Run focused tests and verify GREEN**

Run metadata-scope, project-restoration, project-settings, and SPQC tests. Expected: all pass.

- [ ] **Step 7: Run complete verification**

```powershell
$r='C:\Program Files\R\R-4.4.2\bin\Rscript.exe'
& $r -e "invisible(parse(file='app.R')); cat('app.R parse OK\n')"
$failed=@()
Get-ChildItem -LiteralPath tests -Filter '*.R' | Sort-Object Name | ForEach-Object {
  & $r $_.FullName
  if ($LASTEXITCODE -ne 0) { $failed += $_.Name }
}
if ($failed.Count) { throw ('Failed tests: ' + ($failed -join ', ')) }
```

Expected: parse succeeds and every test script exits 0. Existing explicit skips for unavailable live Shiny or DuckDB dependencies are acceptable only when their skip reason is printed.

- [ ] **Step 8: Final manual acceptance check**

With a restored DuckDB project:

1. note a downstream metadata selector's choices;
2. change exclusions, SPQC assignment mode, Table S1 columns, and two SPQC cells;
3. verify Make metadata previews change but the downstream choices do not;
4. click Apply once;
5. verify all changes appear downstream together;
6. save and reopen the project;
7. verify the applied metadata and scope restore with no dirty status.

- [ ] **Step 9: Record final checkpoint**

Record parse result, number of passing test scripts, skip reasons, and manual acceptance result. Do not claim completion unless all non-skipped checks pass.
