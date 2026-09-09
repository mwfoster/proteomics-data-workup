# Metadata Draft and Apply Boundary

## Purpose

The **Make metadata** tab must act as an editing workspace. Changes made there must not alter metadata used by Protein tables, PCA, CV, statistics, batch correction, feature plots, boxplots, correlation, or exports derived from those analyses until the user clicks **Apply metadata changes**.

This design also fixes the editable SPQC table disappearing after a cell edit.

## User-visible behavior

- The Make metadata tab displays a live draft assembled from its current source files, replacement metadata, exclusion rules, SPQC assignment controls, Table S1 column selection, and direct SPQC cell edits.
- Editing an SPQC cell updates the visible cell and records a draft edit. It does not rebuild or replace the DataTable.
- A status message identifies whether the metadata draft has unapplied changes.
- Clicking **Apply metadata changes** validates and commits the entire metadata-tab draft in one operation.
- Subsequent tabs continue using the last successfully applied metadata snapshot until that operation succeeds.
- If applying fails validation, the existing applied snapshot remains unchanged and the error is shown on the Make metadata tab.
- Loading an existing saved project establishes its restored metadata and Table S1 column scope as the applied snapshot automatically.
- Starting or clearing a project removes both draft and applied session state.
- The existing modified-metadata upload action is labeled **Load modified metadata into draft**. It updates only the draft until the main Apply button is clicked.

## State model

The app will maintain two distinct metadata layers:

1. `draft_metadata`: reactive metadata generated from controls on the Make metadata tab. It may change frequently and is used only by that tab's previews and downloads explicitly labeled as draft.
2. `applied_metadata`: a `reactiveVal` containing the last validated immutable snapshot. All downstream tabs and derived analyses consume this layer.

The app will also maintain:

- draft SPQC cell edits, separate from committed SPQC edits;
- an applied Table S1 metadata-column scope, separate from the current draft selector value;
- an apply revision or equivalent reactive signal so downstream consumers invalidate exactly once after a successful commit;
- a draft-dirty indicator based on changes since the last successful apply.

## Data flow

### Draft editing

Source uploads, replacement metadata, exclusions, SPQC assignment controls, and SPQC cell edits rebuild or amend `draft_metadata`. Make metadata previews use this draft. No downstream reactive reads draft state.

The editable SPQC DataTable performs its normal browser-side cell update. Its edit observer records the changed Sample/Column/Value tuple only. It does not call `replaceData()` for each edit, which avoids destroying or blanking the active widget.

### Apply

When **Apply metadata changes** is clicked, the app:

1. builds the current draft;
2. overlays all draft SPQC cell edits;
3. validates sample identifiers, row uniqueness, selected Table S1 columns, and required metadata structure;
4. stores the resulting data frame in `applied_metadata`;
5. stores the current Table S1 column selection as the applied downstream scope;
6. moves draft SPQC edits into committed project edits;
7. clears the dirty state;
8. invalidates downstream consumers once;
9. requests autosave for an active project.

The operation is atomic from the downstream app's perspective: validation failure changes none of the applied state.

### Project restore

Opening DuckDB, RDS, or legacy ZIP project data populates both draft and applied metadata from the restored project. Restored Table S1 column settings populate both draft and applied column scope. The project begins in a clean, fully applied state.

### Project persistence

Project save and autosave persist only applied metadata and applied Table S1 scope as authoritative analysis state. Unapplied draft changes are not silently written as active analysis metadata. The UI warns that they remain unapplied.

## Component changes

### Metadata builder

The current metadata-building body will become the draft builder. A separate applied-metadata accessor will be the only accessor used outside the Make metadata tab.

### Make metadata outputs

Metadata build notes, exclusion previews, metadata preview, and SPQC preview use draft metadata. The main Table S1 download will use applied metadata; if a draft download is retained, it must be labeled explicitly.

### Downstream consumers

Protein table header mapping, metadata selectors, batch correction, PCA, CV, statistics, feature plots, boxplots, correlation, and analysis exports use applied metadata and the applied Table S1 column scope. They must not depend directly on current metadata-tab inputs.

### Apply controls

There will be one prominent **Apply metadata changes** control for the whole tab. The status beside it reports either **All metadata changes applied** or **Unapplied metadata changes**.

## Error handling

- Apply errors are caught and displayed with a Shiny notification and persistent tab status.
- Failed validation leaves the prior applied metadata and downstream analysis state untouched.
- Applying without usable metadata displays a clear validation message.
- Applying with no changes is allowed and reports that metadata is already current without forcing unnecessary downstream recomputation.
- A project restore clears stale draft edits before establishing the restored snapshot.

## Testing

Tests will verify:

- SPQC cell edits do not call `replaceData()` and do not mutate applied metadata.
- Draft controls can change while applied metadata remains unchanged.
- A successful Apply commits all metadata-tab changes together and increments downstream state once.
- A failed Apply preserves the preceding applied snapshot.
- Applied Table S1 columns, not draft selections, control downstream metadata choices.
- Project restore initializes draft and applied state consistently.
- Project save uses applied metadata and excludes unapplied drafts from authoritative project state.
- Existing metadata, DuckDB, statistics, PCA, CV, and plot-selection regression tests continue to pass.

## Out of scope

- Changing the scientific logic used to infer SPQC samples, exclusions, or metadata joins.
- Changing protein abundance processing, imputation, batch correction, or statistical calculations.
- Implementing browser-local folder access for a future hosted VM version.
