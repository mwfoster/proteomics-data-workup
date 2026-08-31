# Proteomics DuckDB Project Parity Design

## Goal

Bring the proteomics Shiny app's reusable project-state behavior up to date with the phosphoproteomics app while preserving the proteomics-specific processing and visualization workflow. Replace the existing SQLite project catalog with a DuckDB active-project backend and an RDS fallback. Do not add covariate analysis, confounding matrices, PVCA, or other variance-driver analysis.

## Scope

The migration covers:

- DuckDB project save, open, download, and autosave.
- RDS project files using the same save/load interface.
- ZIP project snapshots containing source files, settings, and an embedded DuckDB cache when DuckDB is available.
- Complete restoration of static and dynamically populated controls.
- Persistence of editable SPQC metadata, exclusions, active workflow tab, and analysis settings.
- Selective cache invalidation when metadata or protein source files change.
- Safe metadata overlay after opening a project database.
- User-facing save, restore, autosave, and invalidation status.
- Regression tests for persistence and restoration behavior.

The migration explicitly excludes covariate-adjusted analysis, confounding diagnostics, PVCA, and variance-driver features.

## Storage Architecture

DuckDB is the primary backend. The project file is selected by extension: `.rds` uses the RDS implementation, while `.duckdb` and `.db` use DuckDB. Both implementations expose the same save/load contract so the Shiny server does not contain backend-specific branches beyond backend selection.

The current SQLite multi-project catalog and its external project-files directory are removed from the user workflow. Existing SQLite files are not deleted automatically. The app will explain that legacy SQLite projects must be reopened through the previous app version or migrated separately; silently interpreting them as DuckDB is prohibited.

The project backend stores these logical objects when available:

- Project manifest and format version.
- Complete project settings.
- Combined metadata and editable SPQC metadata edits.
- Source-file manifest.
- Protein tables required by the active workflow.
- Processed/normalized protein matrix and sample map.
- Batch-corrected matrix and sample map.
- Statistics tables.
- Other metadata-independent cached summaries already produced by the proteomics workflow.

DuckDB saves are atomic: the app writes a temporary database in the target directory, closes it successfully, and only then replaces the active project file. A failed save leaves the previous project file intact. RDS saves use the same temporary-file replacement rule.

## User Interface

The existing SQLite panel is replaced with active-project controls modeled on the phosphoproteomics app:

- Open `.duckdb`, `.db`, `.rds`, or project ZIP.
- Set or create an active project path.
- Save the active project immediately.
- Enable or disable debounced autosave.
- Download a DuckDB snapshot.
- Download a ZIP project bundle.
- Display active path, backend, pending autosave reason, last successful save, and last error.

Existing proteomics file uploads and analysis tabs remain. Restored source files remain active internally even though browser security prevents Shiny from repopulating file-input controls. Uploading a new file overrides only that source slot.

The active workflow tab is saved and restored. Restoration occurs after each control's choices are populated so selections are not discarded by later `updateSelect*()` calls.

## Sticky State

All user-controlled analysis choices are project settings unless they are transient action buttons or download triggers. This includes:

- Metadata construction and export columns.
- SPQC assignment settings, batch overrides, and editable metadata cells.
- Sample exclusions and reasons.
- Protein header-label and quantity-column ordering.
- CV groups and plot settings.
- Statistical tables, ordered comparisons, and paired comparisons.
- Batch-correction inputs and options.
- PCA source, groups, palette, labels, shapes, dimensions, and other display controls.
- Volcano, feature, boxplot, identification, GSEA, and correlation controls.
- Active workflow tab.

Dynamic selectors use a two-stage restore. Saved values are remembered immediately, then intersected with the choices available after metadata or processed data loads. A later choice refresh must preserve the remembered values rather than resetting to defaults.

## Metadata Replacement

After a DuckDB/RDS project is restored, the combined cached metadata is the fallback basis because original metadata source files may not all be present. Uploading one replacement metadata source overlays it onto the cached metadata by stable identifiers:

- Condition/run tables match primarily by `Sample` or normalized run label.
- Sample-details data match by `SampleDetailsID` or the supported ID aliases.
- Values supplied by the replacement file win for matched rows.
- Cached columns not present in the replacement remain intact.
- Unmatched or ambiguous identifiers generate a visible warning and audit summary rather than silently rewriting unrelated rows.

Metadata replacement does not empty the complete project cache. It retains settings and metadata-independent processed protein measurements, while invalidating metadata-dependent sample mappings, batch correction, statistics, and grouped visualization state that must be recalculated.

## Source Replacement and Cache Invalidation

Invalidation is implemented through a pure, tested function mapping the changed input to affected cached objects.

- Metadata-source replacement retains settings and metadata-independent processed measurements; it clears sample mappings, batch-corrected results, statistics, and grouped downstream caches.
- Protein-abundance source replacement retains metadata and settings; it clears the corresponding processed table and every result derived from that table.
- A replacement never clears unrelated source slots or the active project path.

The app reports what was retained and which steps must be rerun.

## Autosave

Autosave is optional and enabled by default after an active project path is established. Major state-changing actions queue an autosave reason. A debounce interval consolidates rapid reactive changes into one write after the app becomes idle.

Autosave errors are caught and displayed without terminating the Shiny session. An unsuccessful autosave does not update the last-successful-save timestamp and does not replace the existing project file.

## Error Handling

Open, restore, manual save, autosave, and project-download handlers use consistent `tryCatch` handling and Shiny notifications. Messages identify the failed operation and preserve the underlying error text without exposing internal stack traces in the UI.

Project restore validates the format version and required tables. Optional derived tables may be absent. Corrupt or incompatible files fail cleanly without clearing the currently active session.

## Migration and Compatibility

The existing SQLite UI and runtime dependency are removed. The legacy `.sqlite` file and external files directory are left untouched on disk. README and installation instructions are updated for `DBI`, `duckdb`, and `jsonlite`; RDS remains available without DuckDB.

Current ZIP projects remain supported where their manifest and source-file layout can be read. New ZIP projects add the embedded DuckDB cache and retain the source files required for reproducibility.

## Testing

Tests must cover:

1. DuckDB round-trip of settings, metadata, processed data, mappings, batch-corrected data, statistics, and SPQC edits.
2. Equivalent RDS round-trip.
3. Atomic-save behavior that preserves the previous project after a simulated write failure.
4. Metadata replacement preserving the overall project while invalidating only dependent objects.
5. Run-order and sample-details overlays updating supplied values while retaining unrelated cached columns.
6. Protein-source replacement invalidating its dependent outputs while retaining metadata and settings.
7. Restoration of static controls and dynamically populated selectors across tabs.
8. Restoration of SPQC edits, sample exclusions, and active workflow tab.
9. Debounced autosave state and error reporting.
10. Static parse checks and existing proteomics calculation/visualization regression tests.

The full suite must pass before the SQLite UI is removed and again after documentation/package updates. Tests that require DuckDB may skip only when the package is genuinely unavailable; cache and overlay logic remains testable without DuckDB.

## Acceptance Criteria

- A user can open a saved proteomics DuckDB or RDS project and recover data, settings, edits, selections, and the prior tab.
- Changing one metadata file does not reset the project.
- The status panel clearly identifies invalidated analyses.
- Autosave cannot corrupt the previous active project file.
- No SQLite project controls remain in the current proteomics UI.
- Proteomics processing and visualizations behave as before unless their inputs were intentionally invalidated.
- Covariate analysis and variance-driver functionality are not added.
