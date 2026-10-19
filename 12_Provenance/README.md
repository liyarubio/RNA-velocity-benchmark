# Provenance and validation

- [Figure order and historical numbering](figure_index.csv)
- [Figure-file bindings](file_inventory.csv)
- [Workbook sheets](workbook_sheets.csv)
- [Old-to-new paths and SHA256](path_migration.json)
- [Migration table](path_migration.csv)
- [Current structure verification](restructure_validation.json)
- [Numerical and plot regression results](restructure_regression.json)
- [Executed saved-data checks](saved_data_validation.json)
- [English release changes and hashes](english_translation.json)
- [English release validation](english_validation.json)
- [English release scope](ENGLISH_RELEASE.md)
- [Pre-restructure release manifests](baseline_1019/)
- [Inherited reports](reports/)

Inherited records describe their original release and retain original paths; they are not certificates for the current tree. The migration manifest covers every file tracked at the baseline commit. The initial directory reorganization preserved source-data and artwork bytes, with path adaptations limited to portable tools. The subsequent English release translates documentation, code text and workbook annotations or style names; its separate ledger records the previous and current hashes. Numerical cells, formulas and formal artwork are unchanged. Historical navigation and scripts are retained in the archive with English text.

No new biological analyses, model fits or enrichment requests were performed. Incomplete age provenance, unavailable values and eight documented frozen-versus-current score differences remain explicit.
