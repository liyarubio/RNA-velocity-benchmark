> Historical 19 October release note. Current release and validation: [RELEASE_1021.md](RELEASE_1021.md). The original English validation is preserved in [baseline_1019](baseline_1019/public_1019_english_validation.json).

# English release

The current release uses English for navigation, archived documentation, figure descriptions, code comments and messages, saved notebook messages, and workbook captions and style names. The original bilingual manuscript workspace is outside this repository and remains unchanged. Git history is retained.

## Data and provenance

The [translation ledger](english_translation.json) records hashes before and after language changes, against commit `032a1db845c893a22fd08cb480206906d61622bc`. The path-migration manifest retains its original hashes and marks subsequently translated files explicitly. Current figure-file bindings use current hashes. Inherited manifests still describe their original snapshots.

Workbook numerical values, formulas, sheet names and scientific labels are preserved. A Chinese-only caption is translated; captions that already contain a complete English paragraph retain that paragraph. Chinese font and style names use their English aliases. Original source paths containing non-English characters in table provenance are retained as percent-encoded UTF-8 file URIs; they are historical identifiers, not repository download links. Formal PDFs and preview images are unchanged.

## Legacy code scope

Scientific computation is unchanged. Code comments, log messages and editorial strings are translated. The legacy TFvelo runner under Figure 3 reads its external source directory from `TFVELO_SOURCE`, defaulting to `./TFvelo-main/`; it still requires its historical dependencies and setup. The legacy veloAE metrics output directory is now named `veloAE_metrics_global`.

Historical layout scripts contain English replacement strings. They are provenance utilities, not drop-in editors for the original bilingual manuscript; the original versions remain in Git history. Saved notebook error messages are translated without rerunning notebook cells. No model was retrained and no scientific conclusion was revised.

## Checks

[Validation results](english_validation.json) include a text scan, workbook XML invariants, scientific Python syntax comparisons, saved-data numerical checks and diagnostic-plot comparisons. The reusable checker is `python 09_Core/validation/check_english.py --pdf`. It checks decoded JSON and workbook XML as well as plain text and extracted PDF text. Raster images and vector outlines are not OCR-scanned.
