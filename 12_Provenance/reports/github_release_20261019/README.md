# Release verification — 19 October 2026

[Saved-data validation](saved_data_validation.json) passed all 97 numerical checks when executed from this curated package. The run used the available local Python environment and the NumPy/Pandas dependency fallback in the prior `Manuscript_1010/Code/vendor` directory through `PYTHONPATH`. A fresh installation from `requirements-demo.txt` has not been tested.

```bash
python Code/reproduce.py --out /tmp/rna-velocity-1019-check --no-plots
```

The eight differences between frozen overall scores and current source-derived values are preserved as documented provenance differences, not silently replaced. No models were retrained and no raw sequencing inputs were reprocessed. Historical figure identifiers remain in the checker output.

The saved-data checker skips its optional original-workspace checksum check because the full workspace checksum manifest is excluded. The separate [package validation](../../release_validation.json) checks SHA256 identity for every copied source file, figure order and local navigation:

```bash
# Run from the repository root.
python tools/verify_release.py
```
