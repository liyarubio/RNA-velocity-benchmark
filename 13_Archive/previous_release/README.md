# RNA velocity benchmark — 1019 code and data release

This package organizes the reviewed 2026-10-19 benchmark materials in manuscript order: **7 main figures and 12 Extended Data figures**. It contains saved numerical results, source workbooks, available historical implementations, current figure PDFs/PNGs and provenance records.

- [Browse code, data and figures in manuscript order](Figure_materials/README.md)
- [Per-file inventory and SHA256](Figure_materials/file_inventory.csv)
- [Workbook sheet index](Figure_materials/workbook_sheets.csv)
- [Main and supplementary tables](Figure_materials/03_Tables/README.md)
- [Package manifest](release_manifest.json)
- [Package validation](release_validation.json)
- [Executed saved-data checks and environment limits](Reports/github_release_20261019/README.md)

## Saved-data verification

Use a compatible Python environment with the dependencies in `requirements-demo.txt`:

```bash
python -m pip install -r requirements-demo.txt
python Code/reproduce.py --out /tmp/rna-velocity-1019-check --no-plots
python Code/plot_saved_panels.py --out /tmp/rna-velocity-1019-panels
```

The first command audits saved values. The second produces diagnostic numerical panels, not every publication layout. Neither retrains velocity models or validates new biological claims. Historical ED IDs in these scripts map to current numbering through the figure index.

## Reproduction boundaries

The per-figure READMEs distinguish saved-data checks, historical model/plotting scripts and editorial layout code. Historical scripts can depend on external H5AD files, paths and original environments. Layout scripts can depend on the full local manuscript workspace and earlier releases; they are retained for inspection, not advertised as standalone generators.

Missing values remain missing. Reference comparisons retain their original preprocessing limitations. Frozen overall scores, incomplete age-example provenance and other evidence gaps remain documented. See [inherited evidence gaps](Reports/inherited_1011/evidence_gaps.csv).

The full draft manuscript, editing backups, literature-download caches and duplicate SVG/combined exports are not part of this code/data package. The repository's previous analysis directories are retained separately.

## Attribution and licensing

No new blanket license is assigned by this update. Existing third-party terms and author rights remain applicable; see [LICENSE_STATUS.txt](LICENSE_STATUS.txt). The upstream implementations retain their provenance and historical configuration.
