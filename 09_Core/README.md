# Core workflows

These are shared workflows, separate from per-figure files. Historical code is preserved for inspection; external H5AD inputs, environments and path configuration can still be required. Moving code does not certify a complete model rerun.

- [Simulation: dyngen](simulation/dyngen/): original generation, model runs and evaluation.
- [Simulation: cellDancer](simulation/cellDancer/): original generation, noise comparisons and model runs.
- [Model execution](model_execution/): shared velocity-method implementations.
- [Gene selection](gene_selection/): retained method-specific notebooks.
- [Reference evaluation](reference_evaluation/): shared scLT/scNT workflow.
- [Shared scripts](shared/): cross-figure historical analysis and visualization.
- [Layout](layout/): editorial utilities requiring the original workspace; provenance only.

## Portable saved-data checks

Install the [dependencies](../requirements-demo.txt) in a compatible Python environment. From the repository root:

```bash
python 09_Core/validation/reproduce.py --out /tmp/rna-velocity-check --no-plots
python 09_Core/visualization/plot_saved_panels.py --out /tmp/rna-velocity-panels
python 09_Core/validation/verify_release.py
```

The first command checks saved numerical results; the second draws 17 diagnostic panels, not the formal figure layouts. Neither retrains models. Diagnostic filenames and labels retain historical ED IDs; use the [numbering table](../12_Provenance/figure_index.csv) to find current figures. `--out` and `--no-plots` retain their previous meanings.

[release_paths.py](release_paths.py) resolves historical source names through the [migration manifest](../12_Provenance/path_migration.json). It reads the canonical new locations and does not require the old directory tree or symlinks.
