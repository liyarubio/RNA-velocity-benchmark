# Core workflows and validation

Run commands from the repository root with Python 3.9 or later and `requirements-demo.txt`. The current source resolver reads [current_sources.json](../12_Provenance/current_sources.json); it does not require the manuscript directory.

```bash
python 09_Core/validation/verify_release.py
python 09_Core/validation/reproduce.py --out /tmp/rna-velocity-check --no-plots
python 09_Core/validation/validate_latest_results.py --out /tmp/rna-velocity-1021
python 09_Core/visualization/plot_saved_panels.py --out /tmp/rna-velocity-panels
```

`reproduce.py` checks the saved benchmark and unchanged composite, including the expanded 20-method paired-context table. `validate_latest_results.py` independently recomputes noise means, TopoVelo GDC from saved cluster/edge summaries, and coverage. `plot_saved_panels.py` is now only a dispatcher for 37 independent one-plot files. Use `--list` or `--panel Figure06/b` to locate or run one. The [complete panel audit](../12_Provenance/PANEL_CODE_INDEX.md) distinguishes full numerical redraws, partial diagnostics and retained artwork. Most original diagnostic filenames retain historical IDs; ED4a/b use current IDs to prevent output collisions.

```bash
python -m pip install -r requirements-audit.txt
python 09_Core/validation/validate_phylovelo.py --out /tmp/rna-velocity-1021
python 09_Core/validation/check_english.py --pdf --out /tmp/rna-velocity-english.json
```

The PhyloVelo check reads the two packaged result H5ADs, verifies vector sign and cell identity, and recalculates boundary-cell cosines. It does not fit a model. TopoVelo archived H5AD inputs are external; the default release check verifies saved aggregates, not gene-level predictions or training reproducibility.

## Historical workflows

- [dyngen simulation](simulation/dyngen/) and [cellDancer simulation](simulation/cellDancer/)
- [Model execution](model_execution/) and [gene selection](gene_selection/)
- [Reference evaluation](reference_evaluation/) and [shared scripts](shared/)
- [TopoVelo projection and metric implementation](shared/latest/)
- [Historical layout utilities](layout/)

Historical scripts may contain original filesystem paths and require external inputs, model environments or the manuscript workspace. They are retained for provenance and are not portable public entrypoints. No model was fitted during this release assembly. The PhyloVelo inference was performed earlier on 21 October 2026 and imported here.
