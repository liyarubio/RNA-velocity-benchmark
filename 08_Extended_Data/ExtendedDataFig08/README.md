# ED08 — scLT lineage reference

Current figure: **ED08**. Historical source ID: **ED14**.

![ED08](<figures/Extended_Data_Fig08.png>)

## Evidence and reproduction

1,756 reference cells and 19 configurations; original membership and random-selection records for the five clone examples are incomplete.

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig08.pdf](<figures/Extended_Data_Fig08.pdf>) | Current publication figure or preview |
| Figures | [Extended_Data_Fig08.png](<figures/Extended_Data_Fig08.png>) | Current publication figure or preview |
| Figures | [Extended_Data_Fig08_legend.pdf](<figures/Extended_Data_Fig08_legend.pdf>) | Current publication figure or preview |
| Data | [Extended_Data_Fig14.xlsx](<data/source_data/Extended_Data_Fig14.xlsx>) | Source workbook, saved results or input metadata |
| Data | [scLT_projected_directions.csv](<data/reference_restoration/scLT_projected_directions.csv>) | Source workbook, saved results or input metadata |
| Code | [reproduce.py](<../../09_Core/validation/reproduce.py>) | Numerical consistency checks on saved data; no model fitting |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | Diagnostic plots of saved data; not replacements for publication layouts |
| Code | [restore_reference_figures.py](<../../09_Core/shared/plots/restore_reference_figures.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [evaluate_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/evaluate_reference.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [finalize.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/finalize.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [prepare.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/prepare.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_all.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/run_all.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [test_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/test_reference.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [restore_tfvelo_projection.py](<../../09_Core/shared/restore_tfvelo_projection.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | Current layout, numbering or caption utilities; no new numerical results |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | Current layout, numbering or caption utilities; no new numerical results |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | Current layout, numbering or caption utilities; no new numerical results |
| Data | [scLT_TFvelo_restored_projection.csv](<data/reference_restoration/scLT_TFvelo_restored_projection.csv>) | Source workbook, saved results or input metadata |
| Data | [scLT_method_summary.csv](<data/reference_restoration/scLT_method_summary.csv>) | Source workbook, saved results or input metadata |
| Data | [scLT_per_cell_cosine.csv](<data/reference_restoration/scLT_per_cell_cosine.csv>) | Source workbook, saved results or input metadata |
| Data | [scLT_reference_cells.csv](<data/reference_restoration/scLT_reference_cells.csv>) | Source workbook, saved results or input metadata |
| Data | [scLT_saved_expected_targets.csv](<data/reference_restoration/scLT_saved_expected_targets.csv>) | Source workbook, saved results or input metadata |
| Data | [scLT_timepoints.csv](<data/reference_restoration/scLT_timepoints.csv>) | Source workbook, saved results or input metadata |

## Additional retained files

- [code/legacy/Visualization/Vis_velocity_stream/Vis_scLT.ipynb](<code/legacy/Visualization/Vis_velocity_stream/Vis_scLT.ipynb>)
- [code/legacy/scLT/1_process/1_run_velocyto.sh](<code/legacy/scLT/1_process/1_run_velocyto.sh>)
- [code/legacy/scLT/1_process/2_add_velocyto.R](<code/legacy/scLT/1_process/2_add_velocyto.R>)
- [code/legacy/scLT/1_process/3_build_cell_pairwise.py](<code/legacy/scLT/1_process/3_build_cell_pairwise.py>)
- [code/legacy/scLT/1_process/4_process_and_bulid_clone_based_velocity.ipynb](<code/legacy/scLT/1_process/4_process_and_bulid_clone_based_velocity.ipynb>)
- [code/legacy/scLT/2_Run_velocity/10_velocyto.ipynb](<code/legacy/scLT/2_Run_velocity/10_velocyto.ipynb>)
- [code/legacy/scLT/2_Run_velocity/11_veloAE.ipynb](<code/legacy/scLT/2_Run_velocity/11_veloAE.ipynb>)
- [code/legacy/scLT/2_Run_velocity/12_veloVAE.ipynb](<code/legacy/scLT/2_Run_velocity/12_veloVAE.ipynb>)
- [code/legacy/scLT/2_Run_velocity/13_STT.ipynb](<code/legacy/scLT/2_Run_velocity/13_STT.ipynb>)
- [code/legacy/scLT/2_Run_velocity/1_scvelo_dynamic.ipynb](<code/legacy/scLT/2_Run_velocity/1_scvelo_dynamic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/1_scvelo_stochastic.ipynb](<code/legacy/scLT/2_Run_velocity/1_scvelo_stochastic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/2_TFvelo.py](<code/legacy/scLT/2_Run_velocity/2_TFvelo.py>)
- [code/legacy/scLT/2_Run_velocity/3_celldancer.ipynb](<code/legacy/scLT/2_Run_velocity/3_celldancer.ipynb>)
- [code/legacy/scLT/2_Run_velocity/4_unitvelo.ipynb](<code/legacy/scLT/2_Run_velocity/4_unitvelo.ipynb>)
- [code/legacy/scLT/2_Run_velocity/5_deepveloSA.ipynb](<code/legacy/scLT/2_Run_velocity/5_deepveloSA.ipynb>)
- [code/legacy/scLT/2_Run_velocity/6_deepveloGB.ipynb](<code/legacy/scLT/2_Run_velocity/6_deepveloGB.ipynb>)
- [code/legacy/scLT/2_Run_velocity/7_veloVI.ipynb](<code/legacy/scLT/2_Run_velocity/7_veloVI.ipynb>)
- [code/legacy/scLT/2_Run_velocity/8_latentvelo.ipynb](<code/legacy/scLT/2_Run_velocity/8_latentvelo.ipynb>)
- [code/legacy/scLT/2_Run_velocity/9_dynamo_deterministic.ipynb](<code/legacy/scLT/2_Run_velocity/9_dynamo_deterministic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/9_dynamo_stochastic.ipynb](<code/legacy/scLT/2_Run_velocity/9_dynamo_stochastic.ipynb>)
- [code/legacy/scLT/3_eval_with_clone.ipynb](<code/legacy/scLT/3_eval_with_clone.ipynb>)

## Related figures

[Figure04](<../../04_Figure04/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
