# ED11 — KPTracer and PhyloVelo direction audit

Current release: **Manuscript_1021**, including updates through **22 October 2026**. Historical numerical source ID: **ED17**.

754 cells and 706 reference cells. PhyloVelo uses native backward vectors for display and their negative for metrics. A repeated global gene vector makes ICCoh structurally one: adapted GDC +1 is not a comparable rank.

![ED11](<figures/Extended_Data_Fig11.png>)

## Individual panel code

[One plot per file: current entrypoints and limits](code/README.md).

## Current files

| Role | File |
|---|---|
| results | [figures/Extended_Data_Fig11.pdf](<figures/Extended_Data_Fig11.pdf>) |
| results | [figures/Extended_Data_Fig11.png](<figures/Extended_Data_Fig11.png>) |
| results | [figures/Extended_Data_Fig11.svg](<figures/Extended_Data_Fig11.svg>) |
| results | [figures/legends/Extended_Data_Fig11_legend.pdf](<figures/legends/Extended_Data_Fig11_legend.pdf>) |
| data | [data/source_data/Extended_Data_Fig17.xlsx](<data/source_data/Extended_Data_Fig17.xlsx>) |
| data | [data/kptracer_completion/GDC_display_order.csv](<data/kptracer_completion/GDC_display_order.csv>) |
| data | [data/kptracer_completion/code_manifest.json](<data/kptracer_completion/code_manifest.json>) |
| data | [data/kptracer_completion/completion_audit.json](<data/kptracer_completion/completion_audit.json>) |
| data | [data/kptracer_completion/control.json](<data/kptracer_completion/control.json>) |
| data | [data/kptracer_completion/control_directions.csv](<data/kptracer_completion/control_directions.csv>) |
| data | [data/kptracer_completion/control_unit_scores.csv](<data/kptracer_completion/control_unit_scores.csv>) |
| data | [data/kptracer_completion/execution.json](<data/kptracer_completion/execution.json>) |
| data | [data/kptracer_completion/input_manifest.json](<data/kptracer_completion/input_manifest.json>) |
| data | [data/kptracer_completion/native_topicvelo_transition.npz](<data/kptracer_completion/native_topicvelo_transition.npz>) |
| data | [data/kptracer_completion/tfvelo.json](<data/kptracer_completion/tfvelo.json>) |
| data | [data/kptracer_completion/tfvelo_directions.csv](<data/kptracer_completion/tfvelo_directions.csv>) |
| data | [data/kptracer_completion/tfvelo_gene_mapping.csv](<data/kptracer_completion/tfvelo_gene_mapping.csv>) |
| data | [data/kptracer_completion/tfvelo_resources.json](<data/kptracer_completion/tfvelo_resources.json>) |
| data | [data/kptracer_completion/tfvelo_revision.json](<data/kptracer_completion/tfvelo_revision.json>) |
| data | [data/kptracer_completion/tfvelo_unit_scores.csv](<data/kptracer_completion/tfvelo_unit_scores.csv>) |
| data | [data/kptracer_completion/tfvelo_validation.json](<data/kptracer_completion/tfvelo_validation.json>) |
| data | [data/kptracer_completion/topicvelo.json](<data/kptracer_completion/topicvelo.json>) |
| data | [data/kptracer_completion/topicvelo_directions.csv](<data/kptracer_completion/topicvelo_directions.csv>) |
| data | [data/kptracer_completion/topicvelo_inference_status.json](<data/kptracer_completion/topicvelo_inference_status.json>) |
| data | [data/kptracer_completion/topicvelo_unit_scores.csv](<data/kptracer_completion/topicvelo_unit_scores.csv>) |
| data | [data/kptracer_completion/topicvelo_validation.json](<data/kptracer_completion/topicvelo_validation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/README.md](<../../10_Shared/data/kptracer_phylovelo_20261021/README.md>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/CBDir_independent_cell_scores.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/CBDir_independent_cell_scores.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/CBDir_independent_summary.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/CBDir_independent_summary.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/GDC_cluster_contributions.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/GDC_cluster_contributions.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/ICCoh_neighbor_coverage.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/ICCoh_neighbor_coverage.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/audit.json](<../../10_Shared/data/kptracer_phylovelo_20261021/calculation_audit_20261021/audit.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/evaluation_metrics.py](<../../10_Shared/data/kptracer_phylovelo_20261021/code/evaluation_metrics.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/fetch_inputs.py](<../../10_Shared/data/kptracer_phylovelo_20261021/code/fetch_inputs.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/fetch_ranges.py](<../../10_Shared/data/kptracer_phylovelo_20261021/code/fetch_ranges.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/plot_results.py](<../../10_Shared/data/kptracer_phylovelo_20261021/code/plot_results.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/run_docker.sh](<../../10_Shared/data/kptracer_phylovelo_20261021/code/run_docker.sh>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/code/run_phylovelo.py](<../../10_Shared/data/kptracer_phylovelo_20261021/code/run_phylovelo.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/CBDir_orientation_comparison.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/CBDir_orientation_comparison.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/CBDir_unit_scores.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/CBDir_unit_scores.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/PhyloVelo_direction_comparison.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/PhyloVelo_direction_comparison.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/PhyloVelo_direction_comparison.png](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/PhyloVelo_direction_comparison.png>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/directions_both_conventions.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/directions_both_conventions.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/evaluation_corrected.json](<../../10_Shared/data/kptracer_phylovelo_20261021/direction_audit_20261021/evaluation_corrected.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/evaluation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/evaluation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/external_inputs.json](<../../10_Shared/data/kptracer_phylovelo_20261021/external_inputs.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/inputs/3726_NT_T1_tree.nwk](<../../10_Shared/data/kptracer_phylovelo_20261021/inputs/3726_NT_T1_tree.nwk>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/inputs/cell_type_palette.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/inputs/cell_type_palette.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/inputs/current_GDC.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/inputs/current_GDC.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/CBDir_forward.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/CBDir_forward.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/CBDir_forward_unit_scores.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/CBDir_forward_unit_scores.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_manuscript.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_manuscript.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_manuscript.png](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_manuscript.png>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_panel.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_panel.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_tutorial.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_tutorial.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_tutorial.png](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/PhyloVelo_backward_tutorial.png>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/directions.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/directions.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/evaluation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/evaluation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/panel_style.json](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/panel_style.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/phylovelo_forward_metrics.h5ad](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/phylovelo_forward_metrics.h5ad>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/render_panel.py](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/render_panel.py>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/validation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/presentation_backward_metrics_forward/validation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/code_manifest.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/code_manifest.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/control_validation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/control_validation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/docker_image_id.txt](<../../10_Shared/data/kptracer_phylovelo_20261021/records/docker_image_id.txt>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/download.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/download.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/environment.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/environment.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/input_audit.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/input_audit.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/tutorial_selection_validation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/tutorial_selection_validation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/records/validation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/records/validation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/GDC_with_phylovelo.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/results/GDC_with_phylovelo.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.png](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.png>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.svg](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_manuscript.svg>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.pdf](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.pdf>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.png](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.png>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.svg](<../../10_Shared/data/kptracer_phylovelo_20261021/results/PhyloVelo_KPTracer_tutorial.svg>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/cell_tree_alignment.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/results/cell_tree_alignment.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/evaluation.json](<../../10_Shared/data/kptracer_phylovelo_20261021/results/evaluation.json>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/evaluation_unit_scores.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/results/evaluation_unit_scores.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/gene_velocities.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/results/gene_velocities.csv>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/latent_expression.csv.gz](<../../10_Shared/data/kptracer_phylovelo_20261021/results/latent_expression.csv.gz>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/phylovelo.h5ad](<../../10_Shared/data/kptracer_phylovelo_20261021/results/phylovelo.h5ad>) |
| data | [../../10_Shared/data/kptracer_phylovelo_20261021/results/phylovelo_directions.csv](<../../10_Shared/data/kptracer_phylovelo_20261021/results/phylovelo_directions.csv>) |
| data | [../../10_Shared/data/figure5_1021/DeepVelo(VAE)_plot_vectors.csv](<../../10_Shared/data/figure5_1021/DeepVelo(VAE)_plot_vectors.csv>) |
| data | [../../10_Shared/data/figure5_1021/KPTracer_GDC_with_PhyloVelo.csv](<../../10_Shared/data/figure5_1021/KPTracer_GDC_with_PhyloVelo.csv>) |
| data | [../../10_Shared/data/figure5_1021/PhyloVelo_plot_vectors.csv](<../../10_Shared/data/figure5_1021/PhyloVelo_plot_vectors.csv>) |
| data | [../../10_Shared/data/figure5_1021/UniTVelo_plot_vectors.csv](<../../10_Shared/data/figure5_1021/UniTVelo_plot_vectors.csv>) |
| data | [../../10_Shared/data/figure5_1021/field_style.json](<../../10_Shared/data/figure5_1021/field_style.json>) |
| data | [../../10_Shared/data/figure5_1021/method_colors.json](<../../10_Shared/data/figure5_1021/method_colors.json>) |
| data | [../../10_Shared/data/figure5_1021/panel_map.json](<../../10_Shared/data/figure5_1021/panel_map.json>) |
| code | [code/historical/kptracer_completion/complete_topicvelo.py](<code/historical/kptracer_completion/complete_topicvelo.py>) |
| code | [code/historical/kptracer_completion/config.py](<code/historical/kptracer_completion/config.py>) |
| code | [code/historical/kptracer_completion/evaluate.py](<code/historical/kptracer_completion/evaluate.py>) |
| code | [code/historical/kptracer_completion/evaluation_metrics.py](<code/historical/kptracer_completion/evaluation_metrics.py>) |
| code | [code/historical/kptracer_completion/import_to_manuscript.py](<code/historical/kptracer_completion/import_to_manuscript.py>) |
| code | [code/historical/kptracer_completion/run_tfvelo.py](<code/historical/kptracer_completion/run_tfvelo.py>) |
| code | [code/historical/kptracer_completion/validate_tfvelo.py](<code/historical/kptracer_completion/validate_tfvelo.py>) |
| code | [code/historical/plots/redraw_ed16_reference.py](<code/historical/plots/redraw_ed16_reference.py>) |
| code | [code/historical/plots/ed16_literature_schematic.py](<code/historical/plots/ed16_literature_schematic.py>) |
| code | [code/historical/plots/ed16_aligned_fields.py](<code/historical/plots/ed16_aligned_fields.py>) |
| code | [code/historical/plots/complete_ed16_fields.py](<code/historical/plots/complete_ed16_fields.py>) |
| code | [code/historical/plots/spatial_display_names.py](<code/historical/plots/spatial_display_names.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../../README.md>) · [Release notes](<../../12_Provenance/RELEASE_1021.md>)
