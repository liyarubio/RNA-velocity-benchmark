# Manuscript_1021 public release — 23 October 2026

This release imports the saved results and final figures in Manuscript_1021, including its 22 October updates. It replaces the 19 October public snapshot in the existing figure-organized layout. The manuscript workspace is unchanged. No velocity model was fitted during publication.

## Schematic figure organization

Figures 1 and 2 contain only the current PDF, SVG and PNG plus a navigation README. Figure-specific code, historical backup records and data links were removed on 23 October at the author's request. Shared tables remain in the tables directory. The six current figure files are unchanged.

## Scientific changes and limits

| Material | Current result | Interpretation boundary |
|---|---|---|
| TopoVelo / spatial | Archived ST_small and ST_raw predictions; 608 paired spots, 47 shared genes, 606 valid pairs | Training settings were absent from the archive. Two undefined paired directions remain missing. Public checks use saved aggregates; original input H5ADs are external. |
| PhyloVelo / KPTracer | 754 cells, 706 reference cells; native backward field for display; negative native vectors for metrics | The global gene-wise velocity vector produces ICCoh structurally equal to one. Adapted GDC +1 is diagnostic; comparable GDC and rank are null. Original inference predates assembly. |
| cellDancer / ED4c | 30 scenario/noise conditions; 840 slots, 822 finite and 18 missing; 84 cross-noise means | Means equally weight available dataset medians. STT forward branching has 1/10 coverage in both levels; other groups have 10/10. These are conditions, not independent simulation replicates. Six heatmap scales are independent. |
| Figure 4 | Larger markers with preserved coordinates and nine recovered scatter canvases | Presentation update; no metric recomputation. Historical streamline provenance limits remain. |
| Figure 6 | a runtime; b read-depth; c weighted composite; d VeloGIF | Composite weights are provisionally 0.50 Accuracy, 0.25 Scalability and 0.25 Stability; endpoint inputs remain unchanged. Caption-formula revision is pending with the author. |
| Age-related examples | Retained historical material | Exact 784-gene list and matching run remain unresolved; no new aging conclusion. |

## Current files and historical records

- [Figure index](figure_index.csv) maps current and historical numbering.
- [File inventory](file_inventory.csv) links figures to actual repository files.
- [Current source resolver](current_sources.json) binds original input names to public paths.
- [Source import ledger](release_1021_sources.json) records original/current hashes and explicit exclusions.
- [Checksums](current_checksums.json), [structure validation](release_validation.json), [latest-result validation](latest_results_validation.json), [PhyloVelo validation](phylovelo_validation.json), [saved-data validation](saved_data_validation.json), and [English validation](english_validation.json) describe this release.
- [The 1019 ledgers](baseline_1019/) and inherited reports retain their original scope; their hashes describe the previous tree and are not current certificates.

Final PDFs and PNGs are copied from the current source package. Source-workbook numerical cells and formulas are preserved; inherited Chinese annotations/style names are translated or reduced to their existing English counterpart. Duplicate old ED filenames and obsolete SVG variants are removed from active figure folders. Large TopoVelo input H5ADs, raw PhyloVelo downloads, pickles, upstream mirrors, historical reports and manuscript assembly utilities are excluded; available hashes, sizes and reasons are listed in the source ledger. The two compact PhyloVelo result H5ADs are included for the independent direction audit.

The repository remains English. Public release is not a claim that the manuscript is submission-ready; the source package still records deferred abstract and pending author verification. License status is unchanged.

## Release tooling attribution

The scientific-writing skill was used for evidence and reporting boundaries during release preparation. Tooling reference: Kassis, T., Agarwal, V., He, Y., Patel, D., and Brueckner, A. M. (2026). [Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065). Metadata checked against arXiv on 23 October 2026. This is a tooling attribution, not a scientific evidence source for benchmark results.
