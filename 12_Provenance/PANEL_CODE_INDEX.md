# Panel code audit: current Figures 1-7 and Extended Data 1-12

Every current portable plotting entrypoint below writes one plot. The batch command only dispatches these files. Figure 3 keeps seven separate R plotting entrypoints and now has five separate upstream metric notebooks for a/b/c/d/f. Numerical data, workbooks and formal figure files are unchanged.

The former Python batch implementation was split into 17 independent plotting blocks; their outputs were checked pixel-for-pixel against the pre-split implementation. Additional saved-data redraws expose available numerical panels separately. Diagnostic encodings are not claims of exact formal-layout reconstruction. Historical model pipelines and manuscript layout assemblers remain provenance records; they require original external inputs and do not serve as the current entrypoints.

## Coverage of every figure

| Figure | Current entrypoints / boundary |
|---|---|
| Figure01 | Schematic; artwork only, no code. |
| Figure02 | Schematic; artwork only, no code. |
| Figure03 | [a-h panel guide](../03_Figure03/code/README.md); e is a schematic; upstream metric notebooks are historical and were not executed. |
| [Figure04](../04_Figure04/code/README.md) | Panel a has a GDC-only numerical check. The complete 21-method stream layout and its schematic remain retained artwork; the released scatter NPZ files do not establish a complete field generator. |
| [Figure05](../05_Figure05/code/README.md) | Panel a has a GDC-only numerical check. Its 20-method spatial field layout and schematic remain retained artwork. Panels b-e have separate saved-data redraws. PhyloVelo backward vectors and its non-comparable adapted diagnostic are explicitly preserved. |
| [Figure06](../06_Figure06/code/README.md) | Panels a, b and c now each have their own Python file. Panel d is the VeloGIF framework schematic and has no numerical plotting script. The historical composite generator is now in code/c; the old Figure5c read-depth generator is split into computation and three separate historical plots under code/b. |
| [Figure07](../07_Figure07/code/README.md) | Historical age-group artwork and future-direction schematic are retained. The exact 784-gene list and matching run have not been recovered; there is no validated current numerical renderer to split. |
| [ED01](../08_Extended_Data/ExtendedDataFig01/code/README.md) | Panel a is the retained taxonomy schematic. Panel b is retained simulation artwork; external simulation files and the original workflow are required for a model-level rerun. No current portable renderer was recovered. |
| [ED02](../08_Extended_Data/ExtendedDataFig02/code/README.md) | Panels a (phase portraits) and b (PCA streams) are retained simulation artwork. Historical notebook sources require external simulation data; no current portable renderer was recovered. |
| [ED03](../08_Extended_Data/ExtendedDataFig03/code/README.md) | Panels a (ten dataset embeddings) and b (HuBMAP examples) are retained artwork. Historical notebooks require external inputs and cannot establish a current numerical reproduction from the packaged metadata alone. |
| [ED04](../08_Extended_Data/ExtendedDataFig04/code/README.md) | Separate a, b and c scripts use current saved data. Panel c preserves all 840 slots, including missing results; each of its six heatmaps has its own color scale. |
| [ED05](../08_Extended_Data/ExtendedDataFig05/code/README.md) | Separate a, b and c scripts cover topology means/sample SDs, paired distributions and projected velocity accuracy. The former pooled pairing diagnostic is retained in its own explicitly named file. |
| [ED06](../08_Extended_Data/ExtendedDataFig06/code/README.md) | Separate a and b files check HuBMAP and classical-dataset GDC. These diagnostics use raw-order heatmaps; they do not reconstruct the formal clustering, organ strip or bubble encodings. |
| [ED07](../08_Extended_Data/ExtendedDataFig07/code/README.md) | Separate a, b and c files cover annotation, reference/method fields and median-cosine diagnostics. The c diagnostic shows medians, not the complete formal boxplot distribution. |
| [ED08](../08_Extended_Data/ExtendedDataFig08/code/README.md) | Separate a, c and d files cover timepoints, reference/method fields and median-cosine diagnostics. Panel b clone examples remain retained artwork because original membership/selection records are incomplete. The d diagnostic shows medians rather than formal boxplots. |
| [ED09](../08_Extended_Data/ExtendedDataFig09/code/README.md) | The current field panel has a separate 20-method saved-vector renderer. The archived GDC score check is separately named; it is not a current lettered panel. |
| [ED10](../08_Extended_Data/ExtendedDataFig10/code/README.md) | Separate a and b files draw cell-type medians and all paired spot maps. The former across-method median diagnostic remains separately named; the boxplot is current Figure 5c. |
| [ED11](../08_Extended_Data/ExtendedDataFig11/code/README.md) | Panel a has an annotation-only renderer; its literature schematic remains artwork. The complete panel b field generator has not been recovered from the packaged inputs. The GDC score check is separate from the current field panels; current score display is Figure 5e. |
| [ED12](../08_Extended_Data/ExtendedDataFig12/code/README.md) | Top100 and top200 overlap/enrichment each have a separate diagnostic file. Overlap checks are intersection matrices, not the original chord drawings. The original chord and heatmap artifacts remain packaged. Bottom-row age examples remain retained artwork with unresolved run provenance. |

## Standalone Python files

| Current figure / panel | File | Rendering scope |
|---|---|---|
| Figure03/a | [figure3a_saved_data.py](../03_Figure03/code/a/figure3a_saved_data.py) | saved-data diagnostic; not the formal layout |
| Figure04/a | [figure4a.py](../04_Figure04/code/a/figure4a.py) | score-only diagnostic; original streams and schematic retained as artwork |
| Figure05/a | [figure5a.py](../05_Figure05/code/a/figure5a.py) | score-only diagnostic; original fields and schematic retained as artwork |
| Figure06/a | [figure6a.py](../06_Figure06/code/a/figure6a.py) | saved-data diagnostic; not the formal layout |
| Figure06/b | [figure6b.py](../06_Figure06/code/b/figure6b.py) | saved-data diagnostic; not the formal layout |
| Figure06/c | [figure6c.py](../06_Figure06/code/c/figure6c.py) | saved-data diagnostic; not the formal layout |
| ED04/a | [ed4a.py](../08_Extended_Data/ExtendedDataFig04/code/a/ed4a.py) | saved-data diagnostic; not the formal layout |
| ED04/b | [ed4b.py](../08_Extended_Data/ExtendedDataFig04/code/b/ed4b.py) | saved-data diagnostic; not the formal layout |
| ED05/paired_check | [paired_check.py](../08_Extended_Data/ExtendedDataFig05/code/paired_check/paired_check.py) | pooled diagnostic only; formal a-c are separate |
| ED06/a | [ed6a.py](../08_Extended_Data/ExtendedDataFig06/code/a/ed6a.py) | raw-order heatmap diagnostic; does not reproduce formal clustering or organ strip |
| ED07/c | [ed07c.py](../08_Extended_Data/ExtendedDataFig07/code/c/ed07c.py) | median-only diagnostic; formal boxplots retained as artwork |
| ED08/d | [ed08d.py](../08_Extended_Data/ExtendedDataFig08/code/d/ed08d.py) | median-only diagnostic; formal boxplots retained as artwork |
| ED09/score_check | [ed09_scores.py](../08_Extended_Data/ExtendedDataFig09/code/score_check/ed09_scores.py) | archived score diagnostic; separate from current field panels |
| ED11/score_check | [ed11_scores.py](../08_Extended_Data/ExtendedDataFig11/code/score_check/ed11_scores.py) | archived score diagnostic; separate from current field panels |
| ED10/median_check | [median_check.py](../08_Extended_Data/ExtendedDataFig10/code/median_check/median_check.py) | median-only diagnostic; formal a-b are separate |
| ED12/top100_overlap | [top100_overlap.py](../08_Extended_Data/ExtendedDataFig12/code/top100_overlap/top100_overlap.py) | intersection-matrix diagnostic; original chord artwork retained |
| ED12/top200_overlap | [top200_overlap.py](../08_Extended_Data/ExtendedDataFig12/code/top200_overlap/top200_overlap.py) | intersection-matrix diagnostic; original chord artwork retained |
| Figure04/b | [figure4b.py](../04_Figure04/code/b/figure4b.py) | saved-data redraw; not the formal layout |
| Figure05/b | [figure5b.py](../05_Figure05/code/b/figure5b.py) | saved-data redraw; not the formal layout |
| Figure05/c | [figure5c.py](../05_Figure05/code/c/figure5c.py) | saved-data redraw; not the formal layout |
| Figure05/d | [figure5d.py](../05_Figure05/code/d/figure5d.py) | saved-data redraw; not the formal layout |
| Figure05/e | [figure5e.py](../05_Figure05/code/e/figure5e.py) | saved-data redraw; not the formal layout |
| ED04/c | [ed4c.py](../08_Extended_Data/ExtendedDataFig04/code/c/ed4c.py) | saved-data redraw; not the formal layout |
| ED05/a | [ed5a.py](../08_Extended_Data/ExtendedDataFig05/code/a/ed5a.py) | saved-data redraw; not the formal layout |
| ED05/b | [ed5b.py](../08_Extended_Data/ExtendedDataFig05/code/b/ed5b.py) | saved-data redraw; not the formal layout |
| ED05/c | [ed5c.py](../08_Extended_Data/ExtendedDataFig05/code/c/ed5c.py) | saved-data redraw; not the formal layout |
| ED06/b | [ed6b.py](../08_Extended_Data/ExtendedDataFig06/code/b/ed6b.py) | saved-data redraw; not the formal layout |
| ED07/a | [ed7a.py](../08_Extended_Data/ExtendedDataFig07/code/a/ed7a.py) | saved-data redraw; not the formal layout |
| ED07/b | [ed07b.py](../08_Extended_Data/ExtendedDataFig07/code/b/ed07b.py) | saved-data redraw; not the formal layout |
| ED08/c | [ed08c.py](../08_Extended_Data/ExtendedDataFig08/code/c/ed08c.py) | saved-data redraw; not the formal layout |
| ED08/a | [ed8a.py](../08_Extended_Data/ExtendedDataFig08/code/a/ed8a.py) | saved-data redraw; not the formal layout |
| ED10/a | [ed10a.py](../08_Extended_Data/ExtendedDataFig10/code/a/ed10a.py) | saved-data redraw; not the formal layout |
| ED10/b | [ed10b.py](../08_Extended_Data/ExtendedDataFig10/code/b/ed10b.py) | saved-data redraw; not the formal layout |
| ED12/top100_enrichment | [top100_enrichment.py](../08_Extended_Data/ExtendedDataFig12/code/top100_enrichment/top100_enrichment.py) | saved-data redraw; not the formal layout |
| ED12/top200_enrichment | [top200_enrichment.py](../08_Extended_Data/ExtendedDataFig12/code/top200_enrichment/top200_enrichment.py) | saved-data redraw; not the formal layout |
| ED09/fields | [ed9_fields.py](../08_Extended_Data/ExtendedDataFig09/code/fields/ed9_fields.py) | saved-vector redraw; not the formal layout |
| ED11/a | [ed11a.py](../08_Extended_Data/ExtendedDataFig11/code/a/ed11a.py) | annotation-only redraw; literature schematic retained as artwork |

## Run one panel or all available checks

```bash
python 06_Figure06/code/b/figure6b.py --out /tmp/figure6b
python 09_Core/visualization/plot_saved_panels.py --list
python 09_Core/visualization/plot_saved_panels.py --panel Figure06/b --out /tmp/figure6b
python 09_Core/visualization/plot_saved_panels.py --out /tmp/all-panels
```

Most original diagnostic output stems retain historical IDs for compatibility. ED4a/b use current IDs to avoid collisions with current ED5a/b. The registry and per-output records identify the **current** figure/panel; worksheet names also retain their recorded historical numbering. Missing data remain missing.

## Validation

All 37 current Python entrypoints were executed independently and each produced one PDF/SVG plot plus a source record. All 17 original diagnostic outputs are pixel-identical to the former batch implementation. The three separated historical read-depth plot files were also executed using the packaged pairwise values. R syntax and the six split notebook structures passed; historical R drawing and external-H5AD metric computation were not rerun. All 352 existing data/artwork files retained their hashes. [Machine-readable validation](panel_split_validation.json).

## Formal Mean-only panels (23 October 2026)

The current formal vector entrypoints are [Figure 3a](../03_Figure03/code/a/figure3a_formal.py), [ED4a](../08_Extended_Data/ExtendedDataFig04/code/a/ed4a_formal.py), and [ED4b](../08_Extended_Data/ExtendedDataFig04/code/b/ed4b_formal.py). Each writes one PDF/SVG/PNG panel with only Mean and aligned ranks 1–3. The existing saved-data diagnostics remain separate. [Validation](mean_top3_validation.json).

Figure 3h now has an [independent formal vector entrypoint](../03_Figure03/code/h/figure3h_formal.py) with descending available-case mean order and a color-coded Mean row; the best-performance envelope remains last. [Order validation](figure3h_mean_validation.json) and [Mean-row validation](figure3h_mean_color_validation.json).

ED08a now has an [independent formal entrypoint](../08_Extended_Data/ExtendedDataFig08/code/a/ed8a.py) using the ED08c cell-type palette. The 3,186 point positions and panels b-d are retained. [Color and geometry validation](ed08a_celltype_validation.json).

ED11b now has an [independent formal vector-layout entrypoint](../08_Extended_Data/ExtendedDataFig11/code/b/ed11b.py): 20 equal-size fields, common 754-cell bounds and descending displayed GDC. [Validation](ed11_field_layout_validation.json). PhyloVelo keeps its native backward arrows and non-comparable adapted-score status.
