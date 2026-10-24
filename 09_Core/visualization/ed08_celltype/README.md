# ED8a cell-type colors

The independent renderer preserves all 3,186 original point positions, plotting order, size and opacity, grouped as Day 2.5 and Day 5. Cell IDs join the retained coordinates to annotations extracted from `LSK_lineage.h5ad`, field `obs/cell_type2`, using the same grouping rules as ED8c. All 1,756 shared labels match ED8c exactly.

The seven shared types use the exact ED8c palette. Day 5 additionally contains six Mast cells and two Ccr7_DC cells, displayed together in gray as Other, as requested. Their original annotations remain unchanged in the source CSV. No model fitting, embedding, scores or workbook values are changed. The original sample legends are replaced by a cell-type legend. The formal panels b, c and d and caption text are retained.

Run the panel entry with `--out DIR`; add `--full-figure` for the complete artwork page. The renderer uses the shared `../mean_top3/requirements.txt` dependencies and verifies point-to-cell alignment, colors, geometry and pixels outside panel a after PDF serialization. The frozen vector template is the formal A4 page before this edit.
