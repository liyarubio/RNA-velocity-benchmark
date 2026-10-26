# Figure 6c weighted composite

Provisional weights: **Accuracy 0.50, Scalability 0.25, Stability 0.25**. Each component uses its ascending average rank divided by 19. Overall score is `(2 * accuracy_rank + scalability_rank + stability_rank) / (4 * 19)`, computed with integer weights to preserve ties. Overall ranks use descending minimum rank, as before.

Accuracy retains the eight available-domain aggregates, scalability retains reciprocal runtimes at 1K and 10K cells, and stability retains mean signed velocity cosine at 0.2/0.4/0.6/0.8x reads relative to full depth, excluding self-comparison. The 19 configurations and all endpoint values are unchanged; later endpoint updates are not silently inserted. Historical input/protocol differences remain.

Run `python figure6c.py --out OUTPUT_DIR` to export the formal panel, numerical summary and validation. Requires PyMuPDF, pandas, NumPy, Matplotlib, openpyxl, and the retained shared `mean_top3` vector helper. The renderer preserves the original method order and all content outside the total-score column and panel title, checked after PDF serialization. Top-six rank annotations and the existing color convention are retained.

The manuscript publisher updates only Figure 6 artwork and derived workbook totals/ranks plus their definitions. Captions and manuscript prose remain author-managed: their previous 60:40 formula needs a corresponding editorial update. `assets/summary_before.csv` and `assets/Figure06_before.pdf` preserve the preceding state for verification.
