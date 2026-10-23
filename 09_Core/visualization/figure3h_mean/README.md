# Figure 3h: descending mean GDC

This renderer permutes the current formal vector columns and their method-color strips using means recomputed from the saved numerical table. It preserves the dataset row order, cell colors, missing cells, and the right-hand best-performance envelope. Missing values are excluded from the mean, never replaced with zero; ties retain the earlier display order.

The frozen vector input preserves the current Figure 3a Mean-only/top-three update. The data hash is checked before rendering. This is a presentation edit, not a fresh benchmark run. Use the shared `../mean_top3/requirements.txt` dependencies. The independent panel entry accepts `--out DIR`; `--full-figure` writes the full Figure 3 artwork instead. Both produce `method_order.csv` and validation JSON. The combined `../mean_top3/render_figures.py` pipeline also applies this panel ordering.
