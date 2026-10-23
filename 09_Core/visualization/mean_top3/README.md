# Mean-only formal vector renderer

The three independent panel entrypoints reproduce the published presentation edit from frozen vector templates and the current numerical inputs. They remove the Median display column, shift Mean into its position and display only ranks 1–3 at a fixed horizontal anchor. Ties use minimum rank. The original joint Mean/Median normalization remains in force; Median is calculated only to verify that scale.

The unchanged matrix, annotations and legend are retained as vectors. Mean bar widths are checked against recomputed means, and the original fills and line widths are preserved. The templates are layout inputs, not newly computed benchmark results. Editing the numerical matrix beyond this presentation revision requires rebuilding the full matrix renderer; input/width mismatches fail instead of silently producing stale plots.

Install `requirements.txt`. Each entrypoint accepts `--out DIR` and writes one panel as PDF, SVG and PNG, plus a JSON validation record. The publication compositor combines the updated panels with the retained full figures. Manuscript text, caption text and source spreadsheets are outside this change.

To export both complete artwork pages, run `python 09_Core/visualization/mean_top3/render_figures.py --out demo_output/mean_top3` from the repository root. Companion caption pages are not generated or edited by this command.

The input-hash registry accepts the local and English-release workbook containers, whose two plotted numerical tables were checked for exact equality.
