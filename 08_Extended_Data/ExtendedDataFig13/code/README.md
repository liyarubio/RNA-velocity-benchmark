# Extended Data 13: GDC across read depths

The single panel shows 21 methods at five read-depth fractions (0.2, 0.4, 0.6, 0.8 and 1.0) for ReDeeM `young2_HSPC_Multi`. The 105 GDC values come directly from the five saved `Eval.csv` files; no models or metrics are rerun. Original method names and hashes are retained.

Method colors and legend order are extracted from the current formal Figure 6b vector legend. The script verifies each plotted line's RGB and values, then checks 8-bit RGB in the serialized PDF. TopicVelo retains its main-figure color but uses a dashed line and hollow square to identify the archived base-scVelo proxy, which is not its native integrated-transition result. The full-depth GDC is an evaluated score rather than a self-comparison.

Run `python ed13a.py --out OUTPUT_DIR`. Outputs: portrait A4 PDF, SVG, PNG, plotted values, method-color table and validation JSON. Dependencies: Python, NumPy, pandas, Matplotlib, PyMuPDF and the retained shared `mean_top3` vector-clipping helper. Inputs are located automatically in the manuscript workspace or public release tree.

ED13 is appended to preserve the existing figure numbers. Its first manuscript citation and full caption remain author-managed. The artwork includes the essential proxy interpretation note; no existing caption or manuscript paragraph is rewritten. The new workbook is named `Extended_Data_Fig13_reads_GDC.xlsx` to distinguish it from the legacy-numbered scNT workbook `Extended_Data_Fig13.xlsx`.
