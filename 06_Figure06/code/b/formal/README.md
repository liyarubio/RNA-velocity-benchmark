# Formal Figure 6b axis clarification

Run `python figure6b_label.py --out OUTPUT_DIR` to export the complete Figure 6 artwork and a standalone panel b in PDF, SVG and PNG. Requires PyMuPDF and NumPy.

The y-axis reads **Mean velocity cosine similarity to full-depth prediction**. The frozen artwork in `assets/Figure06_before.pdf` supplies the original curves and all other panels. The script changes the axis text only and verifies that pixels elsewhere are identical after PDF serialization. Numerical values, colors, method order and original caption pages are preserved. TopicVelo retains its documented integrated-transition displacement proxy in the cosine panel.

The separate manuscript-workspace GDC experiment is a local preview and has not replaced the released cosine panel. No models or metrics are rerun by this label update.
