# KPTracer field layout and displayed GDC order

The formal KPTracer field panel (local ED11b, described by the user as Ext Fig15b) retains the native vector paths and 754 fixed cells per method. All 20 panels use the plotting scale and cell bounds of the original 19-method grid. PhyloVelo is aligned to those same bounds, retaining backward arrows; the frozen vectors are transformed rather than recomputed. The independent entry writes one panel in PDF, SVG and PNG, validation JSON and the displayed-order CSV. Add `--full-figure` for the complete artwork page.

Methods are ordered by the existing displayed GDC, descending, with stable ties. PhyloVelo's +1.000 is an adapted diagnostic and its comparable rank remains missing. Display position is not a claim of comparable superiority. No inputs, scores, model fits, caption text or panel-a artwork change.

Python dependencies are in `../mean_top3/requirements.txt`; SVG export also requires the Poppler `pdftocairo` executable. The renderer reuses the retained `ed16_aligned_fields.py` vector-alignment helper, verifies each method's cells against the saved coordinates by cell type, and checks pixels outside panel b after PDF serialization. The local `compact_extended_1022.py` reconstruction pipeline also applies this layout before final export.
