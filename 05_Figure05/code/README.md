# Figure05: individual panel code

Panel a has a GDC-only numerical check. Its 20-method spatial field layout and schematic remain retained artwork. Panels b-e have separate saved-data redraws. PhyloVelo backward vectors and its non-comparable adapted diagnostic are explicitly preserved.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [figure5a.py](a/figure5a.py) | score-only diagnostic; original fields and schematic retained as artwork |
| b | [figure5b.py](b/figure5b.py) | saved-data redraw; not the formal layout |
| c | [figure5c.py](c/figure5c.py) | saved-data redraw; not the formal layout |
| d | [figure5d.py](d/figure5d.py) | saved-data redraw; not the formal layout |
| e | [figure5e.py](e/figure5e.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 05_Figure05/code/a/figure5a.py --out /tmp/figure05-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
