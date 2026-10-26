# Figure06: individual panel code

Panels a, b and c now each have their own Python file. Panel d is the VeloGIF framework schematic and has no numerical plotting script. The old Figure6a composite generator is now [c/historical_composite.py](c/historical_composite.py). The old Figure5c read-depth script has been split into [one computation file and three separate historical plot files](b/historical/README.md).

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [figure6a.py](a/figure6a.py) | saved-data diagnostic; not the formal layout |
| b | [figure6b.py](b/figure6b.py) | saved-data diagnostic; not the formal layout |
| b formal | [figure6b_label.py](b/formal/figure6b_label.py) | formal artwork with explicit full-depth reference on the y-axis |
| c | [figure6c.py](c/figure6c.py) | saved-data weighted-score diagnostic |
| c formal | [figure6c.py](c/formal/figure6c.py) | formal weighted composite; 0.50 / 0.25 / 0.25 |

Run one file from the repository root:

```bash
python 06_Figure06/code/a/figure6a.py --out /tmp/figure06-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

The [formal panel-b entry](b/formal/README.md) requires PyMuPDF and NumPy. Its axis reads “Mean velocity cosine similarity to full-depth prediction”; all saved curves, other panels, and caption text are preserved. It also exports PNG and checks that only the axis-label pixels changed.

[Complete panel audit](../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
