# Figure04: individual panel code

Panel a has a GDC-only numerical check. The complete 21-method stream layout and its schematic remain retained artwork; the released scatter NPZ files do not establish a complete field generator.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [figure4a.py](a/figure4a.py) | score-only diagnostic; original streams and schematic retained as artwork |
| b | [figure4b.py](b/figure4b.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 04_Figure04/code/a/figure4a.py --out /tmp/figure04-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
