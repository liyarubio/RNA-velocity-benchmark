# ED05: individual panel code

Separate a, b and c scripts cover topology means/sample SDs, paired distributions and projected velocity accuracy. The former pooled pairing diagnostic is retained in its own explicitly named file.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| paired_check | [paired_check.py](paired_check/paired_check.py) | pooled diagnostic only; formal a-c are separate |
| a | [ed5a.py](a/ed5a.py) | saved-data redraw; not the formal layout |
| b | [ed5b.py](b/ed5b.py) | saved-data redraw; not the formal layout |
| c | [ed5c.py](c/ed5c.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig05/code/paired_check/paired_check.py --out /tmp/ed05-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
