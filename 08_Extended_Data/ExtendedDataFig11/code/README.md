# ED11: individual panel code

Panel a has an annotation-only renderer; its literature schematic remains artwork. The complete panel b field generator has not been recovered from the packaged inputs. The GDC score check is separate from the current field panels; current score display is Figure 5e.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| score_check | [ed11_scores.py](score_check/ed11_scores.py) | archived score diagnostic; separate from current field panels |
| a | [ed11a.py](a/ed11a.py) | annotation-only redraw; literature schematic retained as artwork |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig11/code/score_check/ed11_scores.py --out /tmp/ed11-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
