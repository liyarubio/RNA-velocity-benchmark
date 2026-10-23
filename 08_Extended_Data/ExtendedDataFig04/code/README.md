# ED04: individual panel code

Separate a, b and c scripts use current saved data. Panel c preserves all 840 slots, including missing results; each of its six heatmaps has its own color scale.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [ed4a.py](a/ed4a.py) | saved-data diagnostic; not the formal layout |
| b | [ed4b.py](b/ed4b.py) | saved-data diagnostic; not the formal layout |
| c | [ed4c.py](c/ed4c.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig04/code/a/ed4a.py --out /tmp/ed04-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)

## Current formal panels a and b (23 October 2026)

| Panel | Formal vector entrypoint |
|---|---|
| a | [ed4a_formal.py](a/ed4a_formal.py) |
| b | [ed4b_formal.py](b/ed4b_formal.py) |

```bash
python -m pip install -r 09_Core/visualization/mean_top3/requirements.txt
python 08_Extended_Data/ExtendedDataFig04/code/a/ed4a_formal.py --out demo_output/mean_top3
python 08_Extended_Data/ExtendedDataFig04/code/b/ed4b_formal.py --out demo_output/mean_top3
```

Each formal entrypoint writes one PDF/SVG/PNG panel and validation JSON. The numerical inputs and original bar normalization are preserved; only Mean and ranks 1–3 remain visible. Rank digits share a fixed horizontal position inside the bars. The unchanged matrix and annotations come from the retained vector template. Panel c is retained exactly. [Renderer scope](../../../09_Core/visualization/mean_top3/README.md).
