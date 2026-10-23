# ED12: individual panel code

Top100 and top200 overlap/enrichment each have a separate diagnostic file. Overlap checks are intersection matrices, not the original chord drawings. The original chord and heatmap artifacts remain packaged. Bottom-row age examples remain retained artwork with unresolved run provenance.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| top100_overlap | [top100_overlap.py](top100_overlap/top100_overlap.py) | intersection-matrix diagnostic; original chord artwork retained |
| top200_overlap | [top200_overlap.py](top200_overlap/top200_overlap.py) | intersection-matrix diagnostic; original chord artwork retained |
| top100_enrichment | [top100_enrichment.py](top100_enrichment/top100_enrichment.py) | saved-data redraw; not the formal layout |
| top200_enrichment | [top200_enrichment.py](top200_enrichment/top200_enrichment.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig12/code/top100_overlap/top100_overlap.py --out /tmp/ed12-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
