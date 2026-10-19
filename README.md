# RNA velocity benchmark

Code, saved data and figures for a multilevel benchmark of RNA velocity inference. The repository is organized by **Main Figures 1–7**, followed by **Extended Data Figures 1–12**. Per-figure folders contain actual files; reusable workflows and shared data have separate locations.

## Main figures

| Figure | Topic | Materials |
|---|---|---|
| 1 | 方法时间线 | [Code, data and figures](<01_Figure01/README.md>) |
| 2 | 评价流程 | [Code, data and figures](<02_Figure02/README.md>) |
| 3 | 模拟与常规数据评价 | [Code, data and figures](<03_Figure03/README.md>) |
| 4 | ReDeeM | [Code, data and figures](<04_Figure04/README.md>) |
| 5 | 空间、运行时间及 reads | [Code, data and figures](<05_Figure05/README.md>) |
| 6 | 冻结综合分数和工具 | [Code, data and figures](<06_Figure06/README.md>) |
| 7 | 年龄示例与未来方向 | [Code, data and figures](<07_Figure07/README.md>) |

![Benchmarking framework — Figure 2](02_Figure02/figures/Figure02.png)

## Extended Data, core workflows and tables

- [Extended Data Figures 1–12](08_Extended_Data/README.md)
- [Core workflows: simulation, model execution, gene selection and checks](09_Core/README.md)
- [Shared data](10_Shared/README.md)
- [Tables](11_Tables/README.md)
- [Provenance, inventories and validation](12_Provenance/README.md)
- [Historical materials](13_Archive/README.md)

## Verification and evidence limits

```bash
python -m pip install -r requirements-demo.txt
python 09_Core/validation/reproduce.py --out /tmp/rna-velocity-check --no-plots
python 09_Core/visualization/plot_saved_panels.py --out /tmp/rna-velocity-panels
python 09_Core/validation/verify_release.py
```

Saved-data checks do not constitute model reruns. Historical scripts may require external inputs and their original environments. Manual schematics and retained historical fields are not claimed to be regenerated. Missing values and documented provenance gaps remain unchanged; the frozen overall score is not silently replaced by later source-derived values. See [verification results](12_Provenance/restructure_regression.json).

Simulation datasets, annotated HuBMAP metadata and processed ReDeeM data are available on [Figshare](https://figshare.com/articles/dataset/28021724). Additional dataset details are in the source tables and per-figure records.

The full manuscript and editorial backups are not included. No new blanket license is assigned; see [LICENSE_STATUS.txt](LICENSE_STATUS.txt). The previous structure remains available in [the baseline commit](https://github.com/liyarubio/RNA-velocity-benchmark/tree/2955a8e2ba4d08845699ae429f7fef2dea85d877).
