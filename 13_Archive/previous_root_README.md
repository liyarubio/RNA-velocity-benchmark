# RNA velocity benchmark

Code and data for evaluating RNA velocity inference across prediction spaces, datasets and biological references.

## Current figure-ordered materials — 19 October 2026

The [1019 code/data release](Manuscript_1019/README.md) provides **7 main figures and 12 Extended Data figures** in first-citation order, with source workbooks, saved numerical results, available implementations and current PDF/PNG figures.

- [Browse code and data by figure](Manuscript_1019/Figure_materials/README.md)
- [Main and supplementary tables](Manuscript_1019/Figure_materials/03_Tables/README.md)
- [File inventory and source hashes](Manuscript_1019/Figure_materials/file_inventory.csv)
- [Verification and reproduction limits](Manuscript_1019/README.md#saved-data-verification)

![Benchmarking framework](Manuscript_1019/figures/main/Figure02.png)

The update preserves the earlier analysis directories below. Saved-data checks are distinct from full model reruns; incomplete provenance and nonuniform historical preprocessing remain documented in the per-figure READMEs. No new blanket license is assigned.

## Earlier analysis directories and dataset access


### Benchmarking framework and datasets
The simulation datasets, the annotated metadata of HuBMAP，and processed Redeem dataset are uploaded on [figshare](https://figshare.com/articles/dataset/28021724)

Other processed real-world datasets are download from published work (more detailed in supplementary table3 of our manuscript), and the visualization of datasets are summarized in [Datasets visualization](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Visualization/Vis_datasets)


### Accuracy in conventional scRNA-seq datasets
#### [Dyngen simulation results](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Simulation_by_dyngen)

#### [HuBAMP datasets results](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/HuBMAP)
#### [Classic 10 datasets results](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Real_Datasets10)

### Accuracy in scRNA-seq with additional information
#### [Redeem Dadaset results](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Redeem)
#### [single-cell metabolic labeling dataset results (scNT-seq)](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/scNT)
#### [single-cell lineage tracing dataset (scLT)](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/scLT)

### Scalability and Stability
#### [Scalability](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Scalability)
#### [Stability](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Simulation_by_cellDancer)

### Overall performance
#### [Overall score](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Visualization/Vis_summary)

### Discussion
#### [Redeem Aged](https://github.com/liyarubio/RNA-velocity-benchmark/tree/main/Redeem/4_Aged)
