setwd("/media/zx/HDD1/cooperation/liyr/analysis/dyngen")
qc_pass_names <- c(
    "zenodo_1443566_real_gold_cellbench-SC1_luyitian",
    "zenodo_1443566_real_gold_cellbench-SC2_luyitian",
    "zenodo_1443566_real_gold_cellbench-SC3_luyitian",
    "zenodo_1443566_real_gold_developing-dendritic-cells_schlitzer",
    "zenodo_1443566_real_gold_psc-astrocyte-maturation-neuron_sloan",
    "zenodo_1443566_real_silver_bone-marrow-mesenchyme-erythrocyte-differentiation_mca",
    "zenodo_1443566_real_silver_cell-cycle_leng",
    "zenodo_1443566_real_silver_embronic-mesenchyme-neuron-differentiation_mca",
    "zenodo_1443566_real_silver_embryonic-mesenchyme-stromal-cell-cxcl14-cxcl12-axis_mca",
    "zenodo_1443566_real_silver_mammary-gland-involution-endothelial-cell-aqp1-gradient_mca",
    "zenodo_1443566_real_silver_mouse-cell-atlas-combination-10",
    "zenodo_1443566_real_silver_mouse-cell-atlas-combination-4",
    "zenodo_1443566_real_silver_mouse-cell-atlas-combination-5",
    "zenodo_1443566_real_silver_mouse-cell-atlas-combination-8",
    "zenodo_1443566_real_silver_olfactory-projection-neurons-DA1_horns",
    "zenodo_1443566_real_silver_placenta-trophoblast-differentiation_mca",
    "zenodo_1443566_real_silver_thymus-t-cell-differentiation_mca",
    "zenodo_1443566_real_silver_trophoblast-stem-cell-trophoblast-differentiation_mca"
)
for (i in qc_pass_names){
system(paste0("wget https://github.com/dynverse/dyngen/raw/data_files/",i,".rds"))
}
