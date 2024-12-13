setwd("/media/zx/HDD1/cooperation/liyr/result/dyngen/1_simulation_datasets")
library(tidyverse)
library(anndata)
library(rlang)
library(dyngen)
#library(dyngen.manuscript)
#library(dynplot2)
library(dyndimred)

###### mkdir
make_directory_function <- function(prefix, postfix = character(0)) {
    function(...) {
        file <- do.call(file.path, as.list(c(prefix, paste0(...), postfix)))
        folder <- gsub("[^/]*$", "", file)
        if (!file.exists(folder)) {
            dir.create(folder, recursive = TRUE)
        }
        file
    }
}

start_analysis <- function(experiment_id) {
    list(
        temporary = make_directory_function(paste0("temporary_files/", experiment_id)),
        result = make_directory_function(paste0("result_files/", experiment_id)),
        dataset_folder = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = ""),
        model_file = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = "model.rds"),
        dataset_file = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = "dataset.rds"),
        velocity_file = function(dataset_id, method_id, params_id) {
            make_directory_function(
                prefix = paste0("temporary_files/", experiment_id, "/velocity"),
                postfix = "velocity.rds"
            )(
                paste0(dataset_id, "-", method_id, "-", params_id)
            )
        }
    )
}

exp <- start_analysis("usecase_rna_velocity")

# file.remove(exp$result("design_datasets.rds"))

############ setup dataset design
design_datasets <-
    crossing(
        seed = 1:3,
        backbone_name = names(list_backbones())) %>%
    mutate(
        id = paste0(backbone_name, "_seed", seed)
    )
write_rds(design_datasets, exp$result("design_datasets.rds"), compress = "gz")

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

# backbones <- dyngen::list_backbones()
# backbones$disconnected <- function() dyngen::backbone_disconnected("linear", "cycle")

################# build datasets

pwalk(design_datasets, function(id, seed, backbone_name, ...) {
    dataset_file <- exp$dataset_file(id)
    plot_file <- gsub("\\.rds$", ".pdf", dataset_file)
    
    if (!file.exists(plot_file)) {
        
        cat("## Generating ", id, "\n", sep = "")
        set.seed(seed)
        
        backbone <- backbones[[backbone_name]]()
        model <-
            initialise_model(
                id = id,
                num_tfs = nrow(backbone$module_info),
                num_targets = 70,
                num_hks = 0,
                backbone = backbone,
                num_cells = 2500,
                simulation_params = simulation_default(
                    burn_time = simtime_from_backbone(backbone, burn = TRUE) * 1.5,
                    census_interval = ifelse(!backbone_name %in% c("linear_simple"), 10, 1),
                    compute_rna_velocity = TRUE,
                    store_reaction_propensities = TRUE,
                    experiment_params = simulation_type_wild_type(
                        num_simulations = 100
                    )
                ),
                verbose = FALSE,
                download_cache_dir = "/media/zx/HDD1/cooperation/liyr/analysis/dyngen/download" ### 这里最好提前下载好一些数据
            )
        generate_dataset(
            model,
            format = "dyno",
            output_dir = exp$dataset_folder(id),
            make_plots = TRUE
        )
        
        
        # add dimred
        dataset <- read_rds(dataset_file)
        
        dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson", ndim = 3)
        
        dataset <- dataset %>%
            dynwrap::add_dimred(dimred = dimred) %>%
            dynwrap::add_waypoints()
        # pl <- dynwrap::plot_dimred(dataset) +
        #     geom_cell_point(aes(color = milestone_percentages), size = 1) +
        #     scale_milestones_colour() +
        #     geom_trajectory_segments(size = 1, color = "#333333") +
        #     geom_milestone_label(aes(label = label), color = "black", fill = "#EEEEEE") +
        #     theme_common(legend.position = "none")
        # ggsave(plot_file, pl, width = 6, height = 5)
        write_rds(dataset, dataset_file, compress = "gz")
        
        gc()
    }
})

################################ transfer to h5ad

data_dir <- "/media/zx/HDD1/cooperation/liyr/result/dyngen/1_simulation_datasets/temporary_files/usecase_rna_velocity/datasets/"
file_list <- list.files(path = data_dir)

for (i in file_list){
    file_path <- paste0(data_dir,i,"/dataset.rds")
    data_use <- readRDS(file = file_path)
    data_use$cell_info$clusters = paste0(as.character(data_use$progressions[["from"]]),"_",as.character(data_use$progressions[["to"]]))
    
    cell_info = data.frame(data_use$cell_info)
    rownames(cell_info) <- data_use$cell_info$cell_id
    feature_info = data.frame(data_use$feature_info)
    rownames(feature_info) <- data_use$feature_info$feature_id
    
    adata <- AnnData(
        X = data_use$counts,
        obs = cell_info,
        var = feature_info,
        layers = list(
            spliced = data_use$counts_spliced,
            unspliced = data_use$counts_unspliced,
            ground_truth_velocity = data_use$rna_velocity
        ),
        obsm = list(
            dimred = data_use$dimred
        ),
        varm = list(),
        uns = list(
            traj_dimred_segments = data_use$dimred_segment_progressions
        )
    )
    write_h5ad(anndata = adata,filename = paste0(data_dir,i,"/anndata.h5ad"))
}