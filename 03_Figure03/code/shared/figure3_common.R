suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

figure3_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  normalizePath(getwd())
}

FIGURE3_DIR <- figure3_script_dir()
OUTPUT_DIR <- Sys.getenv("FIGURE3_OUT", file.path(FIGURE3_DIR, "../../../demo_output/figure3"))
DATA_DIR <- normalizePath(file.path(FIGURE3_DIR, "../../data/figure3_detail"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DATA_DIR)) stop("Local data directory not found: ", DATA_DIR)

method_order <- c(
  "velocyto", "scVelo(sto)", "scVelo(dyn)", "VeloAE", "VeloVAE",
  "dynamo", "UniTVelo", "DeepVelo(VAE)", "DeepVelo(GCN)",
  "LatentVelo", "cellDancer", "veloVI", "STT", "TFvelo", "biVI",
  "DeepKINET", "TopicVelo", "SDEvelo", "cell2fate"
)

method_palette <- c(
  "velocyto" = "#D62728", "scVelo(sto)" = "#CAB2D6",
  "scVelo(dyn)" = "#6A3D9A", "VeloAE" = "#33A02C",
  "VeloVAE" = "#B2DF8A", "dynamo" = "#FFB547",
  "UniTVelo" = "#F0E68C", "DeepVelo(VAE)" = "#0000CD",
  "DeepVelo(GCN)" = "#87CEFA", "LatentVelo" = "#20B2AA",
  "cellDancer" = "#FFB6C1", "veloVI" = "#DB7093",
  "STT" = "#FF8C00", "TFvelo" = "#7B4173", "biVI" = "#1E90FF",
  "DeepKINET" = "#FF6347", "TopicVelo" = "#9ACD32",
  "SDEvelo" = "#4682B4", "cell2fate" = "#CD5C5C"
)

normalize_method <- function(x) {
  key <- tolower(gsub("[ .()-]+", "_", as.character(x)))
  key <- gsub("_+", "_", key)
  key <- gsub("^_|_$", "", key)
  out <- as.character(x)
  out[grepl("velocyto$", key)] <- "velocyto"
  out[grepl("scvelo_stochastic$|^scvelo$", key)] <- "scVelo(sto)"
  out[grepl("scvelo_dynamical$|scvelo_dynamic$", key)] <- "scVelo(dyn)"
  out[grepl("veloae$", key)] <- "VeloAE"
  out[grepl("velovae$", key)] <- "VeloVAE"
  out[grepl("dynamo_stochastic$|^dynamo$", key)] <- "dynamo"
  out[grepl("unitvelo$|unittvelo$|unitevelo$", key)] <- "UniTVelo"
  out[grepl("deepvelo_sa$|deepvelosa$", key)] <- "DeepVelo(VAE)"
  out[grepl("deepvelo_gb$|deepvelogb$", key)] <- "DeepVelo(GCN)"
  out[grepl("latentvelo$", key)] <- "LatentVelo"
  out[grepl("celldancer$", key)] <- "cellDancer"
  out[grepl("velovi$", key)] <- "veloVI"
  out[grepl("stt$", key)] <- "STT"
  out[grepl("tfvelo$", key)] <- "TFvelo"
  out[grepl("bivi$", key)] <- "biVI"
  out[grepl("deepkinet$", key)] <- "DeepKINET"
  out[grepl("topicvelo$", key)] <- "TopicVelo"
  out[grepl("sdevelo$", key)] <- "SDEvelo"
  out[grepl("cell2fate$", key)] <- "cell2fate"
  out
}

read_metric_matrix <- function(path, value_name) {
  x <- fread(path, check.names = FALSE, na.strings = c("", "NA", "NaN"))
  setnames(x, 1L, "Method")
  x[, Method := sub("^[0-9]+_", "", Method)]
  long <- melt(x, id.vars = "Method", variable.name = "Dataset", value.name = value_name)
  long[, Method := normalize_method(Method)]
  long[]
}

save_ggplot <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in")
  ggsave(file.path(OUTPUT_DIR, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
}
