suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath(getwd())
}

source_csv <- file.path(script_dir, "../../data/figure3_detail", "Figure3h_dataset10_global_performance.csv")
output_dir <- Sys.getenv("FIGURE3_OUT", file.path(script_dir, "../../../demo_output/figure3"))
output_table <- file.path(output_dir, "Figure3h_data10_global.csv")
output_pdf <- file.path(output_dir, "Figure3h_reproduced.pdf")
output_png <- file.path(output_dir, "Figure3h_reproduced.png")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_order <- c(
  "velocyto", "scVelo(sto)", "scVelo(dyn)", "VeloAE", "VeloVAE",
  "dynamo", "UniTVelo", "DeepVelo(VAE)", "DeepVelo(GCN)",
  "LatentVelo", "cellDancer", "veloVI", "STT", "TFvelo", "biVI",
  "DeepKINET", "TopicVelo", "SDEvelo", "cell2fate"
)

method_palette <- c(
  "velocyto" = "#E63946FF", "scVelo(sto)" = "#D8BFD8FF",
  "scVelo(dyn)" = "#8A2BE2FF", "VeloAE" = "#228B22FF",
  "VeloVAE" = "#3CB371FF", "dynamo" = "#FFD700FF",
  "UniTVelo" = "#F0E68CFF", "DeepVelo(VAE)" = "#0000CDFF",
  "DeepVelo(GCN)" = "#87CEFAFF", "LatentVelo" = "#20B2AAFF",
  "cellDancer" = "#FFB6C1FF", "veloVI" = "#DB7093FF",
  "STT" = "#FF8C00FF", "TFvelo" = "#556B2FFF",
  "biVI" = "#1E90FFFF", "DeepKINET" = "#FF6347FF",
  "TopicVelo" = "#9ACD32FF", "SDEvelo" = "#4682B4FF",
  "cell2fate" = "#CD5C5CFF"
)

dataset_columns <- c(
  "DentateGyrus (P0 and P5)", "DentateGyrus (P12 and P35)",
  "Pancreatic Endocrinogenesis", "CD34+ Bone Marrow (BM)",
  "Gastrulation", "Gastrulation Erythroid", "Forebrain", "Hindbrain",
  "Organogenesis", "Zebrafish"
)

normalize_method <- function(x) {
  key <- tolower(gsub("[ .-]+", "_", as.character(x)))
  out <- as.character(x)
  out[grepl("^velocyto$", key)] <- "velocyto"
  out[grepl("^scvelo_stochastic$|^scvelo$", key)] <- "scVelo(sto)"
  out[grepl("^scvelo_dynamical$|^scvelo_dynamic$", key)] <- "scVelo(dyn)"
  out[grepl("^veloae$", key)] <- "VeloAE"
  out[grepl("^velovae$", key)] <- "VeloVAE"
  out[grepl("^dynamo_stochastic$|^dynamo$", key)] <- "dynamo"
  out[grepl("^unitvelo$|^unittvelo$", key)] <- "UniTVelo"
  out[grepl("^deepvelo_sa$|^deepvelosa$", key)] <- "DeepVelo(VAE)"
  out[grepl("^deepvelo_gb$|^deepvelogb$", key)] <- "DeepVelo(GCN)"
  out[grepl("^latentvelo$", key)] <- "LatentVelo"
  out[grepl("^celldancer$", key)] <- "cellDancer"
  out[grepl("^velovi$", key)] <- "veloVI"
  out[grepl("^stt$", key)] <- "STT"
  out[grepl("^tfvelo$", key)] <- "TFvelo"
  out[grepl("^bivi$", key)] <- "biVI"
  out[grepl("^deepkinet$", key)] <- "DeepKINET"
  out[grepl("^topicvelo$", key)] <- "TopicVelo"
  out[grepl("^sdevelo$", key)] <- "SDEvelo"
  out[grepl("^cell2fate$", key)] <- "cell2fate"
  out
}

global_scores <- fread(source_csv, na.strings = c("", "NA", "NaN"))
global_scores[, Method := normalize_method(Method)]
global_scores <- global_scores[Method %in% method_order]

required_columns <- c("Method", dataset_columns, "Median", "Mean")
missing_columns <- setdiff(required_columns, names(global_scores))
if (length(missing_columns) > 0) {
  stop("Missing columns in global table: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(global_scores$Method)) {
  stop("Duplicated methods in global table.")
}

missing_methods <- setdiff(method_order, global_scores$Method)
if (length(missing_methods) > 0) {
  stop("Missing methods in global table: ", paste(missing_methods, collapse = ", "))
}

global_scores <- global_scores[match(method_order, Method)]
for (column in c(dataset_columns, "Median", "Mean")) {
  set(global_scores, j = column, value = suppressWarnings(as.numeric(global_scores[[column]])))
}

summary_table <- global_scores[, ..required_columns]
score_cols <- dataset_columns
calculated_median <- apply(summary_table[, ..score_cols], 1, median, na.rm = TRUE)
calculated_mean <- rowMeans(summary_table[, ..score_cols], na.rm = TRUE)
if (max(abs(calculated_median - summary_table$Median), na.rm = TRUE) > 1e-6 ||
    max(abs(calculated_mean - summary_table$Mean), na.rm = TRUE) > 1e-6) {
  stop("Median or Mean is inconsistent with the ten global metrics.")
}

if (nrow(summary_table) != 19L || length(score_cols) != 10L) {
  stop("Expected 19 methods and 10 classic datasets.")
}
if (sum(!is.na(summary_table[Method == "VeloAE", ..score_cols])) != 10L) {
  stop("VeloAE should have 10 valid classic global values.")
}
if (sum(!is.na(summary_table[Method == "dynamo", ..score_cols])) != 10L) {
  stop("dynamo should have 10 valid values in the updated table.")
}
if (sum(!is.na(summary_table[Method == "TopicVelo", ..score_cols])) != 8L) {
  stop("TopicVelo should have 8 valid classic global values in the updated table.")
}
# Order methods by available-case mean across the ten dataset columns.
# Keep original display order for ties, with wholly missing methods last.
mean_order <- order(-calculated_mean, seq_along(calculated_mean), na.last = TRUE)
summary_table <- summary_table[mean_order]
method_order <- summary_table$Method
fwrite(summary_table, output_table, na = "NA")

heatmap_matrix <- as.matrix(summary_table[, ..score_cols])
rownames(heatmap_matrix) <- summary_table$Method
heatmap_matrix <- t(heatmap_matrix)
heatmap_matrix <- cbind(
  heatmap_matrix,
  "Best performance" = apply(heatmap_matrix, 1, max, na.rm = TRUE)
)
# Mean uses the same fixed GDC scale and colors as the dataset cells.
# The reference envelope is not a method and has no Mean summary.
method_means <- calculated_mean[mean_order]
heatmap_matrix <- rbind(heatmap_matrix, Mean = c(method_means, NA_real_))
palette_table <- read.csv(file.path(script_dir, "../../../09_Core/visualization/figure3h_mean/assets/gdc_palette.csv"))

annotation_col <- data.frame(
  Method = factor(
    c(method_order, NA_character_),
    levels = method_order
  )
)
rownames(annotation_col) <- colnames(heatmap_matrix)

plot_heatmap <- function(filename, device) {
  if (device == "pdf") {
    pdf(filename, width = 10, height = 4, bg = "white")
  } else {
    png(filename, width = 3000, height = 1200, res = 300, bg = "white")
  }

  p <- pheatmap(
    heatmap_matrix,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = annotation_col,
    annotation_colors = list(Method = method_palette),
    color = palette_table$hex,
    breaks = seq(-1, 1, length.out = nrow(palette_table) + 1L),
    display_numbers = FALSE,
    gaps_row = length(dataset_columns),
    na_col = "white",
    silent = TRUE
  )
  grid.newpage()
  grid.draw(p$gtable)
  dev.off()
}

plot_heatmap(output_pdf, "pdf")
plot_heatmap(output_png, "png")

message("Wrote: ", output_table)
message("Wrote: ", output_pdf)
message("Wrote: ", output_png)
message(
  "TopicVelo classic-data missing values retained: ",
  paste(
    score_cols[as.logical(is.na(summary_table[Method == "TopicVelo", ..score_cols]))],
    collapse = ", "
  )
)
