args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else normalizePath(getwd())
source(file.path(script_dir, "../shared/figure3_common.R"))
suppressPackageStartupMessages(library(corrplot))

metric_files <- c(
  Global = "Figure3f_global_metric.csv",
  Incluster = "Figure3f_incluster_metric.csv",
  Direction = "Figure3f_direction_metric.csv",
  velocity = "Figure3a_high_dim_cosine_similarity.csv",
  velocity_pca = "Figure3bc_low_dim_cosine_similarity.csv",
  trans_matrix = "Figure3f_transition_matrix_cosine_similarity.csv"
)

metric_tables <- lapply(names(metric_files), function(metric) {
  filename <- metric_files[[metric]]
  x <- fread(file.path(DATA_DIR, filename), check.names = FALSE)
  setnames(x, 1, "Method")
  x <- x[!grepl("Dynamo_deterministic", Method, ignore.case = TRUE)]
  x <- melt(
    x,
    id.vars = "Method",
    variable.name = "Dataset",
    value.name = metric
  )
  stopifnot(!anyDuplicated(x[, .(Method, Dataset)]))
  x
})
metric_table <- Reduce(
  function(x, y) merge(x, y, by = c("Method", "Dataset"), all = TRUE, sort = FALSE),
  metric_tables
)
metric_columns <- names(metric_files)
stopifnot(nrow(metric_table) == 798L)
fwrite(metric_table, file.path(OUTPUT_DIR, "Figure3f_aligned_metrics.csv"))

correlation <- cor(
  as.data.frame(metric_table[, ..metric_columns]),
  use = "pairwise.complete.obs"
)
fwrite(as.data.table(correlation, keep.rownames = "Metric"),
       file.path(OUTPUT_DIR, "Figure3f_correlation.csv"))
correlation_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)

draw_correlation <- function(filename, device) {
  if (device == "pdf") {
    pdf(filename, width = 9, height = 7, bg = "white")
  } else {
    png(filename, width = 2700, height = 2100, res = 300, bg = "white")
  }
  corrplot.mixed(
    correlation,
    lower = "number",
    upper = "circle",
    tl.pos = "d",
    diag = "n",
    tl.col = "red",
    tl.srt = 0,
    number.cex = 0.8,
    number.font = 2,
    upper.col = correlation_colors,
    lower.col = correlation_colors
  )
  dev.off()
}

draw_correlation(file.path(OUTPUT_DIR, "Figure3f_reproduced.pdf"), "pdf")
draw_correlation(file.path(OUTPUT_DIR, "Figure3f_reproduced.png"), "png")

expected <- c(Global_Incluster = 0.17, Global_Direction = 0.77,
              velocity_trans_matrix = 0.73)
observed <- c(correlation["Global", "Incluster"],
              correlation["Global", "Direction"],
              correlation["velocity", "trans_matrix"])
stopifnot(isTRUE(all.equal(unname(round(observed, 2)), unname(expected))))
message("Wrote Figure3f_reproduced; key correlations validated")
