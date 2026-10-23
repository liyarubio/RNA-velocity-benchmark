args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else normalizePath(getwd())
source(file.path(script_dir, "figure3_common.R"))

pca_file <- file.path(DATA_DIR, "Figure3d_pca2_to_30_cosine_similarity.csv")
high_file <- file.path(DATA_DIR, "Figure3a_high_dim_cosine_similarity.csv")
pca <- fread(pca_file, check.names = FALSE)
high <- fread(high_file, check.names = FALSE)
setnames(high, 1L, "Method")

cosine_vector <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!length(x) || denominator == 0) return(NA_real_)
  sum(x * y) / denominator
}

sheet_names <- unique(pca$Dataset)
stopifnot(length(sheet_names) == 42L)
dataset_results <- vector("list", length(sheet_names))
all_low <- vector("list", 29L)
all_high <- vector("list", 29L)

for (sheet_index in seq_along(sheet_names)) {
  dataset <- sheet_names[[sheet_index]]
  low <- as.data.frame(
    pca[pca$Dataset == dataset, setdiff(names(pca), "Dataset"), with = FALSE]
  )
  stopifnot(nrow(low) == 29L)
  methods <- names(low)[-1]
  methods <- methods[!grepl("Dynamo_deterministic", methods, ignore.case = TRUE)]
  high_values <- high[[dataset]][match(methods, high$Method)]

  similarities <- vapply(seq_len(nrow(low)), function(row_index) {
    low_values <- as.numeric(low[row_index, methods, drop = TRUE])
    all_low[[row_index]] <<- c(all_low[[row_index]], low_values)
    all_high[[row_index]] <<- c(all_high[[row_index]], high_values)
    cosine_vector(low_values, high_values)
  }, numeric(1))

  dataset_results[[sheet_index]] <- data.table(
    Dataset = dataset,
    Dataset_type = sub("_seed[0-9]+$", "", dataset),
    Dim = as.integer(sub("^dim_", "", low$Dim)),
    Cosine_similarity = similarities
  )
}

dataset_data <- rbindlist(dataset_results)
plot_data <- dataset_data[, .(
  Cosine_similarity = mean(Cosine_similarity, na.rm = TRUE)
), by = .(Dataset_type, Dim)]

total_data <- data.table(
  Dataset_type = "total",
  Dim = 2:30,
  Cosine_similarity = vapply(1:29, function(i) {
    cosine_vector(all_low[[i]], all_high[[i]])
  }, numeric(1))
)
plot_data <- rbind(plot_data, total_data)
stopifnot(nrow(plot_data) == 15L * 29L)

dataset_order <- c(
  "bifurcating_converging", "bifurcating_cycle", "bifurcating_loop",
  "bifurcating", "binary_tree", "branching", "consecutive_bifurcating",
  "converging", "cycle", "cycle_simple", "disconnected", "linear",
  "linear_simple", "trifurcating", "total"
)
plot_data[, Dataset_type := factor(Dataset_type, levels = dataset_order)]

dataset_colors <- c(
  "#F8766D", "#E78600", "#C49A00", "#A3A500", "#39B600",
  "#00BA38", "#00C08B", "#00BFC4", "#00B0F6", "#00A6FF",
  "#9590FF", "#E76BF3", "#F564E3", "#FF62BC", "#000000"
)
names(dataset_colors) <- dataset_order

panel_d <- ggplot(plot_data, aes(Dim, Cosine_similarity, colour = Dataset_type)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = dataset_colors, drop = FALSE) +
  scale_x_continuous(breaks = 2:30) +
  labs(title = "Cosine similarity", x = "Dim", colour = "Dataset") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 7)
  )

save_ggplot(panel_d, "Figure3d_reproduced", width = 10, height = 6)
fwrite(plot_data, file.path(OUTPUT_DIR, "Figure3d_plot_data.csv"))
message("Wrote Figure3d_reproduced")
