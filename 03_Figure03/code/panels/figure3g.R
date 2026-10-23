args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else normalizePath(getwd())
source(file.path(script_dir, "figure3_common.R"))

input_csv <- file.path(DATA_DIR, "Figure3g_hubmap_global_performance.csv")
wide <- fread(input_csv, check.names = FALSE, na.strings = c("", "NA", "NaN"))
id_column <- names(wide)[1]
method_columns <- setdiff(names(wide), c(id_column, "Dynamo_deterministic"))

long <- melt(
  wide,
  id.vars = id_column,
  measure.vars = method_columns,
  variable.name = "Raw_method",
  value.name = "Global"
)
long[, Method := normalize_method(Raw_method)]
long <- long[Method %in% method_order]

method_medians <- long[, .(Median = median(Global, na.rm = TRUE)), by = Method]
setorder(method_medians, -Median)
plot_order <- method_medians$Method
stopifnot(uniqueN(long[[id_column]]) == 176L, length(plot_order) == 18L)
long[, Method := factor(Method, levels = plot_order)]

set.seed(2024) # Reproducible jitter; numerical summaries are unchanged.
panel_g <- ggplot(long, aes(Method, Global, fill = Method)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.35, colour = "black") +
  scale_fill_manual(values = method_palette, guide = "none") +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  labs(
    title = "Global Performance Comparison (Sorted by Median)",
    x = NULL,
    y = "Global Metric Score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_ggplot(panel_g, "Figure3g_reproduced", width = 14, height = 7)
fwrite(long, file.path(OUTPUT_DIR, "Figure3g_long_data.csv"), na = "NA")
fwrite(method_medians, file.path(OUTPUT_DIR, "Figure3g_method_medians.csv"))
message("Wrote Figure3g_reproduced")
