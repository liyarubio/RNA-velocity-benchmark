args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath(getwd())
}

input_csv <- file.path(script_dir, "../../data/figure3_detail", "Figure3a_high_dim_cosine_similarity.csv")
output_dir <- Sys.getenv("FIGURE3_OUT", file.path(script_dir, "../../../demo_output/figure3"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# This local file contains the complete original plot_summary() implementation.
source(file.path(script_dir, "figure3a_plot_core.R"))

plot_summary(
  input_csv,
  file.path(output_dir, "Figure3a_reproduced.pdf"),
  title_text = "Cosine similarity of high-dim RNA velocity on cell level",
  panel_label = "a"
)

message("Wrote Figure3a_reproduced.pdf")
