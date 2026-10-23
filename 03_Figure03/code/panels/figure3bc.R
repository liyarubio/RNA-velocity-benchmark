args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else normalizePath(getwd())
source(file.path(script_dir, "figure3_common.R"))

high <- read_metric_matrix(
  file.path(DATA_DIR, "Figure3a_high_dim_cosine_similarity.csv"),
  "High"
)
low <- read_metric_matrix(
  file.path(DATA_DIR, "Figure3bc_low_dim_cosine_similarity.csv"),
  "Low"
)

paired_all <- merge(high, low, by = c("Method", "Dataset"), all = TRUE)

# Panel b used check_method_name() in the original script, which removed the
# deterministic Dynamo run before method-level summarization.
paired <- paired_all[!grepl("Dynamo_deterministic", Method, ignore.case = TRUE)]
paired[, Method := factor(Method, levels = method_order)]

summary_data <- paired[, .(
  HighMean = mean(High, na.rm = TRUE),
  HighSD = sd(High, na.rm = TRUE),
  LowMean = mean(Low, na.rm = TRUE),
  LowSD = sd(Low, na.rm = TRUE)
), by = Method][Method != "VeloAE" & !is.na(Method)]

panel_b <- ggplot(summary_data, aes(HighMean, LowMean, colour = Method)) +
  geom_errorbar(aes(ymin = LowMean - LowSD, ymax = LowMean + LowSD)) +
  ggplot2::geom_errorbarh(aes(xmin = HighMean - HighSD, xmax = HighMean + HighSD)) +
  geom_point(size = 3) +
  scale_colour_manual(values = method_palette, drop = TRUE) +
  labs(x = "High-dim Cosine Similarity", y = "Low-dim Cosine Similarity") +
  guides(colour = guide_legend(ncol = 2, byrow = FALSE)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom", legend.title = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black")
  )

save_ggplot(panel_b, "Figure3b_reproduced", width = 7, height = 8)

# Panel c was built before check_method_name() in the original script. Keep the
# deterministic Dynamo row here to preserve the published violin distribution.
panel_c_pairs <- copy(paired_all)
panel_c_pairs[, Pair := paste(Method, Dataset, sep = ".")]
panel_c_data <- melt(
  panel_c_pairs[, .(Pair, High, Low)], id.vars = "Pair",
  variable.name = "Dim", value.name = "Cosine_sim"
)

panel_c <- ggplot(panel_c_data, aes(Dim, Cosine_sim, fill = Dim)) +
  geom_violin(linewidth = 1) +
  geom_line(aes(group = Pair), colour = "grey80", linewidth = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(High = "#08306B", Low = "#C6DBEF")) +
  labs(x = "Dim", y = "Cosine_sim") +
  ggpubr::theme_classic2(base_size = 20) +
  theme(legend.position = "none")

save_ggplot(panel_c, "Figure3c_reproduced", width = 4, height = 6)
fwrite(paired_all, file.path(OUTPUT_DIR, "Figure3bc_paired_data.csv"), na = "NA")

message("Wrote Figure3b_reproduced and Figure3c_reproduced")
