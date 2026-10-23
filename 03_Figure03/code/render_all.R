args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else normalizePath(getwd())

scripts <- c(
  "a/figure3a.R",
  "b/figure3b.R",
  "c/figure3c.R",
  "d/figure3d.R",
  "f/figure3f.R",
  "g/figure3g.R",
  "h/figure3h.R"
)

for (script in scripts) {
  message("\n===== Running ", script, " =====")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(script_dir, script)))
  )
  if (!identical(status, 0L)) stop(script, " failed with status ", status)
}

message("\nSaved-data redraws for panels a-d and f-h written to: ", Sys.getenv("FIGURE3_OUT", file.path(script_dir, "../../demo_output/figure3")))
message("Panel e is retained manual artwork in the formal Figure03; no numerical generator applies.")
