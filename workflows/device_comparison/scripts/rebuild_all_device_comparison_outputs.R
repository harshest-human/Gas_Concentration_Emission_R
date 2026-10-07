script_dir <- {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args_all[grepl(file_arg, args_all)])
  if (length(script_path) > 0) dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_bin)) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
}

run_step <- function(script_name, extra_args = character()) {
  script_path <- file.path(script_dir, script_name)
  cat("RUNNING:", script_path, "\n")
  status <- system2(rscript_bin, c(script_path, extra_args))
  if (!identical(status, 0L)) {
    stop("Script failed: ", script_name, call. = FALSE)
  }
}

run_step("build_pronova_campaign_hourly.R")
run_step("build_cubic_campaign_hourly.R")
run_step("build_otice_campaign_hourly.R")
run_step("build_analyzer_comparison_datasets.R")
run_step("plot_campaign_monthly_device_comparison.R")
run_step("summarize_cigr_conference_metrics.R")
run_step("run_device_comparison_period.R")
run_step("plot_pronova_vs_crds_2026.R")
run_step("summarize_dec_2025_mean_sd.R")

cat("Finished rebuilding device comparison outputs.\n")
