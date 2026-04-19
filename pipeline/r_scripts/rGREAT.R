suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
  library(yaml)
})

load_config <- function(config_path = "config.yaml") {
  yaml::read_yaml(config_path)
}

validate_run_mode <- function(run_mode) {
  if (!run_mode %in% c("shared", "not_open")) {
    stop("run_mode must be either 'shared' or 'not_open'")
  }
}

get_comparison_config <- function(cfg, comparison) {
  comp_cfg <- cfg$rgreat$comparisons[[comparison]]

  if (is.null(comp_cfg)) {
    stop(paste("comparison", comparison, "not found in config.yaml"))
  }

  comp_cfg
}

get_input_file <- function(comp_cfg, run_mode) {
  switch(
    run_mode,
    shared = comp_cfg$shared_file,
    not_open = comp_cfg$not_open_file
  )
}

read_peak_sets <- function(infile, bg_file) {
  gr <- import(infile)
  bg_gr <- import(bg_file)
  list(gr = gr, bg_gr = bg_gr)
}

run_great <- function(gr, bg_gr, genome_build) {
  great(
    gr,
    "GO:BP",
    paste0("GREAT:", genome_build),
    background = bg_gr
  )
}

save_great_results <- function(res, label, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tb <- getEnrichmentTable(res)

  write.table(
    tb,
    file = file.path(outdir, paste0(label, "_GO_BP.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(
    res,
    file.path(outdir, paste0(label, "_great_object.rds"))
  )
}

run_great_analysis <- function(infile, bg_file, label, outdir, genome_build) {
  cat("Reading input for:", label, "\n")
  flush.console()

  peak_sets <- read_peak_sets(infile, bg_file)
  gr <- peak_sets$gr
  bg_gr <- peak_sets$bg_gr

  cat(label, "peaks:", length(gr), "\n")
  cat("Background peaks:", length(bg_gr), "\n")
  flush.console()

  cat("Running GREAT for:", label, "\n")
  flush.console()

  res <- run_great(gr, bg_gr, genome_build)
  save_great_results(res, label, outdir)

  cat("Finished:", label, "\n")
  flush.console()

  rm(gr, bg_gr, res, peak_sets)
  gc()
}

run_one_comparison_mode <- function(cfg, comparison, run_mode) {
  validate_run_mode(run_mode)

  comp_cfg <- get_comparison_config(cfg, comparison)
  infile <- get_input_file(comp_cfg, run_mode)
  outdir <- file.path(comp_cfg$base_outdir, run_mode)

  cat("====================================\n")
  cat("Comparison:", comparison, "\n")
  cat("Run mode:", run_mode, "\n")
  cat("Genome build:", comp_cfg$genome_build, "\n")
  cat("Input file:", infile, "\n")
  cat("Background file:", comp_cfg$background_file, "\n")
  cat("Output dir:", outdir, "\n")
  cat("====================================\n")
  flush.console()

  run_great_analysis(
    infile = infile,
    bg_file = comp_cfg$background_file,
    label = run_mode,
    outdir = outdir,
    genome_build = comp_cfg$genome_build
  )
}

run_all_comparison_modes <- function(config_path = "config.yaml", comparisons = NULL) {
  cfg <- load_config(config_path)

  if (is.null(comparisons)) {
    comparisons <- names(cfg$rgreat$comparisons)
  }

  run_modes <- c("shared", "not_open")

  for (comparison in comparisons) {
    for (run_mode in run_modes) {
      run_one_comparison_mode(cfg, comparison, run_mode)
    }
  }

  cat("Done.\n")
  flush.console()
}

# ===============================
# Call at the end
# ===============================
run_all_comparison_modes(
  config_path = "config.yaml",
  comparisons = c("HtM", "MtH")
)