# run with:
# Rscript pipeline/r_scripts/rGREAT.R Human_full
# Rscript pipeline/r_scripts/rGREAT.R Mouse_full
# Rscript pipeline/r_scripts/rGREAT.R HtM_shared
# Rscript pipeline/r_scripts/rGREAT.R HtM_human_specific
# Rscript pipeline/r_scripts/rGREAT.R MtH_shared
# Rscript pipeline/r_scripts/rGREAT.R MtH_mouse_specific

user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
  library(yaml)
})

resolve_path <- function(project_root, path) {
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(project_root, path)
}

load_config <- function(config_path = "config.yaml") {
  if (!file.exists(config_path)) {
    stop("config.yaml not found. Run this script from the project root.")
  }

  cfg <- yaml::read_yaml(config_path)

  if (is.null(cfg$project_root)) {
    cfg$project_root <- getwd()
  }

  cfg$project_root <- normalizePath(cfg$project_root, mustWork = TRUE)

  path_keys <- c(
    "species_1_peak_file",
    "species_2_peak_file",
    "bed_output_dir_htm",
    "bed_output_dir_mth",
    "rgreat_output_dir",
    "rgreat_temp_dir"
  )

  for (key in path_keys) {
    if (is.null(cfg[[key]])) {
      stop(paste("Missing config key:", key))
    }
    cfg[[key]] <- resolve_path(cfg$project_root, cfg[[key]])
  }

  dir.create(cfg$rgreat_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$rgreat_temp_dir, recursive = TRUE, showWarnings = FALSE)

  cfg
}

get_task_info <- function(cfg, task) {
  htm_shared <- file.path(cfg$bed_output_dir_htm, "shared_ocrs.bed")
  htm_human_specific <- file.path(cfg$bed_output_dir_htm, "human_specific_ocrs.bed")

  mth_shared <- file.path(cfg$bed_output_dir_mth, "shared_ocrs.bed")
  mth_mouse_specific <- file.path(cfg$bed_output_dir_mth, "mouse_specific_ocrs.bed")

  if (task == "Human_full") {
    return(list(
      infile = cfg$species_1_peak_file,
      outdir = file.path(cfg$rgreat_output_dir, "Human_full"),
      label = "human_full",
      genome = "hg38"
    ))
  }

  if (task == "Mouse_full") {
    return(list(
      infile = cfg$species_2_peak_file,
      outdir = file.path(cfg$rgreat_output_dir, "Mouse_full"),
      label = "mouse_full",
      genome = "mm10"
    ))
  }

  if (task == "HtM_shared") {
    return(list(
      infile = htm_shared,
      outdir = file.path(cfg$rgreat_output_dir, "HtM_shared"),
      label = "htm_shared",
      genome = "mm10"
    ))
  }

  if (task == "HtM_human_specific") {
    return(list(
      infile = htm_human_specific,
      outdir = file.path(cfg$rgreat_output_dir, "HtM_human_specific"),
      label = "htm_human_specific",
      genome = "mm10"
    ))
  }

  if (task == "MtH_shared") {
    return(list(
      infile = mth_shared,
      outdir = file.path(cfg$rgreat_output_dir, "MtH_shared"),
      label = "mth_shared",
      genome = "hg38"
    ))
  }

  if (task == "MtH_mouse_specific") {
    return(list(
      infile = mth_mouse_specific,
      outdir = file.path(cfg$rgreat_output_dir, "MtH_mouse_specific"),
      label = "mth_mouse_specific",
      genome = "hg38"
    ))
  }

  stop(paste("Unknown task:", task))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste(
      "Missing task.",
      "Valid tasks:",
      "Human_full, Mouse_full, HtM_shared, HtM_human_specific, MtH_shared, MtH_mouse_specific"
    )
  )
}

task <- args[1]

cfg <- load_config()
info <- get_task_info(cfg, task)

if (!file.exists(info$infile)) {
  stop(paste("Input file not found:", info$infile))
}

dir.create(info$outdir, recursive = TRUE, showWarnings = FALSE)

message("Running rGREAT")
message("Task: ", task)
message("Input file: ", info$infile)
message("Genome: ", info$genome)
message("Output directory: ", info$outdir)

gr <- rtracklayer::import(info$infile)

job <- submitGreatJob(
  gr,
  genome = info$genome,
  help = FALSE
)

tbl <- getEnrichmentTables(
  job,
  ontology = "GO Biological Process"
)

outfile <- file.path(info$outdir, paste0(info$label, "_GO_BP.tsv"))

write.table(
  tbl,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Finished rGREAT")
message("Output file: ", outfile)