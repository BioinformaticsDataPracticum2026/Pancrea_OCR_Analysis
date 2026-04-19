suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
})

# ===============================
# 0. Run mode
# ===============================
#run_mode <- "shared"   
run_mode <- "not_open" 


if (!run_mode %in% c("shared", "not_open")) {
  stop("run_mode must be either 'shared' or 'not_open'")
}

# ===============================
# 1. Paths
# ===============================
#genome_build <- "hg38" # human genome for shared
genome_build <- "mm10" # mouse genome for not_open

shared_file <- "/ocean/projects/bio230007p/jji5/output/bed/MtH/shared_ocrs.bed"
not_open_file <- "/ocean/projects/bio230007p/jji5/output/bed/MtH/mapped_not_open_in_human.bed"
#background_file <- "/ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak" #for shared, using human peaks as background since the input is shared peaks that are open in both species
background_file <- "/ocean/projects/bio230007p/jji5/data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.optimal_peak.narrowPeak" #for not_open, using mouse peaks as background since the input is mouse peaks that are not open in human 

base_outdir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results/MtH"

# ===============================
# 2. Helper function
# ===============================
run_great_analysis <- function(infile, bg_file, label, outdir, genome_build) {
  cat("Reading input for:", label, "\n")
  flush.console()

  gr <- import(infile)
  bg_gr <- import(bg_file)

  cat(label, "peaks:", length(gr), "\n")
  cat("Background peaks:", length(bg_gr), "\n")
  flush.console()

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  cat("Running GREAT for:", label, "\n")
  flush.console()

  res <- great(
    gr,
    "GO:BP",
    paste0("GREAT:", genome_build),
    background = bg_gr
  )

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

  cat("Finished:", label, "\n")
  flush.console()

  rm(gr, bg_gr, res, tb)
  gc()
}

# ===============================
# 3. Run by mode
# ===============================
if (run_mode == "shared") {
  run_great_analysis(
    infile = shared_file,
    bg_file = background_file,
    label = "shared",
    outdir = file.path(base_outdir, "shared"),
    genome_build = genome_build
  )
}

if (run_mode == "not_open") {
  run_great_analysis(
    infile = not_open_file,
    bg_file = background_file,
    label = "not_open",
    outdir = file.path(base_outdir, "not_open"),
    genome_build = genome_build
  )
}

cat("Done.\n")
flush.console()
