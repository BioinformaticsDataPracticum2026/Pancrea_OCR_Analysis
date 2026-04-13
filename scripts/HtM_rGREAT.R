suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
})

# ===============================
# 0. Run mode
# ===============================
run_mode <- "shared"   # "shared" or "not_open"

if (!run_mode %in% c("shared", "not_open")) {
  stop("run_mode must be either 'shared' or 'not_open'")
}

# ===============================
# 1. Paths
# ===============================
genome_build <- "mm10" # mouse genome

shared_file <- "/ocean/projects/bio230007p/jji5/output/bed/HtM/shared_ocrs.bed"
not_open_file <- "/ocean/projects/bio230007p/jji5/output/bed/HtM/mapped_not_open_in_mouse.bed"
background_file <- "/ocean/projects/bio230007p/jji5/output/bed/HtM/mapped_mouse_halper.bed"

base_outdir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results/HtM"

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
