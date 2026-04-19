suppressPackageStartupMessages({
  library(rGREAT)
  library(rtracklayer)
  library(GenomeInfoDb)
})

# ===============================
# 0. Parse arguments
# ===============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript OCR_rGREAT.R <human|mouse>")
}

species <- tolower(args[1])

# ===============================
# 1. Settings
# ===============================
if (species == "human") {
  genome_build <- "hg38"
  input_file <- "/ocean/projects/bio230007p/jji5/data/Human_Pancreas_ATAC/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"
  outdir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results/Human_OCR_conservative"
  standard_chr <- paste0("chr", c(1:22, "X", "Y"))
  output_prefix <- "human_OCR"
} else if (species == "mouse") {
  genome_build <- "mm10"
  input_file <- "/ocean/projects/bio230007p/jji5/data/Mouse_Pancreas_ATAC/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"
  outdir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results/Mouse_OCR_conservative"
  standard_chr <- paste0("chr", c(1:19, "X", "Y"))
  output_prefix <- "mouse_OCR"
} else {
  stop("species must be either 'human' or 'mouse'")
}

# ===============================
# 2. Read input peaks
# ===============================
cat("Reading", species, "OCR peaks...\n")
flush.console()

gr <- import(input_file)

keep_chr <- intersect(seqlevels(gr), standard_chr)
gr <- keepSeqlevels(gr, keep_chr, pruning.mode = "coarse")

cat("Keeping seqlevels:\n")
print(keep_chr)

cat(species, "OCR peaks after filtering:", length(gr), "\n")
flush.console()

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ===============================
# 3. Run GREAT
# ===============================
cat("Running GREAT for", species, "OCRs...\n")
flush.console()

res <- great(
  gr,
  "GO:BP",
  paste0("GREAT:", genome_build)
)

tb <- getEnrichmentTable(res)

# ===============================
# 4. Save results
# ===============================
write.table(
  tb,
  file = file.path(outdir, paste0(output_prefix, "_GO_BP.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  res,
  file.path(outdir, paste0(output_prefix, "_great_object.rds"))
)

cat("Finished", species, "OCR GREAT analysis.\n")
flush.console()