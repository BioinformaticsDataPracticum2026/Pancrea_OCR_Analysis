library(rGREAT)
base <- "/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR"
outdir <- file.path(base, "rGREAT_Analysis/outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load_gr <- function(f) {
  d <- read.table(f, sep = "\t", header = FALSE)
  GRanges(seqnames = d[,1], ranges = IRanges(d[,2], d[,3]))
}

inputs <- list(
  list("human_all",      file.path(base, "Mapping/outputs/human_liver.narrowPeak"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_all",      file.path(base, "Mapping/outputs/mouse_liver.narrowPeak"), "TxDb.Mmusculus.UCSC.mm10.knownGene"),
  list("shared",         file.path(base, "PE_classification/output/unique/shared_open.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("human_specific", file.path(base, "PE_classification/output/unique/human_open_mouse_closed.bed"), "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  list("mouse_specific", file.path(base, "PE_classification/output/unique/mouse_open_human_closed.bed"), "TxDb.Mmusculus.UCSC.mm10.knownGene")
)

for (inp in inputs) {
  nm <- inp[[1]]; fp <- inp[[2]]; gn <- inp[[3]]
  cat("Running:", nm, "\n")
  gr <- load_gr(fp)
  job <- great(gr, gene_sets = "GO:BP", tss_source = gn)
  res <- getEnrichmentTable(job)
  sig <- res[res$p_adjust < 0.05, ]
  cat(" ", nrow(sig), "significant terms\n")
  write.csv(res, file.path(outdir, paste0(nm, "_GOBP.csv")), row.names = FALSE)
  write.csv(sig, file.path(outdir, paste0(nm, "_GOBP_sig.csv")), row.names = FALSE)
}
