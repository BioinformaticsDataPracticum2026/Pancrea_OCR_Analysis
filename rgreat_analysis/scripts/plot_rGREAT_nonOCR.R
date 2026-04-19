library(ggplot2)

# ===============================
# 0. Settings
# ===============================
comparison <- "HtM"     # choose from: "HtM", "MtH"
run_mode <- "shared"    # choose from: "shared", "not_open"
top_n <- 10

# ===============================
# 1. Paths
# ===============================
base_dir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results"

input_file <- file.path(
  base_dir, comparison, run_mode,
  paste0(run_mode, "_GO_BP.tsv")
)

output_file <- file.path(
  base_dir, comparison, "plots",
  paste0(run_mode, "_GO_BP_top", top_n, ".pdf")
)

plot_title <- paste(comparison, "-", run_mode, "OCRs")

cat("Comparison:", comparison, "\n")
cat("Run mode:", run_mode, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")

# ===============================
# 2. Read table
# ===============================
tb <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# ===============================
# 3. Clean data
# ===============================
tb$p_adjust[tb$p_adjust == 0] <- 1e-300
tb <- tb[!is.na(tb$p_adjust), , drop = FALSE]

if (nrow(tb) == 0) {
  stop("No valid rows found after filtering p_adjust.")
}

tb$neg_log10_padj <- -log10(tb$p_adjust)

# ===============================
# 4. Select top terms
# ===============================
tb_top <- tb[order(tb$p_adjust, decreasing = FALSE), , drop = FALSE]
tb_top <- head(tb_top, top_n)

if (!"description" %in% colnames(tb_top)) {
  stop("Column 'description' not found in enrichment table.")
}

if (!"fold_enrichment" %in% colnames(tb_top)) {
  stop("Column 'fold_enrichment' not found in enrichment table.")
}

tb_top$description <- factor(tb_top$description, levels = rev(tb_top$description))

# ===============================
# 5. Plot
# ===============================
p <- ggplot(tb_top, aes(x = neg_log10_padj, y = description)) +
  geom_point(aes(size = fold_enrichment, color = neg_log10_padj), alpha = 0.85) +
  labs(
    title = plot_title,
    x = expression(-log[10]("adjusted p-value")),
    y = NULL,
    size = "Fold enrichment",
    color = expression(-log[10]("adj. p"))
  ) +
  theme_bw(base_size = 12)

ggsave(
  output_file,
  plot = p,
  width = 8,
  height = 5
)

cat("Plot done.\n")
