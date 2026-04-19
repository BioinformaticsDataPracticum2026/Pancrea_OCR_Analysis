library(ggplot2)
library(stringr)
library(grid)

# ===============================
# 0. Settings
# ===============================
top_n <- 10
base_dir <- "/ocean/projects/bio230007p/kfang5/rgreat_analysis/results"

# ===============================
# 1. Paths
# ===============================
human_file <- file.path(base_dir, "Human_OCR_conservative", "human_OCR_GO_BP.tsv")
mouse_file <- file.path(base_dir, "Mouse_OCR_conservative", "mouse_OCR_GO_BP.tsv")

output_file <- file.path(
  base_dir,
  "OCR_plots",
  paste0("Human_vs_Mouse_OCR_GO_BP_top", top_n, "_combined_clean.pdf")
)

cat("Human input:", human_file, "\n")
cat("Mouse input:", mouse_file, "\n")
cat("Output file:", output_file, "\n")

# ===============================
# 2. Read tables
# ===============================
human_tb <- read.table(
  human_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

mouse_tb <- read.table(
  mouse_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# ===============================
# 3. Process function
# ===============================
process_tb <- function(tb, group_name, top_n,
                       min_gene_set_size = 10,
                       max_gene_set_size = 500,
                       min_observed_gene_hits = 10) {

  tb <- tb[!is.na(tb$description) &
           !is.na(tb$fold_enrichment) &
           !is.na(tb$gene_set_size) &
           !is.na(tb$observed_gene_hits), , drop = FALSE]

  if (nrow(tb) == 0) {
    stop(paste("No valid rows in", group_name))
  }

  # Remove very small or overly broad GO terms
  tb <- tb[tb$gene_set_size >= min_gene_set_size &
           tb$gene_set_size <= max_gene_set_size &
           tb$observed_gene_hits >= min_observed_gene_hits, , drop = FALSE]

  if (nrow(tb) == 0) {
    stop(paste("No rows left after filtering in", group_name))
  }

  # Rank by fold enrichment first, then observed gene hits
  tb <- tb[order(tb$fold_enrichment, tb$observed_gene_hits, decreasing = TRUE), , drop = FALSE]
  tb <- head(tb, top_n)

  # Wrap long labels so they do not crush the panel
  tb$label <- str_wrap(tb$description, width = 35)

  # Create facet-specific y labels
  tb$term_id <- paste(group_name, seq_len(nrow(tb)), sep = "__")
  tb$term_id <- factor(tb$term_id, levels = rev(tb$term_id))

  tb$group <- group_name
  return(tb)
}

human_top <- process_tb(human_tb, "Human OCR", top_n)
mouse_top <- process_tb(mouse_tb, "Mouse OCR", top_n)

plot_tb <- rbind(human_top, mouse_top)

# Create label lookup for y-axis
label_map <- setNames(plot_tb$label, plot_tb$term_id)

# ===============================
# 4. Plot
# ===============================
p <- ggplot(plot_tb, aes(x = fold_enrichment, y = term_id)) +
  geom_point(aes(size = observed_gene_hits, color = fold_enrichment), alpha = 0.9) +
  facet_wrap(~group, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) label_map[x]) +
  scale_size_continuous(range = c(3, 9)) +
  labs(
    title = "Human vs Mouse OCR GO Biological Process Enrichment",
    x = "Fold enrichment",
    y = NULL,
    size = "Observed gene hits",
    color = "Fold enrichment"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(1.5, "lines")
  )

# ===============================
# 5. Save plot
# ===============================
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(
  output_file,
  plot = p,
  width = 16,
  height = 8
)

cat("Plot done.\n")