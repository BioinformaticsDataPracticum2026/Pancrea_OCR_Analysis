suppressPackageStartupMessages({
  library(ggplot2)
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
  file.path(
    comp_cfg$base_outdir,
    run_mode,
    paste0(run_mode, "_GO_BP.tsv")
  )
}

get_output_file <- function(comp_cfg, run_mode, top_n) {
  file.path(
    comp_cfg$base_outdir,
    "plots",
    paste0(run_mode, "_GO_BP_top", top_n, ".pdf")
  )
}

read_enrichment_table <- function(input_file) {
  read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
}

clean_enrichment_table <- function(tb) {
  tb <- tb[!is.na(tb$p_adjust) & tb$p_adjust > 0, , drop = FALSE]

  if (nrow(tb) == 0) {
    stop("No valid rows found after filtering p_adjust.")
  }

  tb$neg_log10_padj <- -log10(tb$p_adjust)
  tb
}

select_top_terms <- function(tb, top_n) {
  if (!"description" %in% colnames(tb)) {
    stop("Column 'description' not found in enrichment table.")
  }

  if (!"fold_enrichment" %in% colnames(tb)) {
    stop("Column 'fold_enrichment' not found in enrichment table.")
  }

  tb_top <- tb[order(tb$p_adjust, decreasing = FALSE), , drop = FALSE]
  tb_top <- head(tb_top, top_n)
  tb_top$description <- factor(tb_top$description, levels = rev(tb_top$description))

  tb_top
}

make_plot_title <- function(comparison, run_mode) {
  paste(comparison, gsub("_", "-", run_mode), "OCRs")
}

build_go_plot <- function(tb_top, plot_title) {
  ggplot(tb_top, aes(x = neg_log10_padj, y = description)) +
    geom_point(aes(size = fold_enrichment, color = neg_log10_padj), alpha = 0.85) +
    labs(
      title = plot_title,
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      size = "Fold enrichment",
      color = expression(-log[10]("adj. p"))
    ) +
    theme_bw(base_size = 12)
}

save_go_plot <- function(p, output_file, width = 8, height = 5) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  ggsave(
    output_file,
    plot = p,
    width = width,
    height = height
  )
}

plot_one_comparison_mode <- function(cfg, comparison, run_mode, top_n = 10) {
  validate_run_mode(run_mode)

  comp_cfg <- get_comparison_config(cfg, comparison)
  input_file <- get_input_file(comp_cfg, run_mode)
  output_file <- get_output_file(comp_cfg, run_mode, top_n)
  plot_title <- make_plot_title(comparison, run_mode)

  cat("====================================\n")
  cat("Comparison:", comparison, "\n")
  cat("Run mode:", run_mode, "\n")
  cat("Input file:", input_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("====================================\n")

  tb <- read_enrichment_table(input_file)
  tb <- clean_enrichment_table(tb)
  tb_top <- select_top_terms(tb, top_n)
  p <- build_go_plot(tb_top, plot_title)
  save_go_plot(p, output_file)

  cat("Plot done.\n")
}

plot_all_comparison_modes <- function(config_path = "config.yaml",
                                      comparisons = NULL,
                                      run_modes = c("shared", "not_open"),
                                      top_n = 10) {
  cfg <- load_config(config_path)

  if (is.null(comparisons)) {
    comparisons <- names(cfg$rgreat$comparisons)
  }

  for (comparison in comparisons) {
    for (run_mode in run_modes) {
      plot_one_comparison_mode(
        cfg = cfg,
        comparison = comparison,
        run_mode = run_mode,
        top_n = top_n
      )
    }
  }

  cat("All plots done.\n")
}

# ===============================
# Call at the end
# ===============================
plot_all_comparison_modes(
  config_path = "config.yaml",
  comparisons = c("HtM", "MtH"),
  run_modes = c("shared", "not_open"),
  top_n = 10
)