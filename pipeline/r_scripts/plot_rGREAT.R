# run with:
# Rscript pipeline/r_scripts/plot_rGREAT.R
# Rscript pipeline/r_scripts/plot_rGREAT.R all
# Rscript pipeline/r_scripts/plot_rGREAT.R Human_full
# Rscript pipeline/r_scripts/plot_rGREAT.R Mouse_full
# Rscript pipeline/r_scripts/plot_rGREAT.R HtM_shared
# Rscript pipeline/r_scripts/plot_rGREAT.R HtM_human_specific
# Rscript pipeline/r_scripts/plot_rGREAT.R MtH_shared
# Rscript pipeline/r_scripts/plot_rGREAT.R MtH_mouse_specific

user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

required_packages <- c("ggplot2", "yaml")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Required R package not installed: ", pkg,
        "\nInstall it with: install.packages('", pkg, "', repos='https://cloud.r-project.org')"
      )
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
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

  if (is.null(cfg$rgreat_output_dir)) {
    stop("Missing config key: rgreat_output_dir")
  }

  cfg$rgreat_output_dir <- resolve_path(cfg$project_root, cfg$rgreat_output_dir)

  cfg
}

get_rgreat_tasks <- function() {
  list(
    Human_full = list(
      input_subdir = "Human_full",
      input_label = "human_full",
      plot_title = "Human full OCRs"
    ),
    Mouse_full = list(
      input_subdir = "Mouse_full",
      input_label = "mouse_full",
      plot_title = "Mouse full OCRs"
    ),
    HtM_shared = list(
      input_subdir = "HtM_shared",
      input_label = "htm_shared",
      plot_title = "Human-to-Mouse shared OCRs"
    ),
    HtM_human_specific = list(
      input_subdir = "HtM_human_specific",
      input_label = "htm_human_specific",
      plot_title = "Human-specific OCRs"
    ),
    MtH_shared = list(
      input_subdir = "MtH_shared",
      input_label = "mth_shared",
      plot_title = "Mouse-to-Human shared OCRs"
    ),
    MtH_mouse_specific = list(
      input_subdir = "MtH_mouse_specific",
      input_label = "mth_mouse_specific",
      plot_title = "Mouse-specific OCRs"
    )
  )
}

get_input_file <- function(cfg, task_info) {
  file.path(
    cfg$rgreat_output_dir,
    task_info$input_subdir,
    paste0(task_info$input_label, "_GO_BP.tsv")
  )
}

get_output_file <- function(cfg, task_name, top_n) {
  file.path(
    cfg$rgreat_output_dir,
    "plots",
    paste0(task_name, "_GO_BP_top", top_n, ".pdf")
  )
}

read_enrichment_table <- function(input_file) {
  if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
  }

  read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

find_column <- function(tb, candidates, label) {
  matched <- candidates[candidates %in% colnames(tb)]

  if (length(matched) == 0) {
    stop(
      paste0(
        "Could not find column for ", label, ".\n",
        "Tried: ", paste(candidates, collapse = ", "), "\n",
        "Available columns: ", paste(colnames(tb), collapse = ", ")
      )
    )
  }

  matched[1]
}

replace_zero_pvalues <- function(pvals) {
  pvals <- suppressWarnings(as.numeric(pvals))

  positive_pvals <- pvals[!is.na(pvals) & pvals > 0]

  if (length(positive_pvals) == 0) {
    stop("No positive adjusted p-values found. Cannot compute -log10 adjusted p-values.")
  }

  floor_value <- min(positive_pvals) / 10
  pvals[!is.na(pvals) & pvals == 0] <- floor_value

  pvals
}

clean_enrichment_table <- function(tb) {
  padj_col <- find_column(
    tb,
    candidates = c(
      "GO.Biological.Process.Binom_Adjp_BH",
      "GO.Biological.Process.Hyper_Adjp_BH",
      "p_adjust",
      "p.adjust",
      "padj",
      "adj_p",
      "adjusted_p_value",
      "hyper_p_adjust",
      "binom_p_adjust",
      "hyper_FDR",
      "binom_FDR"
    ),
    label = "adjusted p-value"
  )

  raw_padj <- suppressWarnings(as.numeric(tb[[padj_col]]))

  zero_fraction <- mean(raw_padj == 0, na.rm = TRUE)
  tb$use_fold_plot <- zero_fraction >= 0.5

  tb$p_adjust_plot <- replace_zero_pvalues(raw_padj)

  tb <- tb[!is.na(tb$p_adjust_plot) & tb$p_adjust_plot > 0, , drop = FALSE]

  if (nrow(tb) == 0) {
    stop("No valid rows found after filtering adjusted p-values.")
  }

  tb$neg_log10_padj <- -log10(tb$p_adjust_plot)

  tb
}

select_top_terms <- function(tb, top_n) {
  description_col <- find_column(
    tb,
    candidates = c(
      "GO.Biological.Process.name",
      "description",
      "name",
      "term_name",
      "term",
      "GO_term",
      "go_term",
      "ID"
    ),
    label = "GO term description"
  )

  fold_col <- find_column(
    tb,
    candidates = c(
      "GO.Biological.Process.Binom_Fold_Enrichment",
      "GO.Biological.Process.Hyper_Fold_Enrichment",
      "fold_enrichment",
      "foldEnrichment",
      "fold_enrich",
      "hyper_fold_enrichment",
      "binom_fold_enrichment",
      "observed_fold_enrichment"
    ),
    label = "fold enrichment"
  )

  observed_col <- find_column(
    tb,
    candidates = c(
      "GO.Biological.Process.Binom_Observed_Region_Hits",
      "GO.Biological.Process.Hyper_Observed_Gene_Hits",
      "observed",
      "observed_hits",
      "region_hits",
      "gene_hits"
    ),
    label = "observed hits"
  )

  tb$description_plot <- tb[[description_col]]
  tb$fold_enrichment_plot <- suppressWarnings(as.numeric(tb[[fold_col]]))
  tb$observed_hits_plot <- suppressWarnings(as.numeric(tb[[observed_col]]))

  tb <- tb[
    !is.na(tb$description_plot) &
      !is.na(tb$fold_enrichment_plot) &
      !is.na(tb$observed_hits_plot),
    ,
    drop = FALSE
  ]

  if (nrow(tb) == 0) {
    stop("No valid rows found after filtering plot columns.")
  }

  tb_top <- tb[order(tb$p_adjust_plot, decreasing = FALSE), , drop = FALSE]
  tb_top <- head(tb_top, top_n)

  tb_top$description_plot <- factor(
    tb_top$description_plot,
    levels = rev(tb_top$description_plot)
  )

  tb_top
}

build_go_plot <- function(tb_top, plot_title) {
  use_fold_plot <- unique(tb_top$use_fold_plot)[1]

  if (isTRUE(use_fold_plot)) {
    ggplot(tb_top, aes(x = fold_enrichment_plot, y = description_plot)) +
      geom_point(
        aes(size = observed_hits_plot, color = fold_enrichment_plot),
        alpha = 0.85
      ) +
      labs(
        title = plot_title,
        x = "Fold enrichment",
        y = NULL,
        size = "Observed region hits",
        color = "Fold enrichment"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10)
      )
  } else {
    ggplot(tb_top, aes(x = neg_log10_padj, y = description_plot)) +
      geom_point(
        aes(size = fold_enrichment_plot, color = neg_log10_padj),
        alpha = 0.85
      ) +
      labs(
        title = plot_title,
        x = expression(-log[10]("adjusted p-value")),
        y = NULL,
        size = "Fold enrichment",
        color = expression(-log[10]("adj. p"))
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10)
      )
  }
}

save_go_plot <- function(p, output_file, width = 8, height = 5.5) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  ggsave(
    output_file,
    plot = p,
    width = width,
    height = height,
    units = "in"
  )
}

plot_one_task <- function(cfg, task_name, top_n = 10) {
  tasks <- get_rgreat_tasks()

  if (is.null(tasks[[task_name]])) {
    valid_tasks <- paste(c("all", names(tasks)), collapse = ", ")
    stop(paste("Unknown plot task:", task_name, "\nValid tasks:", valid_tasks))
  }

  task_info <- tasks[[task_name]]
  input_file <- get_input_file(cfg, task_info)
  output_file <- get_output_file(cfg, task_name, top_n)

  cat("====================================\n")
  cat("Task:", task_name, "\n")
  cat("Input file:", input_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("====================================\n")

  tb <- read_enrichment_table(input_file)
  tb <- clean_enrichment_table(tb)
  tb_top <- select_top_terms(tb, top_n)

  p <- build_go_plot(tb_top, task_info$plot_title)
  save_go_plot(p, output_file)

  cat("Plot done.\n")
}

plot_rgreat <- function(task = "all", config_path = "config.yaml", top_n = 10) {
  cfg <- load_config(config_path)
  tasks <- get_rgreat_tasks()

  if (task == "all") {
    for (task_name in names(tasks)) {
      plot_one_task(
        cfg = cfg,
        task_name = task_name,
        top_n = top_n
      )
    }

    cat("All rGREAT plots done.\n")
    return(invisible(NULL))
  }

  plot_one_task(
    cfg = cfg,
    task_name = task,
    top_n = top_n
  )

  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  task <- "all"
} else {
  task <- args[1]
}

plot_rgreat(
  task = task,
  config_path = "config.yaml",
  top_n = 10
)