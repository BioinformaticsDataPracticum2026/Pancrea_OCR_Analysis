user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

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
      plot_title = "HtM shared OCRs"
    ),
    HtM_human_specific = list(
      input_subdir = "HtM_human_specific",
      input_label = "htm_human_specific",
      plot_title = "HtM human-specific OCRs"
    ),
    MtH_shared = list(
      input_subdir = "MtH_shared",
      input_label = "mth_shared",
      plot_title = "MtH shared OCRs"
    ),
    MtH_mouse_specific = list(
      input_subdir = "MtH_mouse_specific",
      input_label = "mth_mouse_specific",
      plot_title = "MtH mouse-specific OCRs"
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
    stringsAsFactors = FALSE
  )
}

clean_enrichment_table <- function(tb) {
  if (!"p_adjust" %in% colnames(tb)) {
    stop("Column 'p_adjust' not found in enrichment table.")
  }

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