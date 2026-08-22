make_interaction_forest <- function(results, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to make publication figures.", call. = FALSE)
  }
  plot_data <- results[results$model_status == "ok", , drop = FALSE]
  if (nrow(plot_data) == 0) {
    stop("No successful model rows are available for plotting.", call. = FALSE)
  }
  plot_data$label <- paste(plot_data$cohort, plot_data$outcome_label, plot_data$model_id, sep = " | ")
  plot_data$label <- factor(plot_data$label, levels = rev(unique(plot_data$label)))

  figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = interaction_beta, y = label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(x = "Standardized pTau217 x PVeD interaction coefficient", y = NULL) +
    ggplot2::theme_classic(base_size = 8)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output_path, figure, width = 7, height = max(3, 0.22 * nrow(plot_data) + 1), units = "in", dpi = 600)
  invisible(output_path)
}
