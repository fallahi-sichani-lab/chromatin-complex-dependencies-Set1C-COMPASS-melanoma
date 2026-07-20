## ============================================================
## regression_utils.R
## Utilities for Δ%-based regression analysis (Experiment 1)
## ============================================================
##
## Purpose
## - Standardize treatment labels across pipelines (IF %High vs cell counts)
## - Compute per-cell-line Δ% = (si-CXXC1 - si-NTC) from summary tables
## - Generate regression plots with consistent styling and PDF export
##
## Expected treatment conventions
## - Inputs may contain: "NTC", "si-NTC", "si_NTC", "ntc"  -> standardized to "si-NTC"
## - Inputs may contain: "si-CXXC1", "si_CXXC1"           -> standardized to "si-CXXC1"
##
## Required packages (loaded by caller)
## - dplyr, tidyr, ggplot2, broom
## ============================================================

#' Standardize treatment labels to {si-NTC, si-CXXC1}
#'
#' @param x Character vector of treatment labels.
#' @return Character vector with standardized labels.
standardize_treatment_labels <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("si-NTC", "si_NTC", "NTC", "ntc") ~ "si-NTC",
    x %in% c("si-CXXC1", "si_CXXC1")          ~ "si-CXXC1",
    TRUE ~ x
  )
}

#' Compute Δ% (si-CXXC1 - si-NTC) from a means table
#'
#' @param means_df Data frame with columns: CellLine, Treatment, and value_col.
#' @param value_col Column name containing the numeric value to delta (default "mean_pct").
#' @return Tibble with columns: CellLine, delta.
compute_delta_from_means <- function(means_df, value_col = "mean_pct") {
  req <- c("CellLine", "Treatment", value_col)
  miss <- setdiff(req, names(means_df))
  if (length(miss) > 0) {
    stop("means_df missing required columns: ", paste(miss, collapse = ", "))
  }
  
  df2 <- means_df %>%
    dplyr::mutate(
      Treatment = standardize_treatment_labels(.data$Treatment)
    )
  
  have <- sort(unique(df2$Treatment))
  need <- c("si-NTC", "si-CXXC1")
  if (!all(need %in% have)) {
    stop(
      "Treatment labels after standardization do not include both si-NTC and si-CXXC1.\n",
      "Found: ", paste(have, collapse = ", ")
    )
  }
  
  df2 %>%
    dplyr::select(.data$CellLine, .data$Treatment, !!rlang::sym(value_col)) %>%
    tidyr::pivot_wider(names_from = .data$Treatment, values_from = !!rlang::sym(value_col)) %>%
    dplyr::mutate(delta = .data[["si-CXXC1"]] - .data[["si-NTC"]]) %>%
    dplyr::select(.data$CellLine, .data$delta)
}

#' Build regression input table for Exp1 from marker %High means and cell-count norm_pct
#'
#' @param cxxc1_means Data frame: CellLine, Treatment, mean_pct for CXXC1 %High
#' @param ki67_means  Data frame: CellLine, Treatment, mean_pct for Ki67 %High
#' @param h3k4_means  Data frame: CellLine, Treatment, mean_pct for H3K4me3 %High
#' @param cell_df     Data frame containing per-well normalized cell count with norm_pct
#' @return Tibble with per-cell-line deltas for all metrics
build_exp1_delta_table <- function(cxxc1_means, ki67_means, h3k4_means, cell_df) {
  # Validate marker means tables
  req_means_cols <- c("CellLine", "Treatment", "mean_pct")
  for (nm in c("cxxc1_means", "ki67_means", "h3k4_means")) {
    obj <- get(nm)
    miss <- setdiff(req_means_cols, names(obj))
    if (length(miss) > 0) stop(nm, " is missing columns: ", paste(miss, collapse = ", "))
  }
  
  # Validate cell-count table
  req_cell_cols <- c("CellLine", "Treatment", "Well", "norm_pct")
  miss_cell <- setdiff(req_cell_cols, names(cell_df))
  if (length(miss_cell) > 0) {
    stop("cell-count table is missing columns: ", paste(miss_cell, collapse = ", "))
  }
  
  # Standardize treatment labels in cell-count table
  cell_df2 <- cell_df %>%
    dplyr::mutate(Treatment = standardize_treatment_labels(.data$Treatment))
  
  # Cell count: mean(norm_pct) per CellLine × Treatment, then Δ%
  cell_means <- cell_df2 %>%
    dplyr::group_by(.data$CellLine, .data$Treatment) %>%
    dplyr::summarise(mean_norm_pct = mean(.data$norm_pct, na.rm = TRUE), .groups = "drop")
  
  delta_norm <- compute_delta_from_means(cell_means, value_col = "mean_norm_pct") %>%
    dplyr::rename(delta_norm_cellcount = .data$delta)
  
  # Marker deltas from saved means tables
  delta_cxxc1 <- compute_delta_from_means(cxxc1_means, value_col = "mean_pct") %>%
    dplyr::rename(delta_CXXC1_high = .data$delta)
  
  delta_h3k4  <- compute_delta_from_means(h3k4_means, value_col = "mean_pct") %>%
    dplyr::rename(delta_H3K4me3_high = .data$delta)
  
  delta_ki67  <- compute_delta_from_means(ki67_means, value_col = "mean_pct") %>%
    dplyr::rename(delta_Ki67_high = .data$delta)
  
  delta_cxxc1 %>%
    dplyr::inner_join(delta_h3k4, by = "CellLine") %>%
    dplyr::inner_join(delta_ki67, by = "CellLine") %>%
    dplyr::inner_join(delta_norm, by = "CellLine")
}

#' Save a ggplot to PDF using Cairo (consistent device across systems)
#'
#' @param plot ggplot object
#' @param out_dir Output directory
#' @param filename Output filename
#' @param width,height Dimensions in inches
save_pdf_regression <- function(plot, out_dir, filename, width = 1.5, height = 1.5) {
  ggplot2::ggsave(
    filename = file.path(out_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    units    = "in",
    device   = grDevices::cairo_pdf
  )
}

#' Fit linear model and plot regression with labeled points
#'
#' @param df Data frame containing xvar, yvar, CellLine, Shape
#' @param xvar,yvar Column names (strings)
#' @param xlab,ylab Axis labels (can be expression())
#' @param filename PDF filename
#' @param out_dir Output directory
#' @param ylim Optional numeric vector of length 2 for matched limits
#' @return ggplot object
plot_regression_lm <- function(df, xvar, yvar, xlab, ylab, filename, out_dir, ylim = NULL) {
  stopifnot(xvar %in% names(df), yvar %in% names(df))
  if (!("Shape" %in% names(df))) stop("df must contain a 'Shape' column (numeric shape codes).")
  if (!("CellLine" %in% names(df))) stop("df must contain a 'CellLine' column for labels.")
  
  model <- stats::lm(df[[yvar]] ~ df[[xvar]])
  stats_gl <- broom::glance(model)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]])) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$Shape), size = 2, color = "black", show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = .data$CellLine), vjust = -1, size = 2) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      subtitle = paste0(
        "R² = ", formatC(stats_gl$r.squared, format = "f", digits = 3),
        "   |   p = ", signif(stats_gl$p.value, 2)
      )
    ) +
    ggplot2::scale_shape_identity() +
    ggplot2::theme_classic(base_size = 8)
  
  if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
  
  save_pdf_regression(p, out_dir, filename, width = 1.5, height = 1.5)
  p
}

#' Run the Experiment 1 regression suite
#'
#' @param rds_cxxc1_means RDS path for CXXC1 percentHigh means table
#' @param rds_ki67_means  RDS path for Ki67 percentHigh means table
#' @param rds_h3k4_means  RDS path for H3K4me3 percentHigh means table
#' @param rds_cellcount   RDS path for per-well cell-count table containing norm_pct
#' @param out_dir Output directory for regression PDFs and CSV
#' @param sensitive_lines Character vector of sensitive cell lines (triangles)
#' @return List containing delta_df and plots
run_exp1_regression_suite <- function(
    rds_cxxc1_means,
    rds_ki67_means,
    rds_h3k4_means,
    rds_cellcount,
    out_dir,
    sensitive_lines = c("UACC62", "MALME3M", "LOXIMVI")
) {
  for (p in c(rds_cxxc1_means, rds_ki67_means, rds_h3k4_means, rds_cellcount)) {
    if (!file.exists(p)) stop("Missing: ", p)
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  cxxc1_means <- readRDS(rds_cxxc1_means)
  ki67_means  <- readRDS(rds_ki67_means)
  h3k4_means  <- readRDS(rds_h3k4_means)
  cell_df     <- readRDS(rds_cellcount)
  
  delta_df <- build_exp1_delta_table(
    cxxc1_means = cxxc1_means,
    ki67_means  = ki67_means,
    h3k4_means  = h3k4_means,
    cell_df     = cell_df
  )
  
  # Save regression input table
  utils::write.csv(
    delta_df,
    file.path(out_dir, "delta_summary_regression_input.csv"),
    row.names = FALSE
  )
  
  # Add shapes
  delta_df <- delta_df %>%
    dplyr::mutate(Shape = ifelse(.data$CellLine %in% sensitive_lines, 17, 16))
  
  # Global limits (pairwise consistency)
  ylim_ki67 <- range(delta_df$delta_Ki67_high, na.rm = TRUE)
  ylim_norm <- range(delta_df$delta_norm_cellcount, na.rm = TRUE)
  
  # Run regressions
  p1 <- plot_regression_lm(
    df       = delta_df,
    xvar     = "delta_CXXC1_high",
    yvar     = "delta_norm_cellcount",
    xlab     = expression(Delta * "% CXXC1"["High"]),
    ylab     = expression(Delta * "% normalized cell count"),
    filename = "reg_delta_CXXC1high_vs_delta_normCellCount.pdf",
    out_dir  = out_dir,
    ylim     = ylim_norm
  )
  
  p2 <- plot_regression_lm(
    df       = delta_df,
    xvar     = "delta_H3K4me3_high",
    yvar     = "delta_norm_cellcount",
    xlab     = expression(Delta * "% H3K4me3"["High"]),
    ylab     = expression(Delta * "% normalized cell count"),
    filename = "reg_delta_H3K4me3high_vs_delta_normCellCount.pdf",
    out_dir  = out_dir,
    ylim     = ylim_norm
  )
  
  p3 <- plot_regression_lm(
    df       = delta_df,
    xvar     = "delta_CXXC1_high",
    yvar     = "delta_Ki67_high",
    xlab     = expression(Delta * "% CXXC1"["High"]),
    ylab     = expression(Delta * "% Ki67"["High"]),
    filename = "reg_delta_CXXC1high_vs_delta_Ki67high.pdf",
    out_dir  = out_dir,
    ylim     = ylim_ki67
  )
  
  p4 <- plot_regression_lm(
    df       = delta_df,
    xvar     = "delta_H3K4me3_high",
    yvar     = "delta_Ki67_high",
    xlab     = expression(Delta * "% H3K4me3"["High"]),
    ylab     = expression(Delta * "% Ki67"["High"]),
    filename = "reg_delta_H3K4me3high_vs_delta_Ki67high.pdf",
    out_dir  = out_dir,
    ylim     = ylim_ki67
  )
  
  list(
    delta_df = delta_df,
    plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4),
    ylim = list(ylim_ki67 = ylim_ki67, ylim_norm = ylim_norm)
  )
}