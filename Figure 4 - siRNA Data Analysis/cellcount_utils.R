## ============================================================
## cellcount_utils.R
## Single-cell IF (Experiment 1): per-well cell count utilities
## ============================================================
##
## Purpose
## - Standardize per-well cell count tables (labels, factor order)
## - Normalize to si-NTC baseline per cell line
## - Compute summary statistics and p-values (replicate unit = Well)
## - Generate a consistent bar plot (mean ± SE) for normalized cell counts
##
## Expected input format (per-well counts)
## - Required columns:
##   CellLine   : character/factor
##   Treatment  : character/factor (expected: "si-NTC" and "si-CXXC1")
##   Well       : character (replicate unit)
##   n_cells    : numeric (cell count per well)
##
## Normalized output format
## - Adds:
##   baseline_ntc : numeric (mean n_cells in si-NTC wells for that cell line)
##   norm_pct     : numeric (100 * n_cells / baseline_ntc)
##
## Notes
## - Utilities only: no file I/O here
## - Replicate unit is Well throughout (consistent with signal_utils)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

## ------------------------------------------------------------
## Standardize treatment labels for cell counts
## Keeps the "cell-count convention": si-NTC and si-CXXC1
## ------------------------------------------------------------
normalize_treatment_cellcount <- function(x) {
  x <- as.character(x)
  
  out <- dplyr::case_when(
    x %in% c("si-NTC", "si_NTC", "NTC", "ntc") ~ "si-NTC",
    x %in% c("si-CXXC1", "si_CXXC1", "siCXXC1", "si-cxxc1", "siCxxc1") ~ "si-CXXC1",
    TRUE ~ x
  )
  
  factor(out, levels = c("si-NTC", "si-CXXC1"))
}

## ------------------------------------------------------------
## Enforce required columns, clean types, apply cell-line ordering
## ------------------------------------------------------------
standardize_cellcount_df <- function(df, cellline_order = NULL) {
  req <- c("CellLine", "Treatment", "Well", "n_cells")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  df <- df %>%
    mutate(
      CellLine  = as.character(CellLine),
      Treatment = normalize_treatment_cellcount(Treatment),
      Well      = as.character(Well),
      n_cells   = as.numeric(n_cells)
    ) %>%
    filter(!is.na(CellLine), !is.na(Treatment), !is.na(Well), is.finite(n_cells))
  
  if (!is.null(cellline_order)) {
    df <- df %>%
      filter(CellLine %in% cellline_order) %>%
      mutate(CellLine = factor(CellLine, levels = cellline_order))
  } else {
    df <- df %>%
      mutate(CellLine = factor(CellLine, levels = unique(CellLine)))
  }
  
  df
}

## ------------------------------------------------------------
## Normalize per-well counts to si-NTC baseline within each cell line
## - baseline_ntc is computed as mean(n_cells) over si-NTC wells
## - norm_pct = 100 * n_cells / baseline_ntc
## ------------------------------------------------------------
normalize_cell_counts <- function(df, baseline_treatment = "si-NTC") {
  # basic check: n_cells must exist
  if (!("n_cells" %in% names(df))) stop("normalize_cell_counts: df must contain 'n_cells'")
  
  df %>%
    group_by(CellLine) %>%
    mutate(
      baseline_ntc = mean(n_cells[Treatment == baseline_treatment], na.rm = TRUE),
      norm_pct     = 100 * n_cells / baseline_ntc
    ) %>%
    ungroup()
}

## ------------------------------------------------------------
## Summarize mean ± SE of normalized counts per CellLine × Treatment
## Also compute Welch t-test per cell line on norm_pct (replicate unit = Well)
## ------------------------------------------------------------
summarise_norm_cell_counts <- function(df_norm) {
  if (!("norm_pct" %in% names(df_norm))) {
    stop("summarise_norm_cell_counts: expected column 'norm_pct'. Did you run normalize_cell_counts()?")
  }
  
  means_df <- df_norm %>%
    group_by(CellLine, Treatment) %>%
    summarise(
      mean_pct = mean(norm_pct, na.rm = TRUE),
      se_pct   = sd(norm_pct, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups  = "drop"
    )
  
  pvals_df <- df_norm %>%
    group_by(CellLine) %>%
    summarise(
      p_value = {
        dat <- dplyr::pick(norm_pct, Treatment)
        if (n_distinct(dat$Treatment) < 2 || nrow(dat) < 3) NA_real_
        else tryCatch(t.test(norm_pct ~ Treatment, data = dat)$p.value,
                      error = function(e) NA_real_)
      },
      max_pct = max(norm_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(
        is.na(p_value),
        "p = NA",
        paste0("p = ", format.pval(p_value, digits = 2, eps = 1e-300))
      ),
      y_pos = pmin(max_pct + 5, 0.95 * 120)  # safe default; overridden by y_max in plot if needed
    )
  
  list(means_df = means_df, pvals_df = pvals_df)
}

## ------------------------------------------------------------
## Bar plot: normalized cell count (mean ± SE) with per-cell-line p-values
## Styling matches signal_utils barplots
## ------------------------------------------------------------
plot_norm_cellcount_bar <- function(means_df,
                                    pvals_df,
                                    pal_treat,
                                    ylab = "% normalized cell count (vs si-NTC)",
                                    y_max = 120) {
  # keep CellLine ordering consistent with means_df factor levels
  if (!is.factor(means_df$CellLine)) {
    means_df <- means_df %>% mutate(CellLine = factor(CellLine, levels = unique(CellLine)))
  }
  pvals_df <- pvals_df %>%
    mutate(
      CellLine = factor(CellLine, levels = levels(means_df$CellLine)),
      y_pos    = pmin(y_pos, 0.92 * y_max)
    )
  
  ggplot(means_df, aes(x = CellLine, y = mean_pct, fill = Treatment)) +
    geom_col(
      position  = position_dodge(width = 0.6),
      width     = 0.5,
      colour    = "black",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
      position  = position_dodge(width = 0.6),
      width     = 0.2,
      linewidth = 0.3
    ) +
    geom_text(
      data = pvals_df,
      aes(x = CellLine, y = y_pos, label = p_label),
      inherit.aes = FALSE,
      vjust = 0,
      size  = 3
    ) +
    labs(x = NULL, y = ylab) +
    scale_fill_manual(values = pal_treat) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

## ------------------------------------------------------------
## Convenience: save a ggplot to PDF with consistent defaults
## ------------------------------------------------------------
save_pdf <- function(plot, out_dir, filename, width = 3, height = 2.5) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    filename    = file.path(out_dir, filename),
    plot        = plot,
    width       = width,
    height      = height,
    units       = "in",
    useDingbats = FALSE
  )
}