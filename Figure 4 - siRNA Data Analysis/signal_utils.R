## ============================================================
## signal_utils.R
## Single-cell IF utilities: standardization, summary statistics, plotting
## ============================================================
##
## Expected input format (both experiments)
## - Required columns:
##   CellLine   : character/factor
##   Treatment  : character/factor (normalized to "NTC" and "si-CXXC1")
##   Well       : character (replicate unit for statistics)
## - Optional columns (left untouched):
##   Replicate, Site, Cell, etc.
## - Signal columns expected by downstream analysis:
##   <marker>_scaled : numeric (scaled intensity; e.g., x*(2^16-1))
##   <marker>_ln     : numeric (natural log of scaled intensity)
##
## Notes
## - Utilities only (no file I/O)
## - Designed to be reused across experiments

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

## ------------------------------------------------------------
## Treatment label normalization
## ------------------------------------------------------------
normalize_treatment <- function(x) {
  x <- as.character(x)
  
  out <- dplyr::case_when(
    x %in% c("si-NTC", "si_NTC", "NTC", "ntc") ~ "NTC",
    x %in% c("si-CXXC1", "si_CXXC1", "siCXXC1", "si-cxxc1", "siCxxc1") ~ "si-CXXC1",
    TRUE ~ x
  )
  
  factor(out, levels = c("NTC", "si-CXXC1"))
}

## ------------------------------------------------------------
## Standardize dataset columns and apply requested cell-line order
## ------------------------------------------------------------
standardize_single_cell_df <- function(df, cellline_order = NULL) {
  req <- c("CellLine", "Treatment", "Well")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  df <- df %>%
    mutate(
      CellLine  = as.character(CellLine),
      Treatment = normalize_treatment(Treatment),
      Well      = as.character(Well)
    )
  
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
## Validate required marker columns: <marker>_ln and <marker>_scaled
## ------------------------------------------------------------
validate_marker_columns <- function(df, markers) {
  need_ln     <- paste0(markers, "_ln")
  need_scaled <- paste0(markers, "_scaled")
  
  miss_ln     <- setdiff(need_ln, names(df))
  miss_scaled <- setdiff(need_scaled, names(df))
  
  if (length(miss_ln) > 0) stop("Missing *_ln columns: ", paste(miss_ln, collapse = ", "))
  if (length(miss_scaled) > 0) stop("Missing *_scaled columns: ", paste(miss_scaled, collapse = ", "))
  
  invisible(TRUE)
}

## ------------------------------------------------------------
## Per-well Welch t-test p-values for violin plots (continuous ln signal)
## Replicate unit = Well; test per CellLine: NTC vs si-CXXC1 on per-well mean
## ------------------------------------------------------------
compute_violin_pvals <- function(df, signal_ln_col) {
  per_well <- df %>%
    group_by(CellLine, Treatment, Well) %>%
    summarise(mean_ln = mean(.data[[signal_ln_col]], na.rm = TRUE), .groups = "drop")
  
  pvals <- per_well %>%
    group_by(CellLine) %>%
    summarise(
      p_value = {
        if (n_distinct(Treatment) < 2 || n() < 3) NA_real_
        else tryCatch(t.test(mean_ln ~ Treatment)$p.value, error = function(e) NA_real_)
      },
      .groups = "drop"
    )
  
  ymax <- df %>%
    group_by(CellLine) %>%
    summarise(max_y = max(.data[[signal_ln_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(y_pos = max_y + 0.25)
  
  pvals %>%
    left_join(ymax, by = "CellLine") %>%
    mutate(
      p_label = ifelse(
        is.na(p_value),
        "p = NA",
        paste0("p = ", format.pval(p_value, digits = 3, eps = 1e-300))
      )
    )
}

## ------------------------------------------------------------
## Violin + inner boxplot + p-value labels
## Optional: overlay global ln-threshold line (thr_ln)
## ------------------------------------------------------------
plot_violin_signal <- function(df,
                               signal_ln_col,
                               ylab,
                               pvals_df,
                               pal_treat,
                               thr_ln = NULL) {
  pos_dodge <- position_dodge(width = 0.75)
  
  p <- ggplot(df, aes(x = CellLine, y = .data[[signal_ln_col]], fill = Treatment)) +
    geom_violin(
      position  = pos_dodge,
      trim      = TRUE,
      alpha     = 0.8,
      linewidth = 0.3,
      colour    = "black"
    ) +
    geom_boxplot(
      position      = pos_dodge,
      width         = 0.18,
      alpha         = 1,
      outlier.shape = NA,
      linewidth     = 0.3,
      colour        = "black",
      coef          = 0
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
    theme_classic() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  if (!is.null(thr_ln) && is.finite(thr_ln)) {
    p <- p + geom_hline(yintercept = thr_ln, linetype = "dashed", linewidth = 0.3)
  }
  
  p
}

## ------------------------------------------------------------
## Summarize % high across wells + Welch t-test on per-well % high
## high_col is a logical column (TRUE/FALSE per cell)
## ------------------------------------------------------------
summarise_high_replicates <- function(df, high_col) {
  per_well <- df %>%
    group_by(CellLine, Treatment, Well) %>%
    summarise(
      n      = sum(!is.na(.data[[high_col]])),
      n_high = sum(.data[[high_col]], na.rm = TRUE),
      frac   = n_high / n,
      pct    = 100 * frac,
      .groups = "drop"
    )
  
  means_df <- per_well %>%
    group_by(CellLine, Treatment) %>%
    summarise(
      mean_pct = mean(pct, na.rm = TRUE),
      se_pct   = sd(pct,  na.rm = TRUE) / sqrt(dplyr::n()),
      .groups  = "drop"
    )
  
  pvals_df <- per_well %>%
    group_by(CellLine) %>%
    summarise(
      p_value = {
        if (n_distinct(Treatment) < 2 || n() < 3) NA_real_
        else tryCatch(t.test(pct ~ Treatment)$p.value, error = function(e) NA_real_)
      },
      max_pct = max(pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(
        is.na(p_value),
        "p = NA",
        paste0("p = ", format.pval(p_value, digits = 3, eps = 1e-300))
      ),
      y_pos = max_pct + 5
    )
  
  list(per_well = per_well, means_df = means_df, pvals_df = pvals_df)
}

## ------------------------------------------------------------
## Barplot: mean ± SE across wells with p-value label per cell line
## ------------------------------------------------------------
plot_bar_pct <- function(means_df, pvals_df, pal_treat, ylab = "% High cells") {
  lvls <- if (is.factor(means_df$CellLine)) levels(means_df$CellLine) else unique(means_df$CellLine)
  
  means_df <- means_df %>%
    mutate(CellLine = factor(as.character(CellLine), levels = lvls))
  
  pvals_df <- pvals_df %>%
    mutate(
      CellLine = factor(as.character(CellLine), levels = lvls),
      y_pos    = pmin(y_pos, 98)
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
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

## ------------------------------------------------------------
## Save a ggplot to PDF with consistent defaults
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