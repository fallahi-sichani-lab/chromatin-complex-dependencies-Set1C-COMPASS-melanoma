## ============================================================
## density_scatter_utils.R
## 2D density scatter plots with fixed axes, threshold lines, and quadrant fractions
## - Works on *_ln columns already saved in your RDS
## - Reads a thresholds CSV with Marker + thr_ln (or thr_log fallback)
## - Robust to missing treatments and missing thresholds
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(MASS)   # kde2d
})

# ---- safe name normalizer to map thresholds by marker name ----
signal_name_norm <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+", "_")
}

# ---- MATLAB-jet-like palette (for density coloring) ----
jet_colors <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#00007F", "blue", "cyan", "yellow", "red", "#7F0000")
  )(n)
}

# ---- add per-point kernel density estimate (finite-only) ----
add_point_density <- function(df, x_col, y_col, n_grid = 200) {
  xv <- df[[x_col]]
  yv <- df[[y_col]]
  
  dens <- rep(NA_real_, length(xv))
  ok <- is.finite(xv) & is.finite(yv)
  
  if (sum(ok) < 5) {
    df$density <- dens
    return(df)
  }
  
  lims <- c(range(xv[ok]), range(yv[ok]))
  kde  <- MASS::kde2d(xv[ok], yv[ok], n = n_grid, lims = lims)
  
  ix <- findInterval(xv[ok], kde$x, all.inside = TRUE)
  iy <- findInterval(yv[ok], kde$y, all.inside = TRUE)
  dens[ok] <- kde$z[cbind(ix, iy)]
  
  df$density <- dens
  df
}

# ---- sample matched N per Treatment to avoid unequal point counts ----
sample_matched_treat_n <- function(df, cap = 1500L, treat_levels) {
  df <- df %>% dplyr::filter(.data$Treatment %in% treat_levels)
  
  counts <- df %>% dplyr::count(.data$Treatment)
  n_a <- counts$n[counts$Treatment == treat_levels[1]]
  n_b <- counts$n[counts$Treatment == treat_levels[2]]
  
  # If either treatment missing, return empty df (prevents downstream recycling errors)
  if (length(n_a) == 0 || length(n_b) == 0) {
    return(df[0, , drop = FALSE])
  }
  
  n_take <- min(cap, n_a, n_b)
  if (!is.finite(n_take) || n_take < 2) return(df[0, , drop = FALSE])
  
  df %>%
    dplyr::group_by(.data$Treatment) %>%
    dplyr::slice_sample(n = n_take) %>%
    dplyr::ungroup()
}

# ---- formatting helpers ----
fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * x), "%")

pad_limits <- function(lim, pad_frac = 0.04) {
  lim <- as.numeric(lim)
  lim <- lim[is.finite(lim)]
  if (length(lim) != 2) return(c(NA_real_, NA_real_))
  
  rng <- diff(lim)
  if (!is.finite(rng) || rng <= 0) return(lim)
  
  c(lim[1] - pad_frac * rng, lim[2] + pad_frac * rng)
}

# ---- read thresholds CSV into a named vector of thr_ln ----
# expected columns:
#   Marker, thr_ln   (preferred)
# or Marker, thr_log (fallback)
read_thresholds_ln_map <- function(thresholds_csv) {
  thr_df <- readr::read_csv(thresholds_csv, show_col_types = FALSE)
  
  # pick which threshold column exists
  thr_col <- NULL
  if ("thr_ln" %in% names(thr_df)) {
    thr_col <- "thr_ln"
  } else if ("thr_log" %in% names(thr_df)) {
    thr_col <- "thr_log"
  } else {
    stop("Thresholds CSV must have a column named 'thr_ln' (preferred) or 'thr_log'. File: ", thresholds_csv)
  }
  
  thr_df %>%
    dplyr::mutate(Marker = signal_name_norm(.data$Marker)) %>%
    dplyr::transmute(Marker = .data$Marker, thr_ln = as.numeric(.data[[thr_col]])) %>%
    dplyr::distinct(.data$Marker, .data$thr_ln) %>%
    dplyr::filter(is.finite(.data$thr_ln)) %>%
    { stats::setNames(.$thr_ln, .$Marker) }
}

# ---- get threshold on ln scale for a signal name ----
# signal can be "DNA", "DNA_ln", "pRb_807_811", etc.
get_thr_ln <- function(signal, thr_ln_map, overrides_ln = list()) {
  sig2 <- signal %>%
    sub("_ln$", "", .) %>%
    sub("_log$", "", .) %>%
    signal_name_norm()
  
  # override wins
  if (!is.null(overrides_ln[[sig2]]) && is.finite(overrides_ln[[sig2]])) {
    return(as.numeric(overrides_ln[[sig2]]))
  }
  
  v <- thr_ln_map[sig2]
  if (length(v) == 0 || is.na(v)) NA_real_ else as.numeric(v)
}

# ---- core plotter (one cell line, one x/y pair) ----
plot_jet_density_pair <- function(
    df,
    cl_name,
    keep_lines,
    x_col,
    y_col,
    x_lab,
    y_lab,
    thr_ln_map,
    overrides_ln = list(),
    cap = 1500L,
    pad_frac = 0.04,
    treat_levels
) {
  x_sig <- sub("_ln$", "", x_col)
  y_sig <- sub("_ln$", "", y_col)
  
  thr_x <- get_thr_ln(x_sig, thr_ln_map, overrides_ln = overrides_ln)
  thr_y <- get_thr_ln(y_sig, thr_ln_map, overrides_ln = overrides_ln)
  
  # fixed limits across keep_lines + both treatments for this pair
  df_limits <- df %>%
    dplyr::filter(.data$CellLine %in% keep_lines, .data$Treatment %in% treat_levels) %>%
    dplyr::filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  
  if (nrow(df_limits) < 10) {
    return(list(plot = ggplot() + theme_void(), sampled = df[0, , drop = FALSE]))
  }
  
  x_lim <- stats::quantile(df_limits[[x_col]], probs = c(0.005, 0.995), na.rm = TRUE)
  y_lim <- stats::quantile(df_limits[[y_col]], probs = c(0.005, 0.995), na.rm = TRUE)
  
  x_lim <- pad_limits(x_lim, pad_frac = pad_frac)
  y_lim <- pad_limits(y_lim, pad_frac = pad_frac)
  
  df_samp <- df %>%
    dplyr::filter(.data$CellLine == cl_name, .data$Treatment %in% treat_levels) %>%
    dplyr::filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
    sample_matched_treat_n(cap = cap, treat_levels = treat_levels)
  
  # if sampling returns empty (e.g., missing treatment), stop gracefully
  if (nrow(df_samp) < 2) {
    return(list(plot = ggplot() + theme_void(), sampled = df_samp))
  }
  
  df_samp <- df_samp %>%
    dplyr::group_by(.data$Treatment) %>%
    dplyr::group_modify(~ add_point_density(.x, x_col, y_col, n_grid = 200)) %>%
    dplyr::ungroup()
  
  # quadrant annotation only if BOTH thresholds are finite
  ann_df <- NULL
  if (is.finite(thr_x) && is.finite(thr_y)) {
    ann_df <- df_samp %>%
      dplyr::group_by(.data$Treatment) %>%
      dplyr::summarise(
        n  = dplyr::n(),
        LL = mean(.data[[x_col]] <= thr_x & .data[[y_col]] <= thr_y),
        LR = mean(.data[[x_col]] >  thr_x & .data[[y_col]] <= thr_y),
        UL = mean(.data[[x_col]] <= thr_x & .data[[y_col]] >  thr_y),
        UR = mean(.data[[x_col]] >  thr_x & .data[[y_col]] >  thr_y),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        label = paste0(
          "n=", .data$n, "\n",
          "LL ", fmt_pct(.data$LL), " | LR ", fmt_pct(.data$LR), "\n",
          "UL ", fmt_pct(.data$UL), " | UR ", fmt_pct(.data$UR)
        ),
        x = x_lim[1] + 0.02 * diff(x_lim),
        y = y_lim[2] - 0.02 * diff(y_lim)
      )
  }
  
  p <- ggplot2::ggplot(df_samp, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(ggplot2::aes(color = .data$density), size = 0.4, alpha = 0.9) +
    ggplot2::facet_wrap(~Treatment, nrow = 1) +
    ggplot2::scale_color_gradientn(colors = jet_colors(256)) +
    ggplot2::coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    ggplot2::labs(title = NULL, x = x_lab, y = y_lab) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position  = "none",
      plot.margin      = ggplot2::margin(t = 0, r = 2, b = 2, l = 2, unit = "pt"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 0, b = 1)),
      panel.spacing    = grid::unit(0.15, "lines"),
      axis.title       = ggplot2::element_text(size = 10),
      axis.text        = ggplot2::element_text(size = 9)
    )
  
  # add lines only if finite
  if (is.finite(thr_x)) p <- p + ggplot2::geom_vline(xintercept = thr_x, linetype = "dotted", linewidth = 0.4)
  if (is.finite(thr_y)) p <- p + ggplot2::geom_hline(yintercept = thr_y, linetype = "dotted", linewidth = 0.4)
  
  # add annotation only if built
  if (!is.null(ann_df)) {
    p <- p + ggplot2::geom_text(
      data = ann_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3
    )
  }
  
  list(plot = p, sampled = df_samp)
}

# ---- runner/saver for multiple pairs and lines ----
run_save_jet_density_suite <- function(
    df,
    keep_lines,
    plot_lines,
    pairs,
    thresholds_csv,
    out_dir,
    cap_n_scatter = 1500L,
    overrides_ln = list(),
    pad_frac = 0.04,
    treat_levels
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  thr_ln_map <- read_thresholds_ln_map(thresholds_csv)
  
  df_sc <- df %>%
    dplyr::filter(.data$CellLine %in% keep_lines) %>%
    dplyr::mutate(
      CellLine  = as.character(.data$CellLine),
      Treatment = as.character(.data$Treatment)
    )
  
  for (pair in pairs) {
    for (cl in plot_lines) {
      
      out <- plot_jet_density_pair(
        df           = df_sc,
        cl_name      = cl,
        keep_lines   = keep_lines,
        x_col        = pair$x,
        y_col        = pair$y,
        x_lab        = pair$xlab,
        y_lab        = pair$ylab,
        thr_ln_map   = thr_ln_map,
        overrides_ln = overrides_ln,
        cap          = cap_n_scatter,
        pad_frac     = pad_frac,
        treat_levels = treat_levels
      )
      
      pdf_name <- file.path(out_dir, paste0("jetDensityScatter_", pair$stub, "_", cl, ".pdf"))
      
      grDevices::pdf(pdf_name, width = 4.0, height = 2.0, useDingbats = FALSE)
      print(out$plot)
      grDevices::dev.off()
      
      readr::write_csv(
        out$sampled,
        file.path(out_dir, paste0("jetDensityScatter_", pair$stub, "_sampledMatchedN_", cl, ".csv"))
      )
    }
  }
  
  invisible(NULL)
}
