
## ============================================================
## gmm_utils.R
## Gaussian mixture modeling utilities for single-cell IF
## ============================================================
##
## Expected input format
## - Signal columns should already exist as:
##   <marker>_scaled : numeric (e.g., x*(2^16-1))
##   <marker>_ln     : numeric (natural log of scaled intensity)
## - Metadata columns used for per-cell-line thresholds:
##   CellLine : character/factor
##
## Notes
## - Utilities only: no file I/O here (except optional ggsave in runner)
## - Uses mclust with exactly G = 2 components
## - Threshold definition: posterior crossing (more robust than midpoint)
## - Includes robust column validation to prevent length-0 assignment errors
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(mclust)
  library(ggplot2)
})

## ------------------------------------------------------------
## Internal helpers
## ------------------------------------------------------------

# Validate that required marker columns exist
validate_marker_columns <- function(df, markers, require = c("ln"), verbose = TRUE) {
  require <- match.arg(require, choices = c("ln", "scaled", "either"), several.ok = FALSE)
  
  need_ln     <- paste0(markers, "_ln")
  need_scaled <- paste0(markers, "_scaled")
  
  has_ln     <- need_ln %in% names(df)
  has_scaled <- need_scaled %in% names(df)
  
  if (verbose) {
    msg <- tibble(
      Marker = markers,
      has_ln = has_ln,
      has_scaled = has_scaled
    )
    message("[GMM] Marker column check:")
    print(msg)
  }
  
  if (require == "ln") {
    miss <- markers[!has_ln]
    if (length(miss) > 0) {
      stop("Missing required *_ln columns for markers: ", paste(miss, collapse = ", "))
    }
  } else if (require == "scaled") {
    miss <- markers[!has_scaled]
    if (length(miss) > 0) {
      stop("Missing required *_scaled columns for markers: ", paste(miss, collapse = ", "))
    }
  } else { # either
    miss <- markers[!(has_ln | has_scaled)]
    if (length(miss) > 0) {
      stop("Missing both *_ln and *_scaled columns for markers: ", paste(miss, collapse = ", "))
    }
  }
  
  invisible(TRUE)
}

# Get ln-scale vector for a marker (prefer *_ln, else log(*_scaled))
get_marker_ln <- function(df, marker) {
  col_ln     <- paste0(marker, "_ln")
  col_scaled <- paste0(marker, "_scaled")
  
  if (col_ln %in% names(df)) {
    return(df[[col_ln]])
  }
  if (col_scaled %in% names(df)) {
    x <- df[[col_scaled]]
    return(log(x))
  }
  
  stop("Marker '", marker, "' missing both columns: ", col_ln, " and ", col_scaled)
}

# Get scaled-scale vector for a marker (prefer *_scaled, else exp(*_ln))
get_marker_scaled <- function(df, marker) {
  col_ln     <- paste0(marker, "_ln")
  col_scaled <- paste0(marker, "_scaled")
  
  if (col_scaled %in% names(df)) {
    return(df[[col_scaled]])
  }
  if (col_ln %in% names(df)) {
    x <- df[[col_ln]]
    return(exp(x))
  }
  
  stop("Marker '", marker, "' missing both columns: ", col_ln, " and ", col_scaled)
}

## ------------------------------------------------------------
## Core: Fit 2-component GMM + posterior-crossing threshold
## ------------------------------------------------------------

fit_gmm_threshold <- function(x_ln, marker_name = "") {
  x_ln <- x_ln[is.finite(x_ln)]
  if (length(x_ln) < 20) stop("Not enough finite values for GMM: ", marker_name)
  
  m <- Mclust(x_ln, G = 2, verbose = FALSE)
  
  mu      <- as.numeric(m$parameters$mean)
  sigmasq <- m$parameters$variance$sigmasq
  if (is.matrix(sigmasq)) sigmasq <- as.numeric(sigmasq)
  sig  <- sqrt(as.numeric(sigmasq))
  prop <- as.numeric(m$parameters$pro)
  
  if (length(mu) != 2L) stop("mclust did not return 2 components for: ", marker_name)
  
  ord  <- order(mu) # low -> high mean
  grid <- seq(min(x_ln), max(x_ln), length.out = 1000)
  
  dens_low  <- prop[ord[1]] * dnorm(grid, mean = mu[ord[1]], sd = sig[ord[1]])
  dens_high <- prop[ord[2]] * dnorm(grid, mean = mu[ord[2]], sd = sig[ord[2]])
  
  tot <- dens_low + dens_high
  post_low  <- dens_low  / tot
  post_high <- dens_high / tot
  
  diff_post <- post_low - post_high
  idx <- which(diff_post[-1] * diff_post[-length(diff_post)] < 0)
  
  thr_ln <- if (length(idx) > 0) {
    i  <- idx[1]
    x1 <- grid[i]; x2 <- grid[i + 1]
    y1 <- diff_post[i]; y2 <- diff_post[i + 1]
    x1 - y1 * (x2 - x1) / (y2 - y1)
  } else {
    grid[which.min(abs(diff_post))]
  }
  
  list(
    model      = m,
    grid       = grid,
    ord        = ord,
    mu         = mu,
    sig        = sig,
    prop       = prop,
    thr_ln     = thr_ln,
    thr_scaled = exp(thr_ln),
    dens1      = dens_low,
    dens2      = dens_high
  )
}

## ------------------------------------------------------------
## Plot a 1D mixture on ln-scale
## ------------------------------------------------------------

plot_gmm_mixture_1d <- function(x_ln, fit, xlab, title = "") {
  x_ln <- x_ln[is.finite(x_ln)]
  
  dens_df <- tibble(
    x       = c(fit$grid, fit$grid),
    density = c(fit$dens1, fit$dens2),
    comp    = factor(rep(c("Low", "High"), each = length(fit$grid)),
                     levels = c("Low", "High"))
  )
  
  ggplot() +
    geom_histogram(
      data   = tibble(x = x_ln),
      aes(x = x, y = after_stat(density)),
      bins   = 60,
      fill   = "grey85",
      colour = "grey85",
      alpha  = 0.6
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = density, colour = comp),
      linewidth = 1
    ) +
    scale_colour_manual(values = c("Low" = "#8770E0", "High" = "#85D484")) +
    geom_vline(
      xintercept = fit$thr_ln,
      linetype   = "dashed",
      linewidth  = 1,
      colour     = "red"
    ) +
    labs(x = xlab, y = "Density", title = title) +
    theme_classic() +
    theme(legend.position = "none")
}

## ------------------------------------------------------------
## Threshold tables
## ------------------------------------------------------------

# Global (pooled) thresholds per marker
compute_global_thresholds <- function(df, markers, verbose = TRUE) {
  validate_marker_columns(df, markers, require = "either", verbose = verbose)
  
  out <- vector("list", length(markers))
  names(out) <- markers
  
  for (mk in markers) {
    x_ln <- get_marker_ln(df, mk)
    fit  <- fit_gmm_threshold(x_ln, marker_name = mk)
    
    out[[mk]] <- tibble(
      Marker      = mk,
      thr_ln      = fit$thr_ln,
      thr_scaled  = fit$thr_scaled,
      mu_low      = sort(fit$mu)[1],
      mu_high     = sort(fit$mu)[2]
    )
  }
  
  bind_rows(out)
}

# Per-cell-line thresholds per marker (fits on ln-scale per cell line)
compute_cellline_thresholds <- function(df, markers, cellline_col = "CellLine", min_n = 20, verbose = TRUE) {
  if (!cellline_col %in% names(df)) stop("Missing column: ", cellline_col)
  validate_marker_columns(df, markers, require = "either", verbose = verbose)
  
  need_ln <- paste0(markers, "_ln")
  need_scaled <- paste0(markers, "_scaled")
  
  # Keep whichever exists; pivot longer across existing columns
  keep_cols <- c()
  for (mk in markers) {
    if (paste0(mk, "_ln") %in% names(df)) {
      keep_cols <- c(keep_cols, paste0(mk, "_ln"))
    } else if (paste0(mk, "_scaled") %in% names(df)) {
      keep_cols <- c(keep_cols, paste0(mk, "_scaled"))
    }
  }
  
  df %>%
    mutate(.cellline = as.character(.data[[cellline_col]])) %>%
    select(.cellline, all_of(keep_cols)) %>%
    pivot_longer(cols = - .cellline, names_to = "ColName", values_to = "x") %>%
    mutate(
      Marker = sub("_(ln|scaled)$", "", ColName),
      Scale  = sub("^.*_(ln|scaled)$", "\\1", ColName),
      x_ln   = dplyr::if_else(Scale == "ln", as.numeric(x), log(as.numeric(x)))
    ) %>%
    group_by(.cellline, Marker) %>%
    summarise(
      n_finite = sum(is.finite(x_ln)),
      thr_ln = {
        if (n_finite[1] < min_n) NA_real_
        else fit_gmm_threshold(x_ln, marker_name = paste(.cellline[1], Marker[1]))$thr_ln
      },
      .groups = "drop"
    ) %>%
    mutate(
      thr_scaled = ifelse(is.na(thr_ln), NA_real_, exp(thr_ln))
    ) %>%
    rename(CellLine = .cellline)
}

## ------------------------------------------------------------
## Assign GMM calls (prevents your length-0 assignment crash)
## ------------------------------------------------------------

# Adds logical columns: <marker>_GMM_high
# Uses scaled if available, else ln (so it never tries to compare NULL > threshold)
assign_global_gmm_calls <- function(df,
                                    markers,
                                    out_col_suffix = "_GMM_high",
                                    thr_ln_overrides = list(),
                                    verbose = TRUE) {
  
  validate_marker_columns(df, markers, require = "either", verbose = verbose)
  
  for (mk in markers) {
    # Fit on ln-scale (or use override)
    x_ln <- get_marker_ln(df, mk)
    
    if (!is.null(thr_ln_overrides[[mk]])) {
      thr_ln <- as.numeric(thr_ln_overrides[[mk]])
      if (!is.finite(thr_ln)) stop("thr_ln_overrides for '", mk, "' is not finite.")
      fit <- fit_gmm_threshold(x_ln, marker_name = mk)
      fit$thr_ln <- thr_ln
      fit$thr_scaled <- exp(thr_ln)
    } else {
      fit <- fit_gmm_threshold(x_ln, marker_name = mk)
    }
    
    mk_scaled <- paste0(mk, "_scaled")
    mk_ln     <- paste0(mk, "_ln")
    
    out_col <- paste0(mk, out_col_suffix)
    
    if (mk_scaled %in% names(df)) {
      df[[out_col]] <- df[[mk_scaled]] > fit$thr_scaled
    } else if (mk_ln %in% names(df)) {
      df[[out_col]] <- df[[mk_ln]] > fit$thr_ln
    } else {
      stop("Unexpected: marker '", mk, "' missing both ln and scaled after validation.")
    }
    
    if (verbose) {
      message(sprintf("[GMM] Assigned %s using %s", out_col,
                      if (mk_scaled %in% names(df)) mk_scaled else mk_ln))
    }
  }
  
  df
}

## ------------------------------------------------------------
## Runner: global fits + mixture plots + optional call assignment
## ------------------------------------------------------------

run_global_gmm_with_plots <- function(df,
                                      markers,
                                      out_dir,
                                      file_prefix = "",
                                      plot_width = 2,
                                      plot_height = 2,
                                      assign_calls = TRUE,
                                      out_col_suffix = "_GMM_high",
                                      thr_ln_overrides = list(),
                                      verbose = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Global threshold table + fits/plots
  global_tbl <- compute_global_thresholds(df, markers, verbose = verbose)
  
  fits  <- list()
  plots <- list()
  
  for (mk in markers) {
    x_ln <- get_marker_ln(df, mk)
    fit  <- fit_gmm_threshold(x_ln, marker_name = mk)
    
    # Apply override (for consistent reproduction)
    if (!is.null(thr_ln_overrides[[mk]])) {
      thr_ln <- as.numeric(thr_ln_overrides[[mk]])
      if (!is.finite(thr_ln)) stop("thr_ln_overrides for '", mk, "' is not finite.")
      fit$thr_ln <- thr_ln
      fit$thr_scaled <- exp(thr_ln)
    }
    
    fits[[mk]] <- fit
    
    p <- plot_gmm_mixture_1d(
      x_ln  = x_ln,
      fit   = fit,
      xlab  = paste0(mk, " ln(scaled intensity)"),
      title = ""
    )
    plots[[mk]] <- p
    
    ggsave(
      filename = file.path(out_dir, paste0(file_prefix, mk, "_mixture.pdf")),
      plot = p,
      width = plot_width,
      height = plot_height,
      units = "in",
      useDingbats = FALSE
    )
  }
  
  # Optionally assign per-cell calls
  df_out <- df
  if (assign_calls) {
    df_out <- assign_global_gmm_calls(
      df = df_out,
      markers = markers,
      out_col_suffix = out_col_suffix,
      thr_ln_overrides = thr_ln_overrides,
      verbose = verbose
    )
  }
  
  list(
    df = df_out,
    global_thresholds = global_tbl,
    fits = fits,
    plots = plots
  )
}
