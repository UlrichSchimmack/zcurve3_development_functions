plot.zcurve <- function(
  x,
  histogram_col = "blue3",
  curve_col = "red",
  histogram_width = 0.2,
  local_interval = 0.5,
  xlim = c(0, 6),
  ylim = c(0, 0.6),
  ylab = "Density",
  show_text = TRUE,
  show_local_power = TRUE,
  line_width = 4,
  text_size = 1,
  show_ci_band = TRUE,        # pointwise bootstrap density band
  show_edr_band = FALSE,      # opt-in: fitted curves at the EDR 2.5/97.5% reps
  ci_level = 0.95,
  ci_band_col = "grey70",
  edr_band_col = "darkorange",
  ...
) {
  values <- x$val_input
  if (is.null(values) && !is.null(x$plot_data$values)) {
    values <- x$plot_data$values
  }
  if (is.null(values)) {
    stop("The zcurve object does not contain val_input.")
  }

  point_values <- function(object) {
    if (is.matrix(object) || is.data.frame(object)) {
      return(as.numeric(object[1, , drop = TRUE]))
    }
    as.numeric(object)
  }

  weights <- point_values(x$w_all)
  if (!length(weights) || !any(is.finite(weights))) {
    stop("The zcurve object does not contain finite w_all estimates.")
  }
  weights[!is.finite(weights)] <- 0
  weights <- weights / sum(weights)

  ncp <- if (is.null(x$ncp)) seq_along(weights) - 1 else point_values(x$ncp)
  z_sd <- if (is.null(x$z_sd)) rep(1, length(weights)) else point_values(x$z_sd)
  if (length(ncp) != length(weights) || length(z_sd) != length(weights)) {
    stop("w_all, ncp, and z_sd must have the same number of components.")
  }

  critical_value <- if (!is.null(x$plot_data$critical_value)) {
    as.numeric(x$plot_data$critical_value)
  } else {
    stats::qnorm(0.975)
  }
  selection_start <- if (!is.null(x$plot_data$selection_interval[["lower"]])) {
    as.numeric(x$plot_data$selection_interval[["lower"]])
  } else {
    critical_value
  }

  values <- values[is.finite(values)]
  if (!length(values)) stop("There are no finite values to plot.")

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mar = c(6, 4, 2, 1),
    family = "sans",
    cex.axis = 1,
    cex.lab = 1,
    font.axis = 1,
    font.lab = 1,
    xpd = FALSE
  )

  breaks <- seq(xlim[1], xlim[2], by = histogram_width)
  if (tail(breaks, 1) < xlim[2]) breaks <- c(breaks, xlim[2])

  # The small shift reproduces the alignment of the original zcurve3 plot,
  # where the significant histogram begins at 2 when the critical value is 1.96.
  visible_values <- values[values > xlim[1] & values < xlim[2] - 0.04] + 0.04
  histogram <- graphics::hist(
    visible_values,
    breaks = breaks,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )

  graphics::plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = "",
    ylab = ylab,
    axes = FALSE,
    ...
  )

  significant_bin <- histogram$mids >= round(selection_start, 1)
  bar_col <- ifelse(
    significant_bin,
    grDevices::adjustcolor(histogram_col, alpha.f = 0.30),
    grDevices::adjustcolor(histogram_col, alpha.f = 0.20)
  )
  graphics::rect(
    histogram$breaks[-length(histogram$breaks)],
    0,
    histogram$breaks[-1],
    histogram$density,
    col = bar_col,
    border = "white"
  )

  graphics::axis(2)
  graphics::axis(1, line = 2)
  graphics::abline(h = 0)

  grid_step <- 0.01
  density_grid <- seq(xlim[1], xlim[2], by = grid_step)

  visible <- values > xlim[1] & values < xlim[2]
  observed_selected_mass <- if (any(visible)) {
    mean(values[visible] > selection_start)
  } else {
    NA_real_
  }

  # Fold the reconstructed density only for a nondirectional (folded) analysis.
  # A directional fit has signed components and must NOT be mirrored about 0,
  # or the curve is forced symmetric and sits over the wrong center.
  folded <- if (is.null(x$plot_data$folded)) TRUE else isTRUE(x$plot_data$folded)

  # Reconstruct one fitted density from a set of component parameters and
  # scale it to the observed selected mass. This is the exact recipe the
  # point estimate uses, reused per bootstrap replicate so every band curve
  # shares the point estimate's vertical scale.
  reconstruct_scaled <- function(comp_w, comp_ncp, comp_sd) {
    comp_w <- comp_w / sum(comp_w)
    cd <- vapply(
      seq_along(comp_w),
      function(j) {
        d <- stats::dnorm(density_grid, mean = comp_ncp[j], sd = comp_sd[j])
        if (folded) {
          d <- d + stats::dnorm(-density_grid, mean = comp_ncp[j], sd = comp_sd[j])
        }
        d
      },
      numeric(length(density_grid))
    )
    cd <- t(cd)
    cd <- cd / (rowSums(cd) * grid_step)
    dens <- colSums(cd * comp_w)
    model_mass <- sum(dens[density_grid >= selection_start] * grid_step)
    if (is.finite(observed_selected_mass) && is.finite(model_mass) &&
        model_mass > 0) {
      dens <- dens * (observed_selected_mass / model_mass)
    }
    dens
  }

  fitted_density <- reconstruct_scaled(weights, ncp, z_sd)

  # ---- Bootstrap density bands (drawn behind the point-estimate curve) ----
  # Two views of curve uncertainty, in different colours so their difference
  # is visible: (1) a POINTWISE band -- at each z the 2.5/97.5% quantile of the
  # reconstructed density across replicates; (2) the EDR-EXTREME fits -- the two
  # replicate curves whose EDR equals the EDR CI bounds (real, internally
  # consistent fits, but only two of them).
  drew_ci_band <- FALSE
  drew_edr_band <- FALSE
  boot_curves <- x$plot_data$boot_curves
  if (!is.null(boot_curves) && (isTRUE(show_ci_band) || isTRUE(show_edr_band))) {
    to_mat <- function(z) {
      if (is.null(dim(z))) matrix(z, nrow = length(weights)) else z
    }
    boot_w <- to_mat(boot_curves$w_all)
    boot_ncp <- to_mat(boot_curves$ncp)
    boot_sd <- to_mat(boot_curves$z_sd)
    n_boot <- ncol(boot_w)
    if (!is.null(n_boot) && n_boot >= 2L) {
      boot_dens <- vapply(
        seq_len(n_boot),
        function(b) reconstruct_scaled(boot_w[, b], boot_ncp[, b], boot_sd[, b]),
        numeric(length(density_grid))
      )
      finite_col <- apply(boot_dens, 2, function(cc) all(is.finite(cc)))
      boot_dens <- boot_dens[, finite_col, drop = FALSE]
      boot_edr <- boot_curves$EDR[finite_col]
      probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)

      if (isTRUE(show_ci_band) && ncol(boot_dens) >= 2L) {
        band <- apply(boot_dens, 1, stats::quantile, probs = probs, na.rm = TRUE)
        graphics::polygon(
          c(density_grid, rev(density_grid)),
          c(band[1, ], rev(band[2, ])),
          col = grDevices::adjustcolor(ci_band_col, alpha.f = 0.5),
          border = NA
        )
        drew_ci_band <- TRUE
      }

      if (isTRUE(show_edr_band) && ncol(boot_dens) >= 2L &&
          any(is.finite(boot_edr))) {
        i_lo <- which.min(abs(boot_edr - stats::quantile(boot_edr, probs[1], na.rm = TRUE)))
        i_hi <- which.min(abs(boot_edr - stats::quantile(boot_edr, probs[2], na.rm = TRUE)))
        graphics::lines(density_grid, boot_dens[, i_lo],
                        col = edr_band_col, lwd = 2, lty = 2)
        graphics::lines(density_grid, boot_dens[, i_hi],
                        col = edr_band_col, lwd = 2, lty = 2)
        drew_edr_band <- TRUE
      }
    }
  }

  below <- density_grid < selection_start
  selected <- density_grid >= selection_start
  graphics::lines(
    density_grid[below],
    fitted_density[below],
    col = curve_col,
    lty = 3,
    lwd = line_width
  )
  graphics::lines(
    density_grid[selected],
    fitted_density[selected],
    col = curve_col,
    lty = 1,
    lwd = line_width
  )

  critical_index <- which.min(abs(density_grid - critical_value))
  graphics::segments(
    critical_value,
    0,
    critical_value,
    fitted_density[critical_index],
    col = curve_col,
    lwd = 3
  )

  # Local power must condition on the latent NCP, not merely on membership in
  # a broad marginal-Z component. Expand continuous components before forming
  # the posterior mixture at each plotted z value.
  if (any(z_sd > 1 + 1e-10)) {
    if (!exists("expand_ncp_components", mode = "function")) {
      stop("Source zcurv3_package.R before plotting continuous components.")
    }
    local_components <- expand_ncp_components(
      ncp = ncp,
      weight = weights,
      z_sd = z_sd
    )
    local_ncp <- local_components$ncp
    local_weights <- local_components$weight
  } else {
    local_ncp <- ncp
    local_weights <- weights
  }

  local_component_power <- stats::pnorm(
    -critical_value,
    mean = local_ncp,
    sd = 1
  ) + stats::pnorm(
    critical_value,
    mean = local_ncp,
    sd = 1,
    lower.tail = FALSE
  )

  local_component_density <- vapply(
    seq_along(local_weights),
    function(j) {
      stats::dnorm(density_grid, mean = local_ncp[j], sd = 1) +
        stats::dnorm(-density_grid, mean = local_ncp[j], sd = 1)
    },
    numeric(length(density_grid))
  )
  local_component_density <- t(local_component_density)

  if (isTRUE(show_local_power) && local_interval > 0) {
    local_breaks <- seq(xlim[1], xlim[2], by = local_interval)
    if (tail(local_breaks, 1) < xlim[2]) local_breaks <- c(local_breaks, xlim[2])
    local_bins <- cut(
      density_grid,
      breaks = local_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
    mixture_density <- colSums(local_component_density * local_weights)
    local_numerator <- colSums(
      local_component_density * (local_weights * local_component_power)
    )
    local_at_grid <- local_numerator / mixture_density
    local_power <- vapply(
      levels(local_bins),
      function(bin) {
        index <- which(local_bins == bin)
        sum(local_at_grid[index] * mixture_density[index]) /
          sum(mixture_density[index])
      },
      numeric(1)
    )
    local_midpoints <- head(local_breaks, -1) + diff(local_breaks) / 2
    graphics::mtext(
      paste0(round(100 * local_power), "%"),
      side = 1,
      line = 0.5,
      at = local_midpoints,
      cex = 1,
      las = 1
    )
  }

  estimate_percent <- function(values) {
    estimate <- suppressWarnings(as.numeric(values[1]))
    if (!is.finite(estimate)) return("NA")
    out <- paste0(round(100 * estimate), " %")
    # Append the bootstrap CI in brackets when it is available (lb/ub are
    # NA for boot_iter = 0, in which case only the point estimate is shown).
    if (length(values) >= 3) {
      lb <- suppressWarnings(as.numeric(values[2]))
      ub <- suppressWarnings(as.numeric(values[3]))
      if (is.finite(lb) && is.finite(ub)) {
        out <- paste0(out, " [", round(100 * lb), ", ", round(100 * ub), "]")
      }
    }
    out
  }

  if (isTRUE(show_text)) {
    graphics::par(xpd = NA)
    ymax <- ylim[2]
    line_height <- ymax * 0.03
    left_x <- xlim[1]
    right_x <- xlim[2]

    version <- x$metadata$version
    if (is.null(version) || !nzchar(version)) version <- ""
    version <- sub("^[Vv]", "", version)
    version_text <- if (nzchar(version)) paste0("zcurve(", version, ")") else "zcurve3"
    date_text <- x$metadata$date
    if (is.null(date_text)) date_text <- ""

    line <- 0
    graphics::text(
      left_x,
      ymax - line_height * line,
      "Schimmack, Bartos, & Brunner",
      pos = 4,
      cex = text_size
    )
    line <- line + 1.5
    graphics::text(left_x, ymax - line_height * line, version_text, pos = 4, cex = text_size)
    line <- line + 1.5
    graphics::text(left_x, ymax - line_height * line, date_text, pos = 4, cex = text_size)

    n_tests <- length(values)
    n_significant <- sum(values > critical_value)
    n_not_shown <- sum(values > xlim[2])

    line <- 0
    n_clusters <- x$metadata$n_clusters
    if (!is.null(n_clusters) && is.finite(n_clusters)) {
      graphics::text(
        right_x,
        ymax - line_height * line,
        paste0(format(n_clusters, big.mark = ","), " clusters"),
        pos = 2,
        cex = text_size
      )
      line <- line + 1.5
    }
    graphics::text(
      right_x,
      ymax - line_height * line,
      paste0(format(n_tests, big.mark = ","), " tests"),
      pos = 2,
      cex = text_size
    )
    line <- line + 1.5
    graphics::text(
      right_x,
      ymax - line_height * line,
      paste0(format(n_significant, big.mark = ","), " significant"),
      pos = 2,
      cex = text_size
    )
    line <- line + 1.5
    graphics::text(
      right_x,
      ymax - line_height * line,
      paste0(n_not_shown, " z > ", xlim[2]),
      pos = 2,
      cex = text_size
    )

    line <- line + 3

    # Aligned estimate block. %3d space-padding only lines up in a fixed-width
    # font, so this block is drawn in "mono": labels, point estimates, and CI
    # brackets fall into neat columns. Right-anchored at right_x (pos = 2), so
    # the closing "]" of every row aligns with the counts above.
    fmt_est <- function(label, v, with_ci) {
      pe <- suppressWarnings(as.numeric(v[1]))
      if (!is.finite(pe)) return(sprintf("%s   NA", label))
      if (with_ci) {
        lb <- suppressWarnings(as.numeric(v[2]))
        ub <- suppressWarnings(as.numeric(v[3]))
        if (is.finite(lb) && is.finite(ub)) {
          return(sprintf(
            "%s %3d %% [%3d, %3d]", label,
            round(100 * pe), round(100 * lb), round(100 * ub)
          ))
        }
      }
      sprintf("%s %3d %%", label, round(100 * pe))
    }

    est_labels <- c("ODR:", "EDR:", "ERR:")
    est_values <- list(x$ODR, x$EDR, x$ERR)
    if (!is.null(x$FDR)) {
      est_labels <- c(est_labels, "FDR:")
      est_values <- c(est_values, list(x$FDR))
    }

    has_ci <- length(x$EDR) >= 3 &&
      is.finite(suppressWarnings(as.numeric(x$EDR[2]))) &&
      is.finite(suppressWarnings(as.numeric(x$EDR[3])))

    graphics::par(family = "mono")
    if (has_ci) {
      conf <- x$metadata$confidence_level
      conf_pct <- if (is.numeric(conf) && length(conf) && is.finite(conf)) {
        round(100 * conf)
      } else {
        95
      }
      # Center the header over the CI bracket columns (not the full row). The
      # bracket "[%3d, %3d]" is fixed-width in the monospace font and right-
      # anchored at right_x, so its center is right_x - bracket_width/2.
      bracket_w <- graphics::strwidth("[100, 100]", cex = text_size)
      graphics::text(
        right_x - bracket_w / 2, ymax - line_height * line,
        paste0(conf_pct, "% CI"), cex = text_size
      )
      line <- line + 1.5
    }
    for (r in seq_along(est_labels)) {
      graphics::text(
        right_x, ymax - line_height * line,
        fmt_est(est_labels[r], est_values[[r]], has_ci),
        pos = 2, cex = text_size
      )
      line <- line + 2
    }
    graphics::par(family = "sans")
  }

  # Legend distinguishing the point-estimate curve from the two band views.
  if (drew_ci_band || drew_edr_band) {
    items <- "fitted (PE)"
    cols  <- curve_col
    ltys  <- 1
    lwds  <- line_width
    if (drew_ci_band) {
      items <- c(items, paste0(round(100 * ci_level), "% density band"))
      cols  <- c(cols, ci_band_col)
      ltys  <- c(ltys, 1)
      lwds  <- c(lwds, 6)
    }
    if (drew_edr_band) {
      items <- c(items, "EDR 2.5/97.5% fit")
      cols  <- c(cols, edr_band_col)
      ltys  <- c(ltys, 2)
      lwds  <- c(lwds, 2)
    }
    graphics::legend(
      "top", legend = items, col = cols, lty = ltys, lwd = lwds,
      bty = "n", cex = 0.75 * text_size, seg.len = 2
    )
  }

  invisible(x)
}



###################################################################
###################################################################
###################################################################



zcurve_forest <- function(
  object,
  min_effect = NULL,
  z_min = 4,
  es_lb_min = .2,
  es_lb_max = Inf,
  k_max = 30,
  sort_by = c(
    "lower_bound", "prob_meaningful", "adjusted_effect", "z", "study_id"
  ),
  decreasing = TRUE,
  show_raw = TRUE,
  xlim = NULL,
  digits = 2,
  point_col = "blue3",
  raw_col = "grey60",
  min_effect_col = "red3",
  labels = NULL,
  label_chars = 50,
  show_id = TRUE,
  stat = c("power", "z"),
  label_max_frac = 0.5,
  forest_frac = 0.26,
  text_size = 1,
  panel_widths = c(1.9, 5.7, 1.9),
  ...
) {
  if (!inherits(object, "zcurve")) {
    stop("object must inherit from class 'zcurve'.")
  }
  if (is.null(object$individual_effects)) {
    stop("Individual effect estimates are unavailable. Supply yi and sei to zcurve().")
  }
  if (all(is.na(object$individual_effects$effect_adjusted_lb))) {
    stop("The forest plot requires a bootstrapped fit: run zcurve() with boot_iter > 0 so the adjusted-effect intervals reflect the uncertainty in the fitted population, not just each study's own standard error.")
  }

  sort_by <- match.arg(sort_by)

  if (is.null(min_effect)) min_effect <- object$metadata$min_effect
  if (length(min_effect) != 1L || !is.finite(min_effect) || min_effect < 0) {
    stop("min_effect must be one nonnegative finite number.")
  }
  if (length(z_min) != 1L || !is.finite(z_min) || z_min < 0) {
    stop("z_min must be one nonnegative finite number.")
  }
  if (length(text_size) != 1L || !is.finite(text_size) || text_size <= 0) {
    stop("text_size must be one positive number.")
  }
  if (length(panel_widths) != 3L || any(!is.finite(panel_widths)) ||
      any(panel_widths <= 0)) {
    stop("panel_widths must be three positive numbers (labels, forest, estimates).")
  }
  if (length(k_max) != 1L || !is.finite(k_max) || k_max < 1) {
    stop("k_max must be one positive number.")
  }
  if (length(es_lb_min) != 1L || !is.finite(es_lb_min)) {
    stop("es_lb_min must be one finite number.")
  }
  if (length(es_lb_max) != 1L || is.na(es_lb_max) || es_lb_max < es_lb_min) {
    stop("es_lb_max must be one number >= es_lb_min (Inf to disable the upper cap).")
  }
  if (!is.null(label_chars) &&
      (length(label_chars) != 1L || is.na(label_chars) || label_chars < 1)) {
    stop("label_chars must be NULL (full label) or a positive number.")
  }
  if (length(label_max_frac) != 1L || !is.finite(label_max_frac) ||
      label_max_frac <= 0 || label_max_frac >= 1) {
    stop("label_max_frac must be one number in (0, 1).")
  }
  if (length(forest_frac) != 1L || !is.finite(forest_frac) ||
      forest_frac <= 0 || forest_frac >= 1) {
    stop("forest_frac must be one number in (0, 1).")
  }
  stat <- match.arg(stat)
  if (identical(stat, "power") && is.null(object$individual_effects$power_adjusted)) {
    stat <- "z"  # older fit without adjusted power; fall back to observed z
  }

  forest_data <- object$individual_effects

  # Recalculate the probability when the plotting threshold differs from the
  # threshold used in the original fit. The posterior effect distribution does
  # not need to be refitted.
  fitted_min_effect <- unique(forest_data$min_effect)
  if (length(fitted_min_effect) != 1L ||
      !isTRUE(all.equal(min_effect, fitted_min_effect))) {
    posterior <- object$individual_posterior
    if (is.null(posterior)) {
      stop("The fitted object does not contain individual posterior probabilities.")
    }
    theta <- outer(posterior$sei, posterior$component_ncp, "*")
    forest_data$prob_meaningful <- rowSums(
      posterior$probabilities * (abs(theta) > min_effect),
      na.rm = FALSE
    )
    forest_data$min_effect <- min_effect
  }

  # Row labels. A user-supplied vector replaces the study_id / row-number
  # column: either named by study_id, or one label per study in fit order.
  if (is.null(labels)) {
    forest_data$label <- as.character(forest_data$study_id)
  } else if (!is.null(names(labels))) {
    matched <- as.character(labels)[
      match(as.character(forest_data$study_id), names(labels))
    ]
    forest_data$label <- ifelse(
      is.na(matched), as.character(forest_data$study_id), matched
    )
  } else if (length(labels) == nrow(forest_data)) {
    forest_data$label <- as.character(labels)
  } else {
    stop(
      "labels must be named by study_id or have one entry per study (length ",
      nrow(forest_data), ", the number of studies fit with yi/sei)."
    )
  }


  forest_data <- forest_data[
    is.finite(forest_data$z) & abs(forest_data$z) > z_min,
    ,
    drop = FALSE
  ]
  if (!nrow(forest_data)) {
    stop("No studies exceed z_min.")
  }

  # Keep studies whose adjusted lower bound falls in the window
  # [es_lb_min, es_lb_max]. Default es_lb_max = Inf leaves it a simple floor.
  keep <- is.finite(forest_data$effect_adjusted_lb) &
    forest_data$effect_adjusted_lb >= es_lb_min &
    forest_data$effect_adjusted_lb <= es_lb_max
  forest_data <- forest_data[keep, , drop = FALSE]
  if (!nrow(forest_data)) {
    stop("No studies fall in the lower-bound window [es_lb_min, es_lb_max].")
  }

  order_value <- switch(
    sort_by,
    lower_bound = forest_data$effect_adjusted_lb,
    prob_meaningful = forest_data$prob_meaningful,
    adjusted_effect = abs(forest_data$effect_adjusted),
    z = abs(forest_data$z),
    study_id = seq_len(nrow(forest_data))
  )
  forest_data <- forest_data[
    order(order_value, decreasing = decreasing, na.last = TRUE),
    ,
    drop = FALSE
  ]
  rownames(forest_data) <- NULL

  # Cap to at most k_max studies AFTER sorting, so it's the top k_max by sort_by.
  if (nrow(forest_data) > k_max) {
    forest_data <- forest_data[seq_len(k_max), , drop = FALSE]
  }

  # Shorten long labels for display: cut at the first ":" (keep the main title,
  # drop the subtitle), then keep as many WHOLE words as fit within label_chars
  # characters. Budgeting by characters rather than words means short tokens
  # like "-->" don't burn a word slot. Supply full titles / author lists as-is;
  # label_chars = NULL leaves the label untouched.
  if (!is.null(label_chars)) {
    main <- sub(":.*$", "", trimws(forest_data$label))
    forest_data$label <- vapply(main, function(s) {
      w <- strsplit(trimws(s), "\\s+")[[1]]
      w <- w[nzchar(w)]
      if (!length(w)) return("")
      out <- w[1]; len <- nchar(w[1])
      for (j in seq_along(w)[-1]) {
        nl <- len + 1L + nchar(w[j])
        if (nl > label_chars) break
        out <- c(out, w[j]); len <- nl
      }
      paste(out, collapse = " ")
    }, character(1), USE.NAMES = FALSE)
  }

  interval_text <- sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    forest_data$effect_adjusted,
    forest_data$effect_adjusted_lb,
    forest_data$effect_adjusted_ub
  )
  forest_data$effect_adjusted_interval <- interval_text

  if (is.null(xlim)) {
    finite_limits <- c(
      0, min_effect, forest_data$effect_adjusted_lb,
      forest_data$effect_adjusted_ub,
      if (isTRUE(show_raw)) forest_data$effect_observed else NULL
    )
    finite_limits <- finite_limits[is.finite(finite_limits)]
    padding <- max(diff(range(finite_limits)), .2) * .06
    xlim <- range(finite_limits) + c(-padding, padding)
    if (object$metadata$effect_scale == "absolute") xlim[1] <- min(0, xlim[1])
  }

  n_studies <- nrow(forest_data)
  y <- rev(seq_len(n_studies))
  row_cex <- max(.55, min(.95, 22 / n_studies)) * text_size
  header_cex <- .85 * text_size
  point_cex <- .85 * text_size
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(old_par)
  }, add = TRUE)

  # ---- Column layout: id | Study labels | forest | estimates -------------
  # The id, label and estimate columns hold text and are sized in absolute cm
  # (lcm) to fit their content; the forest is the single RELATIVE column and
  # fills whatever is left, so the text columns can never squeeze it out. The
  # label column is additionally capped at label_max_frac of the figure width;
  # when a label would exceed the cap it is re-fit to whole words within it.
  strw_in <- function(s, cex) {
    tryCatch(max(graphics::strwidth(s, units = "inches", cex = cex), na.rm = TRUE),
             error = function(e) NA_real_)
  }
  fig_w <- tryCatch(graphics::par("din")[1], error = function(e) NA_real_)

  # estimates column content: the chosen statistic (adjusted power or z) + CI
  if (stat == "power") {
    stat_txt  <- ifelse(is.finite(forest_data$power_adjusted),
                        formatC(100 * forest_data$power_adjusted, format = "f", digits = 0), "")
    stat_head <- "Power"
  } else {
    stat_txt  <- formatC(abs(forest_data$z), format = "f", digits = 2)
    stat_head <- "z"
  }
  eff_head <- paste0("Effect [", round(100 * object$metadata$confidence_level), "% CI]")

  show_id_col <- isTRUE(show_id) && !is.null(forest_data$study_id) &&
    any(nzchar(as.character(forest_data$study_id)))
  id_txt <- as.character(forest_data$study_id)
  id_cm  <- if (show_id_col) {
    iw <- strw_in(id_txt, row_cex); if (is.finite(iw)) (iw + 0.08) * 2.54 else NA_real_
  } else NA_real_

  # Estimates column: statistic (a_est wide) + a small gap + effect [CI] (b_est
  # wide), both right-aligned. The stat sub-column is sized to the VALUES only,
  # not the header: a wide header like "Power" (vs "100") is right-aligned to the
  # same edge and overhangs left into the empty strip above the rows, so no width
  # is wasted on either side of the numbers. gap_est is the ONLY space between
  # the statistic and the effect.
  a_est <- strw_in(stat_txt, row_cex)
  if (!is.finite(a_est)) a_est <- 0
  b_est <- max(strw_in(interval_text, row_cex), strw_in(eff_head, header_cex), na.rm = TRUE)
  if (!is.finite(b_est)) b_est <- 0
  gap_est <- 0.10
  est_in  <- a_est + gap_est + b_est + 0.06
  est_cm  <- est_in * 2.54
  x_stat  <- 1 - (b_est + gap_est) / est_in   # right edge of the stat sub-column

  # Study-label column. Labels are capped so the forest keeps at least
  # forest_frac of the figure width: whatever the auto-sized id and estimate
  # columns leave is split between labels and forest, forest floored first.
  # Short labels leave the forest even wider (it is the relative column), so the
  # cap only bites when titles are genuinely long.
  label_in <- max(strw_in(forest_data$label, row_cex), strw_in("Study", header_cex), na.rm = TRUE)
  id_in    <- if (is.finite(id_cm)) id_cm / 2.54 else 0
  cap_in   <- if (is.finite(fig_w)) {
    room <- fig_w - id_in - est_in - forest_frac * fig_w
    max(0.6, min(label_max_frac * fig_w, room))
  } else NA_real_
  if (is.finite(label_in) && is.finite(cap_in) && label_in > cap_in) {
    fit_w <- function(s) {
      wv <- strsplit(trimws(s), "\\s+")[[1]]; wv <- wv[nzchar(wv)]
      if (!length(wv)) return("")
      keep <- min(2L, length(wv))          # never emit a one-word stub ("The")
      out  <- wv[seq_len(keep)]
      for (j in seq_along(wv)[-seq_len(keep)]) {
        cand <- paste(c(out, wv[j]), collapse = " ")
        if (strw_in(cand, row_cex) > cap_in - 0.10) break
        out <- c(out, wv[j])
      }
      paste(out, collapse = " ")
    }
    forest_data$label <- vapply(forest_data$label, fit_w, character(1), USE.NAMES = FALSE)
    label_in <- max(strw_in(forest_data$label, row_cex), strw_in("Study", header_cex), na.rm = TRUE)
  }
  label_cm <- if (is.finite(label_in)) (label_in + 0.10) * 2.54 else NA_real_

  use_id <- show_id_col && is.finite(id_cm) && is.finite(label_cm)
  lab_w  <- if (is.finite(label_cm)) graphics::lcm(label_cm) else panel_widths[1]
  est_w  <- if (is.finite(est_cm))   graphics::lcm(est_cm)   else panel_widths[3]

  if (use_id) {
    graphics::layout(matrix(1:4, nrow = 1),
                     widths = c(graphics::lcm(id_cm), lab_w, 1, est_w))
    graphics::par(mar = c(5, 0, 3, 0), xpd = NA)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(.5, n_studies + 1.2))
    graphics::text(1, y, id_txt, pos = 2, cex = row_cex)
  } else {
    graphics::layout(matrix(1:3, nrow = 1), widths = c(lab_w, 1, est_w))
  }

  # Study-label column. Bottom margin 5 matches the forest/estimates panels so
  # every column maps the shared ylim onto the same height (row alignment).
  graphics::par(mar = c(5, if (use_id) 0.2 else 0, 3, 0.2), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(.5, n_studies + 1.2))
  graphics::text(0, n_studies + .8, "Study \n", pos = 4, font = 2, cex = header_cex)
  graphics::text(0, y, forest_data$label, pos = 4, cex = row_cex)

  graphics::par(mar = c(5, 0.6, 3, 0.3), xpd = FALSE)
  graphics::plot(
    NA_real_, NA_real_,
    xlim = xlim,
    ylim = c(.5, n_studies + 1.2),
    xlab = if (object$metadata$effect_scale == "absolute") {
      "Adjusted absolute effect size |d|"
    } else {
      "Adjusted effect size"
    },
    ylab = "",
    yaxt = "n",
    bty = "n",
    ...
  )

  graphics::abline(v = 0, col = "grey35", lwd = 1.5)
  graphics::abline(v = min_effect, col = min_effect_col, lty = 2, lwd = 1.5)
  graphics::segments(
    forest_data$effect_adjusted_lb, y,
    forest_data$effect_adjusted_ub, y,
    col = point_col, lwd = 1.5
  )
  if (isTRUE(show_raw)) {
    graphics::points(
      forest_data$effect_observed, y,
      pch = 1, col = raw_col, cex = point_cex
    )
  }
  graphics::points(
    forest_data$effect_adjusted, y,
    pch = 16, col = point_col, cex = point_cex * .9
  )
  graphics::mtext(
    paste0("Minimum effect = ", format(min_effect, trim = TRUE)),
    side = 3, line = .4, at = min_effect, adj = .5,
    col = min_effect_col, cex = .5 * text_size
  )

  # Estimates column: the statistic (adjusted Power or observed z) on the left,
  # the adjusted effect + CI on the right, both right-aligned within the
  # auto-sized column (est_cm = stat + gap + effect).
  graphics::par(mar = c(5, 0.15, 3, 0.15), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(.5, n_studies + 1.2))
  graphics::text(x_stat, n_studies + .8, paste0(stat_head, " \n"), pos = 2, font = 2, cex = header_cex)
  graphics::text(1, n_studies + .8, paste0(eff_head, " \n"), pos = 2, font = 2, cex = header_cex)
  graphics::text(x_stat, y, stat_txt, pos = 2, cex = row_cex)
  graphics::text(1, y, interval_text, pos = 2, cex = row_cex)

  invisible(forest_data)
}


forest.zcurve <- function(x, ...) {
  zcurve_forest(x, ...)
}
