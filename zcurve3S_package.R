

#rm(list = ls())
#options(scipen = 999)  # Disable scientific notation




### INSTALL PACKAGES (only once – manually run if needed)
if (1 == 2) {  # This block is ignored unless manually changed to (1 == 1)
  install.packages("KernSmooth")
  install.packages("parallel")
} # END install block

### LOAD LIBRARIES
library(parallel)
library(KernSmooth)



########################################################################
### DRAFT: new zcurve() package interface + main function
### NOT YET intEGRATED into zcurve.V3.90.R -- for review/testing only.
########################################################################



#######################################################################
### FITTING CONTROL
#######################################################################

zcurve_control <- function(
  boot_iter            = 0,       # bootstrap iterations (500+ recommended for final models)
  parallel             = TRUE,    # parallelize bootstrap/EM where available
  cores                = NULL,    # NULL chooses a conservative automatic worker count
  seed                 = NULL,    # optional reproducible bootstrap seed

  tol_criterion        = 1e-6,    # EM convergence threshold
  max_iter             = 2000,    # max EM iterations
  max_iter_boot        = 500,     # max EM iterations per bootstrap replicate

  ci_alpha             = .05,     # CI level for reported intervals (default 95%)

  augment              = TRUE,    # bias correction at the lower bound (OF method)
  n_bars               = 512,     # density-estimation resolution (OF method)
  density_step         = 0.01,    # numerical density-grid step for the point estimate
  bootstrap_density_step = 0.05,  # coarser numerical grid for faster bootstrap fits
  bw_est               = 0.05,    # kernel bandwidth used for fitting (lower = less smoothing)
  spline_width         = .5,      # blending of truncated and normal KD (OF method)

  floor_correct        = TRUE,    # remove the homogeneous-null floor from es_tau_sig
  floor_reps           = 5,       # homogeneous-null refits averaged to estimate that floor

  es_null_share        = 0.5      # INDIVIDUAL-effects prior only: give the null (ncp 0)
                                  # at least this share of the low-power (ncp<=1) mass, so
                                  # genuinely-null studies can shrink to 0. Applied to w_all
                                  # for the individual posterior only; the fit's w_inp and
                                  # thus EDR/ERR/FDR are untouched. 0 disables.
) {
  if (length(density_step) != 1L || !is.finite(density_step) || density_step <= 0) {
    stop("density_step must be one positive finite number.")
  }
  if (length(bootstrap_density_step) != 1L ||
      !is.finite(bootstrap_density_step) || bootstrap_density_step <= 0) {
    stop("bootstrap_density_step must be one positive finite number.")
  }
  if (length(ci_alpha) != 1L || !is.finite(ci_alpha) ||
      ci_alpha <= 0 || ci_alpha >= 1) {
    stop("ci_alpha must be one number strictly between 0 and 1.")
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.")
  }
  if (!is.null(cores) &&
      (length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != as.integer(cores))) {
    stop("cores must be NULL or one positive integer.")
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
       seed > .Machine$integer.max || seed != as.integer(seed))) {
    stop("seed must be NULL or one nonnegative integer.")
  }
  if (!is.logical(floor_correct) || length(floor_correct) != 1L ||
      is.na(floor_correct)) {
    stop("floor_correct must be TRUE or FALSE.")
  }
  if (length(floor_reps) != 1L || !is.finite(floor_reps) ||
      floor_reps < 1 || floor_reps != as.integer(floor_reps)) {
    stop("floor_reps must be one positive integer.")
  }
  if (length(es_null_share) != 1L || !is.finite(es_null_share) ||
      es_null_share < 0 || es_null_share > 1) {
    stop("es_null_share must be one number in [0, 1].")
  }

  control <- list(
    boot_iter = boot_iter, parallel = parallel, cores = cores, seed = seed,
    tol_criterion = tol_criterion, max_iter = max_iter, max_iter_boot = max_iter_boot,
    ci_alpha = ci_alpha,
    augment = augment, n_bars = n_bars,
    density_step = density_step,
    bootstrap_density_step = bootstrap_density_step,
    bw_est = bw_est, spline_width = spline_width,
    floor_correct = floor_correct, floor_reps = as.integer(floor_reps),
    es_null_share = es_null_share
  )

  # Allows the main function to detect the legacy placement explicitly.
  attr(control, "boot_iter_supplied") <- !missing(boot_iter)
  control
} ### EOF zcurve_control


#######################################################################
### COMMON LATENT-NCP COMPONENT REPRESENTATION
#######################################################################

# Convert a fitted mixture to the common representation used by the
# downstream estimation code. At this first stage every fitted component is
# discrete, so it contributes exactly one latent-NCP point. Continuous
# components will later contribute several quadrature points through this same
# interface; the power and effect-size code will not need a separate branch.
expand_ncp_components <- function(
    ncp, weight, z_sd = NULL, quadrature_points = 31L) {
  component_count <- length(ncp)

  if (is.null(z_sd)) z_sd <- rep(1, component_count)
  if (length(z_sd) == 1L) z_sd <- rep(z_sd, component_count)

  if (length(weight) != component_count || length(z_sd) != component_count) {
    stop("ncp, weight, and z_sd must have compatible lengths.")
  }
  if (any(!is.finite(ncp)) || any(!is.finite(weight)) ||
      any(!is.finite(z_sd)) || any(weight < 0) || any(z_sd < 1 - 1e-10)) {
    stop("Component locations, weights, and SDs must be finite and valid.")
  }
  if (sum(weight) <= 0) stop("Component weights must have a positive sum.")
  if (length(quadrature_points) != 1L || !is.finite(quadrature_points) ||
      quadrature_points < 3 || quadrature_points != as.integer(quadrature_points)) {
    stop("quadrature_points must be one integer of at least 3.")
  }

  # Golub-Welsch quadrature for expectations under a standard normal. The
  # eigenvalues are the nodes and the squared first-row eigenvectors are the
  # probability weights. This avoids adding another package dependency.
  normal_quadrature <- function(point_count) {
    jacobi <- matrix(0, nrow = point_count, ncol = point_count)
    off_diagonal <- sqrt(seq_len(point_count - 1L))
    jacobi[cbind(seq_len(point_count - 1L), 2:point_count)] <- off_diagonal
    jacobi[cbind(2:point_count, seq_len(point_count - 1L))] <- off_diagonal
    eig <- eigen(jacobi, symmetric = TRUE)
    component_order <- order(eig$values)
    list(
      node = eig$values[component_order],
      weight = eig$vectors[1, component_order]^2
    )
  }

  quadrature <- normal_quadrature(as.integer(quadrature_points))
  expanded <- vector("list", component_count)

  for (component in seq_len(component_count)) {
    # z_sd is the marginal SD of the observed statistic. Sampling error has
    # unit variance, leaving this SD for the latent NCP distribution.
    ncp_sd <- sqrt(max(z_sd[component]^2 - 1, 0))

    if (ncp_sd <= 1e-10) {
      node <- 0
      node_weight <- 1
    } else {
      node <- quadrature$node
      node_weight <- quadrature$weight
    }

    expanded[[component]] <- data.frame(
      component_id = rep(component, length(node)),
      ncp = ncp[component] + ncp_sd * node,
      ncp_sd = rep(ncp_sd, length(node)),
      z_sd = rep(1, length(node)),
      component_z_sd = rep(z_sd[component], length(node)),
      component_weight = rep(weight[component], length(node)),
      quadrature_weight = node_weight,
      weight = weight[component] * node_weight,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, expanded)
}


#######################################################################
### MAIN FITTING FUNCTION
#######################################################################

zcurve <- function(
  ### DATA INPUT (supply exactly one consistent combination)
  zval         = NULL,
  tval         = NULL,
  df           = NULL,
  pval         = NULL,
  n            = NULL,
  yi           = NULL,
  sei          = NULL,
  study_id     = NULL,
  cluster_id   = NULL,

  ### CORE MODEL CHOICES
  est_method   = "OF",           # "SQP", "OF", "EM", "ML_gamma"
  directional  = FALSE,          # FALSE: sign of test_stat doesn't matter (default meta-science use)
  folded       = TRUE,           # TRUE: input is |z| or |t| (matches directional = FALSE almost always)

  alpha        = 0.05,           # significance threshold
  crit         = qnorm(1 - alpha / 2),   # two-sided critical value, derived from alpha (z-scale default)
  int_beg      = crit,                   # selection window start, derived from crit
  int_end      = 6,              # right-truncation boundary for the modeling grid (test_stat > int_end: power = 1)
  edr_ci_adj   = .05,
  err_ci_adj   = .03,
  min_effect   = .20,             # smallest effect size of practical interest

  ### MODEL STRUCTURE (component grid)
  ncp          = 0:6,             # component locations
  k_ncp        = length(ncp),     # number of components
  z_sd         = rep(1, k_ncp),   # one SD per component
  z_sd_fixed   = TRUE,            # hold component SDs fixed (EM method)
  ncp_fixed    = TRUE,            # hold component locations fixed (EM method)
  w_fixed      = FALSE,           # hold component weights fixed (EM method)

  x_lim_min    = 0,               # lower bound of the modeling grid
  x_lim_max    = 6,               # upper bound of the modeling grid (plot x-axis)
  int_loc      = 1,               # bin width for local-power/local-ES reporting (0 disables)

  ### OUTPUT / DISPLAY
  version      = "3S",            # set to NULL or "" to suppress
  date         = "2026.07.22",    # set to NULL or "" to suppress
  show_plot    = TRUE,            # draw the fitted density automatically
  show_summary = FALSE,           # print summary output immediately when TRUE
  boot_iter    = 0,               # bootstrap iterations (500+ for final models)

  ### FITTING CONTROL
  control      = zcurve_control()
) {


boot_iter_supplied <- !missing(boot_iter)
control_boot_iter_supplied <- isTRUE(attr(control, "boot_iter_supplied"))

# Recognize control objects created before the explicit-use attribute existed.
if (is.null(attr(control, "boot_iter_supplied")) &&
    !is.null(control$boot_iter) && !identical(control$boot_iter, 0)) {
  control_boot_iter_supplied <- TRUE
}

if (boot_iter_supplied && control_boot_iter_supplied) {
  stop(
    "Supply boot_iter either directly to zcurve() or inside zcurve_control(), not both."
  )
}
if (!boot_iter_supplied && control_boot_iter_supplied) {
  boot_iter <- control$boot_iter
}
if (length(boot_iter) != 1L || !is.finite(boot_iter) ||
    boot_iter < 0 || boot_iter != as.integer(boot_iter)) {
  stop("boot_iter must be one nonnegative integer.")
}
boot_iter <- as.integer(boot_iter)
control$boot_iter <- boot_iter

if (length(min_effect) != 1L || !is.finite(min_effect) || min_effect < 0) {
  stop("min_effect must be one nonnegative finite number.")
}

### --- component grid check -------------------------------------------------
# Supplying your own ncp / z_sd is fully supported (e.g. components with varying
# means and SDs). They must simply describe the SAME components. A length
# mismatch otherwise produces a silently wrong fit rather than an error.
if (length(ncp) < 1L || any(!is.finite(ncp))) {
  stop("ncp must be a finite numeric vector with at least one component.")
}
if (length(z_sd) != length(ncp)) {
  stop(sprintf(paste0(
    "z_sd must have the same length as ncp (got length(z_sd) = %d, ",
    "length(ncp) = %d).\n  Note: writing z_sd = rep(1, length(ncp)) inside the ",
    "zcurve() call evaluates length(ncp) in the CALLING environment, so it uses ",
    "your own `ncp` object rather than the ncp argument. Define ncp first, then ",
    "pass ncp = ncp, z_sd = rep(1, length(ncp))."),
    length(z_sd), length(ncp)))
}
if (any(!is.finite(z_sd)) || any(z_sd < 1 - 1e-10)) {
  stop("z_sd must be finite and >= 1 (the marginal SD of z cannot fall below the sampling SD of 1).")
}

# EXPERIMENTAL: choose the source of the significant-results ES estimates
# (es_mean_sig, es_tau_sig). FALSE = NCP mixture pathway only
# (get_es_estimates). TRUE = override with fit_folded_effect_model when any
# component SD > 1 (the previous default). Set to FALSE to study the mixture
# pathway on its own.
use_folded_es_model <- FALSE

# EXPERIMENTAL: use the study-level normal-RE fit_tau_sig() for es_tau_sig.
# FALSE = keep the distribution-free ncp mixture for es_tau_sig (the default;
# fit_tau_sig re-introduces a normal-effect assumption that fails on bimodal /
# high-FDR data, so it is off).
use_fit_tau_sig <- FALSE

# es_tau_sig floor correction (set via zcurve_control): subtract a plug-in
# homogeneous-null floor with a variance-share blend, so tau ~ 0 collapses while
# genuine (high-tau) heterogeneity is left essentially untouched. NULL-guarded so
# control objects built before these options existed still work.
floor_correct <- if (is.null(control$floor_correct)) TRUE else isTRUE(control$floor_correct)
floor_reps    <- if (is.null(control$floor_reps)) 5L else as.integer(control$floor_reps)


zcurve_call <- match.call()


### begin of actual function


##############################################
### Get Densities
##############################################

build_dens <- function(d_x, ncp, z_sd, df = NULL, 
                       curve_type = "z", folded = TRUE) {
  
  xmat  <- matrix(rep(d_x, each = length(ncp)), nrow = length(ncp))
  mumat <- matrix(rep(ncp, times = length(d_x)), nrow = length(ncp))
  
  if (curve_type == "z") {
    
    if (is.null(z_sd)) {
      z_sd <- rep(1, length(ncp))
    }
    
    if (length(z_sd) == 1) {
      z_sd <- rep(z_sd, length(ncp))
    }
    
    sdmat <- matrix(rep(z_sd, times = length(d_x)), nrow = length(ncp))
    

    if (folded) {
      
      # Standard z-curve: folded normal density.
      dens <- dnorm(xmat, mean = mumat, sd = sdmat) +
              dnorm(-xmat, mean = mumat, sd = sdmat)
      
    } else {

      # mass of component k inside [int_beg, int_end]
      denom <- pnorm(int_end, mean = mumat, sd = sdmat) -
               pnorm(int_beg, mean = mumat, sd = sdmat)
      denom <- pmax(denom, .Machine$double.xmin)

      dens <- dnorm(xmat, mean = mumat, sd = sdmat) / denom

     
    }
    
  } else {
    
    if (is.null(df)) {
      stop("df must be supplied when curve type is not 'z'.")
    }
    
    if (length(df) == 1) {
      df <- rep(df, length(ncp))
    }
    
    dfmat <- matrix(rep(df, times = length(d_x)), nrow = length(ncp))
    
    if (folded) {
      
      # Standard z-curve: folded noncentral t density.
      dens <- dt(xmat, df = dfmat, ncp = mumat) +
              dt(-xmat, df = dfmat, ncp = mumat)
      
    } else {
      
      # Directional model: signed positive noncentral t density,
      # truncated at zero.
      denom <- pt(0, df = dfmat, ncp = mumat, lower.tail = FALSE)
      denom <- pmax(denom, .Machine$double.xmin)
      
      dens <- dt(xmat, df = dfmat, ncp = mumat) / denom
      
      # Only relevant if d_x accidentally contains negative values.
      dens[xmat < 0] <- 0
    }
  }

    
  dens

}

###

#ddd

get_densities = function(
    int, bw = 0.20, from = 1.96, to = 6, width = 1,
    augment = TRUE, density_step = 0.01) {

  #from = 0; to = 6;bw = .05
  #
  #summary(int)

  n <- max(2L, length(seq(from, to, by = density_step)))
  grid <- seq(from, to, length.out = n)
  dens <- numeric(n)

  if (augment) {

    bnd_width <- width

    splice <- from + bnd_width
    bnd_idx <- which(grid < splice)
    int_idx <- which(grid >= splice)
    
    ## ---- Truncated normal kernel for boundary zone ----
    int_s <- int - from
    for (i in bnd_idx) {
      z_s <- grid[i] - from
      raw <- dnorm(int_s, z_s, bw)
      trunc_corr <- pnorm(z_s / bw)
      dens[i] <- mean(raw) / trunc_corr
    }
 
    ## ---- bkde for interior ----
    bkde_out <- bkde(int, bandwidth = bw, range.x = c(from, to), gridsize = n)
    dens[int_idx] <- bkde_out$y[int_idx]
    
    ## ---- Calibrate: scale bkde to match truncated normal at splice ----
    tn_at_splice <- dens[max(bnd_idx)]
    bkde_at_splice <- bkde_out$y[min(int_idx)]
    if (bkde_at_splice > 0) {
      scale_factor <- tn_at_splice / bkde_at_splice
      dens[int_idx] <- dens[int_idx] * scale_factor
    
    
      ## ---- Blend over transition ----
      blend_n <- min(5, length(bnd_idx))
      for (j in 1:blend_n) {
        idx <- max(bnd_idx) - blend_n + j
        w <- j / (blend_n + 1)
        dens[idx] <- (1 - w) * dens[idx] + w * bkde_out$y[idx] * scale_factor
      }
    
      dens <- pmax(0, dens)
    }

    
  } else {

    bkde_out <- bkde(int, bandwidth = bw, range.x = c(from, to), gridsize = n)
    dens <- bkde_out$y

  }
  
  dens_df <- data.frame(zx = grid, zy = dens)
  dens_bw <- dens_df$zx[2] - dens_df$zx[1]
  dens_df$zy <- dens_df$zy / sum(dens_df$zy * dens_bw)

  #plot(dens_df$zx, dens_df$zy)

  return(dens_df)

} # EOF Densities 


# Random-effects effect-size model with a fixed selection weight derived
# from ODR and EDR. This estimates the latent effect-size heterogeneity
# directly; it is not the SD of empirical-Bayes posterior means.
fit_fixed_weight_selection <- function(
    yi, sei, ODR, EDR, crit = 1.96,
    selection = c("one_sided", "two_sided"),
    init_mu = NULL, init_tau = NULL) {

  selection <- match.arg(selection)
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]
  sei <- sei[keep]
  if (length(yi) < 2L) {
    return(c(selection_weight = NA_real_, mu = NA_real_, tau = NA_real_))
  }

  eps <- 1e-10
  ODR <- pmin(pmax(ODR, eps), 1 - eps)
  EDR <- pmin(pmax(EDR, eps), ODR)
  selection_weight <- EDR * (1 - ODR) / (ODR * (1 - EDR))
  selection_weight <- pmin(pmax(selection_weight, eps), 1)

  if (is.null(init_mu)) {
    init_mu <- sum(yi / sei^2) / sum(1 / sei^2)
  }
  if (is.null(init_tau)) {
    init_tau <- stats::sd(yi)
    if (!is.finite(init_tau) || init_tau <= 0) init_tau <- 0.10
  }

  neg_log_likelihood <- function(par) {
    mu <- par[1]
    tau <- exp(par[2])
    sampling_sd <- sqrt(sei^2 + tau^2)

    if (selection == "one_sided") {
      significant <- yi / sei > crit
      p_sig <- 1 - pnorm(crit * sei, mean = mu, sd = sampling_sd)
    } else {
      significant <- abs(yi / sei) > crit
      p_sig <- pnorm(-crit * sei, mean = mu, sd = sampling_sd) +
        1 - pnorm(crit * sei, mean = mu, sd = sampling_sd)
    }

    normalizer <- pmax(p_sig + selection_weight * (1 - p_sig), eps)
    log_weight <- ifelse(significant, 0, log(selection_weight))
    log_likelihood <- dnorm(
      yi, mean = mu, sd = sampling_sd, log = TRUE
    ) + log_weight - log(normalizer)

    -sum(log_likelihood)
  }

  fit <- tryCatch(
    optim(
      par = c(init_mu, log(init_tau)),
      fn = neg_log_likelihood,
      method = "BFGS",
      control = list(maxit = 1000)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || any(!is.finite(fit$par))) {
    return(c(selection_weight = selection_weight, mu = NA_real_, tau = NA_real_))
  }

  c(
    selection_weight = selection_weight,
    mu = fit$par[1],
    tau = exp(fit$par[2])
  )
}


# Folded random-effects model on the effect-size scale. Unlike an NCP-scale
# component, this model conditions on each study's actual SE, so variation in
# precision is not mistaken for effect-size heterogeneity.
fit_folded_effect_model <- function(
    yi, sei, ODR, EDR, crit = 1.96, init_mu = NULL, init_tau = NULL) {

  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  y <- abs(yi[keep])
  sei <- sei[keep]
  if (length(y) < 2L) {
    return(c(mu = NA_real_, tau = NA_real_, mean_sig = NA_real_,
             tau_sig = NA_real_))
  }

  eps <- 1e-10
  ODR <- pmin(pmax(ODR, eps), 1 - eps)
  EDR <- pmin(pmax(EDR, eps), ODR)
  selection_weight <- EDR * (1 - ODR) / (ODR * (1 - EDR))
  selection_weight <- pmin(pmax(selection_weight, eps), 1)

  if (is.null(init_mu) || !is.finite(init_mu)) init_mu <- median(y)
  if (is.null(init_tau) || !is.finite(init_tau)) init_tau <- stats::sd(y)
  init_mu <- max(init_mu, 0)
  init_tau <- max(init_tau, 0)

  neg_log_likelihood <- function(par) {
    mu <- par[1]
    tau <- par[2]
    marginal_sd <- sqrt(sei^2 + tau^2)
    folded_density <-
      dnorm(y, mean = mu, sd = marginal_sd) +
      dnorm(-y, mean = mu, sd = marginal_sd)
    p_sig <-
      pnorm(-crit * sei, mean = mu, sd = marginal_sd) +
      1 - pnorm(crit * sei, mean = mu, sd = marginal_sd)
    normalizer <- pmax(p_sig + selection_weight * (1 - p_sig), eps)
    significant <- y / sei > crit
    log_weight <- ifelse(significant, 0, log(selection_weight))
    -sum(log(pmax(folded_density, .Machine$double.xmin)) +
           log_weight - log(normalizer))
  }

  upper_mu <- max(5, max(y) * 2)
  upper_tau <- max(5, stats::sd(y) * 5)
  fit <- tryCatch(
    optim(
      par = c(init_mu, init_tau),
      fn = neg_log_likelihood,
      method = "L-BFGS-B",
      lower = c(0, 0),
      upper = c(upper_mu, upper_tau),
      control = list(maxit = 1000)
    ),
    error = function(error) NULL
  )
  if (is.null(fit) || any(!is.finite(fit$par))) {
    return(c(mu = NA_real_, tau = NA_real_, mean_sig = NA_real_,
             tau_sig = NA_real_))
  }

  mu <- fit$par[1]
  tau <- fit$par[2]
  significant <- y / sei > crit
  y_sig <- y[significant]
  sei_sig <- sei[significant]
  if (!length(y_sig)) {
    return(c(mu = mu, tau = tau, mean_sig = NA_real_, tau_sig = NA_real_))
  }

  if (tau <= 1e-10) {
    posterior_mean_abs <- rep(abs(mu), length(y_sig))
    posterior_second <- rep(mu^2, length(y_sig))
  } else {
    marginal_sd <- sqrt(sei_sig^2 + tau^2)
    positive_weight <- dnorm(y_sig, mean = mu, sd = marginal_sd)
    negative_weight <- dnorm(-y_sig, mean = mu, sd = marginal_sd)
    positive_weight <- positive_weight /
      pmax(positive_weight + negative_weight, .Machine$double.xmin)

    posterior_variance <- 1 / (1 / tau^2 + 1 / sei_sig^2)
    posterior_sd <- sqrt(posterior_variance)
    posterior_mean_positive <- posterior_variance *
      (mu / tau^2 + y_sig / sei_sig^2)
    posterior_mean_negative <- posterior_variance *
      (mu / tau^2 - y_sig / sei_sig^2)

    folded_normal_mean <- function(mean, sd) {
      sd * sqrt(2 / pi) * exp(-mean^2 / (2 * sd^2)) +
        mean * (2 * pnorm(mean / sd) - 1)
    }
    mean_abs_positive <- folded_normal_mean(
      posterior_mean_positive, posterior_sd
    )
    mean_abs_negative <- folded_normal_mean(
      posterior_mean_negative, posterior_sd
    )
    posterior_mean_abs <-
      positive_weight * mean_abs_positive +
      (1 - positive_weight) * mean_abs_negative
    posterior_second <-
      positive_weight * (posterior_variance + posterior_mean_positive^2) +
      (1 - positive_weight) *
        (posterior_variance + posterior_mean_negative^2)
  }

  mean_sig <- mean(posterior_mean_abs)
  tau_sig <- sqrt(max(mean(posterior_second) - mean_sig^2, 0))
  c(mu = mu, tau = tau, mean_sig = mean_sig, tau_sig = tau_sig)
}


# Heterogeneity of true effects among SIGNIFICANT studies, estimated directly
# on the effect-size scale. Unlike the ncp mixture, this conditions on each
# study's own sei, so variation in precision (se) is not mistaken for effect
# heterogeneity -- the confound that floors the mixture's tau_sig at high mean.
#   theta ~ N(mu, tau^2);  y | theta ~ N(theta, sei^2)
#   significance selection |y| > crit*sei is built into the likelihood
#   (truncated normal), so the estimator is consistent rather than inflated.
# tau_sig = SD(|theta|) over the significant subset, via posterior moments.
fit_tau_sig <- function(yi, sei, crit = 1.96, folded = TRUE) {
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  y <- yi[keep]; s <- sei[keep]
  sig <- abs(y) / s > crit
  y <- y[sig]; s <- s[sig]; a <- abs(y); k <- length(a)
  if (k < 3L) return(c(mu = NA_real_, tau = NA_real_,
                       mean_sig = NA_real_, tau_sig = NA_real_))
  cs <- crit * s; tiny <- .Machine$double.xmin

  # P(|Y| > crit*s) for Y ~ N(mu, sqrt(tau^2 + s^2))
  p_sig <- function(mu, sig_i) pnorm((mu - cs) / sig_i) + pnorm((-cs - mu) / sig_i)

  nll <- function(par) {
    mu <- par[1]; tau <- exp(par[2]); sig_i <- sqrt(tau^2 + s^2)
    dens <- if (folded) dnorm(a, mu, sig_i) + dnorm(-a, mu, sig_i)
            else        dnorm(y, mu, sig_i)
    ll <- log(pmax(dens, tiny)) - log(pmax(p_sig(mu, sig_i), tiny))
    if (!all(is.finite(ll))) return(1e10)
    -sum(ll)
  }

  mu0 <- mean(a); tau0 <- max(stats::sd(a) / 2, 0.05)
  fit <- tryCatch(
    stats::optim(c(mu0, log(tau0)), nll, method = "L-BFGS-B",
                 lower = c(0, log(1e-4)), upper = c(5, log(5)),
                 control = list(maxit = 500)),
    error = function(e) NULL)
  if (is.null(fit) || !all(is.finite(fit$par)))
    return(c(mu = NA_real_, tau = NA_real_, mean_sig = NA_real_, tau_sig = NA_real_))

  mu <- fit$par[1]; tau <- exp(fit$par[2])

  # posterior moments of theta given each observed |y| (folded, normal-normal)
  sig_i <- sqrt(tau^2 + s^2); v <- 1 / (1 / tau^2 + 1 / s^2); sd_p <- sqrt(v)
  m_pos <- v * (mu / tau^2 + a / s^2)
  m_neg <- v * (mu / tau^2 - a / s^2)
  if (folded) {
    w_pos <- dnorm(a, mu, sig_i); w_neg <- dnorm(-a, mu, sig_i)
    w_pos <- w_pos / pmax(w_pos + w_neg, tiny)
  } else { w_pos <- rep(1, k); m_neg <- m_pos }
  fold_mean <- function(m, sd) sd * sqrt(2 / pi) * exp(-m^2 / (2 * sd^2)) +
    m * (2 * pnorm(m / sd) - 1)
  e_abs <- w_pos * fold_mean(m_pos, sd_p) + (1 - w_pos) * fold_mean(m_neg, sd_p)
  e_sq  <- w_pos * (v + m_pos^2)          + (1 - w_pos) * (v + m_neg^2)

  mean_sig <- mean(e_abs)
  tau_sig  <- sqrt(max(mean(e_sq) - mean_sig^2, 0))
  c(mu = mu, tau = tau, mean_sig = mean_sig, tau_sig = tau_sig)
}


#######################################################
### End of Get Densities
#######################################################




#####################################################################
### Compute Power Function (New: Discrete & Continuous
#####################################################################


get_estimates = function(
    est_inp, test_stat, yi, sei, int_beg, int_end,
    x_lim_min, x_lim_max) {


  #-------------------------------------------------------------------
  # helper: estimate a sampling-error function sei(z) over the z grid
  #-------------------------------------------------------------------

  estimate_sei_function <- function(
      test_stat, sei, grid, int_beg, int_end, min_k = 10) {

    dat <- data.frame(test_stat = abs(test_stat), sei = sei)
    dat <- dat[is.finite(dat$test_stat) & is.finite(dat$sei) & dat$sei > 0, ]

    dat_full <- dat[dat$test_stat <= int_end, ]

    if (nrow(dat_full) < min_k) {
      return(data.frame(grid = grid, sei_grid = rep(NA_real_, length(grid))))
    }

    # full-range SE function (fallback / non-selecting case)
    fit_full <- lm(sei ~ test_stat, data = dat_full)
    b1       <- coef(fit_full)[["test_stat"]]
    if (b1 > 0) b1 <- 0
    sei_out  <- coef(fit_full)[["(Intercept)"]] + b1 * grid

    # selecting case: fit on the significant interval and extrapolate
    if (int_beg > 0) {
      dat_int <- dat[dat$test_stat >= int_beg & dat$test_stat <= int_end, ]
    if (nrow(dat_int) >= min_k) {
      fit_int <- lm(sei ~ test_stat, data = dat_int)
      b1      <- coef(fit_int)[["test_stat"]]
      if (b1 > 0) b1 <- 0
      sei_out <- coef(fit_int)[["(Intercept)"]] + b1 * grid
      }
    }

    sei_out <- pmax(sei_out, min(dat$sei))   # floor at smallest observed SE
    data.frame(grid = grid, sei_grid = sei_out)
  }



  #-------------------------------------------------------------------
  # helper: empirical-Bayes effect-size estimates
  #-------------------------------------------------------------------

  get_es_estimates <- function(
      grid,
      pow_grid,
      sei_grid,
      w_grid,
      ncp_ext,
      w_all_ext,
      test_stat,
      yi,
      directional,
      folded,
      crit,
      int_beg,
      int_end,
      x_lim_min,
      x_lim_max,
      int_loc = 0.2
     ) {

    xz <- grid[2] - grid[1]

    ok       <- is.finite(grid) & is.finite(w_grid) & is.finite(sei_grid)
    stopifnot(length(grid) == length(w_grid), length(grid) == length(sei_grid))
    grid     <- grid[ok]
    w_grid   <- w_grid[ok]
    sei_grid <- sei_grid[ok]

    # normalize reconstructed density to probability weights
    w_grid <- w_grid / sum(w_grid, na.rm = TRUE)

    # ---- Empirical-Bayes posterior-mean ncp given observed z ---------------
    #   ncp_hat(z) = sum_k w_k phi(z - mu_k) mu_k / sum_k w_k phi(z - mu_k)
    if (!folded) {
      lik <- outer(grid, ncp_ext, function(x, mu) dnorm(x - mu))
    } else {
      lik <- outer(grid, ncp_ext, function(x, mu) dnorm(x - mu) + dnorm(x + mu))
    }

    # In a nondirectional/folded analysis the estimand is effect magnitude.
    # Continuous components can contain both positive and negative latent NCPs,
    # so compute E(|NCP| | |z|), not |E(NCP | |z|)|. The latter cancels to
    # nearly zero for a symmetric component centered at zero.
    posterior_ncp_value <- if (directional) ncp_ext else abs(ncp_ext)
    num <- as.numeric(lik %*% (w_all_ext * posterior_ncp_value))
    den <- as.numeric(lik %*% w_all_ext)

    ncp_second_grid <- as.numeric(lik %*% (w_all_ext * ncp_ext^2)) / den
    ncp_second_grid[!is.finite(ncp_second_grid)] <- NA_real_

    posterior_ncp <- sweep(lik, 2, w_all_ext, "*")
    posterior_ncp <- sweep(posterior_ncp, 1, den, "/")
    posterior_ncp[!is.finite(posterior_ncp)] <- 0

    ncp_grid <- num / den
    ncp_grid[!is.finite(ncp_grid)] <- NA_real_

    es_grid <- ncp_grid * sei_grid

    es_second_grid <- ncp_second_grid * sei_grid^2

    # ---- effect-size heterogeneity (weighted SD of the shrunken ES) --------
    m      <- sum(w_grid * es_grid)
    es_tau <- sqrt(sum(w_grid * (es_grid - m)^2))

    # ---- Extreme-value (tail) correction -----------------------------------
    if (directional) {
      es_ext_neg <- mean(yi[test_stat < -int_end])
    } else {
      es_ext_neg <- abs(mean(yi[test_stat < -int_end]))
    }
    w_ext_neg <- mean(test_stat < -int_end)

    es_ext_pos <- mean(yi[test_stat > int_end])
    w_ext_pos  <- mean(test_stat > int_end)

    grid_share <- 1 - w_ext_neg - w_ext_pos

    grid_ext    <- grid
    es_grid_ext <- es_grid
    es_second_grid_ext <- es_second_grid
    w_grid_ext  <- w_grid * grid_share
    if (directional) sig_grid <- grid > crit else sig_grid <- abs(grid) > crit
    sig_ext <- sig_grid

    if (w_ext_neg > 0) {
      grid_ext    <- c(min(grid) - xz, grid_ext)
      es_grid_ext <- c(es_ext_neg,      es_grid_ext)
      es_second_grid_ext <- c(es_ext_neg^2, es_second_grid_ext)
      w_grid_ext  <- c(w_ext_neg,       w_grid_ext)
      sig_ext     <- c(TRUE,            sig_ext)   # extreme tail is always significant
    }
    if (w_ext_pos > 0) {
      grid_ext    <- c(grid_ext,    max(grid) + xz)
      es_grid_ext <- c(es_grid_ext, es_ext_pos)
      es_second_grid_ext <- c(es_second_grid_ext, es_ext_pos^2)
      w_grid_ext  <- c(w_grid_ext,  w_ext_pos)
      sig_ext     <- c(sig_ext,     TRUE)
    }

    # ---- weighted MEAN over all reconstructed results ----------------------
    es_mean_all <- sum(es_grid_ext * w_grid_ext, na.rm = TRUE)

    # ---- weighted mean over reconstructed SIGNIFICANT results only ---------
    if (any(sig_ext) && sum(w_grid_ext[sig_ext], na.rm = TRUE) > 0) {
      w_sig       <- w_grid_ext[sig_ext] / sum(w_grid_ext[sig_ext], na.rm = TRUE)
      es_mean_sig <- sum(es_grid_ext[sig_ext] * w_sig, na.rm = TRUE)
      es_second_sig <- sum(es_second_grid_ext[sig_ext] * w_sig, na.rm = TRUE)
      es_tau_sig <- sqrt(max(es_second_sig - es_mean_sig^2, 0))
    } else {
      es_mean_sig <- NA_real_
      es_tau_sig <- NA_real_
    }


    # ---- posterior quantiles for SIGNIFICANT results -----------------------
    # Use the complete posterior distribution of latent effects. Quantiles of
    # posterior means alone would make heterogeneity look artificially small.
    theta_grid <- outer(sei_grid, ncp_ext, "*")
    if (!directional) theta_grid <- abs(theta_grid)
    joint_weight <- posterior_ncp * (w_grid * grid_share)

    sig_values <- as.numeric(theta_grid[sig_grid, , drop = FALSE])
    sig_weights <- as.numeric(joint_weight[sig_grid, , drop = FALSE])
    if (w_ext_neg > 0) {
      sig_values <- c(sig_values, es_ext_neg)
      sig_weights <- c(sig_weights, w_ext_neg)
    }
    if (w_ext_pos > 0) {
      sig_values <- c(sig_values, es_ext_pos)
      sig_weights <- c(sig_weights, w_ext_pos)
    }

    weighted_quantile <- function(values, weights, probs) {
      keep <- is.finite(values) & is.finite(weights) & weights > 0
      values <- values[keep]
      weights <- weights[keep]
      if (!length(values) || sum(weights) <= 0) return(rep(NA_real_, length(probs)))
      ord <- order(values)
      values <- values[ord]
      cum_weight <- cumsum(weights[ord]) / sum(weights)
      vapply(probs, function(prob) values[which(cum_weight >= prob)[1]], numeric(1))
    }

    es_sig_quantiles <- weighted_quantile(
      sig_values, sig_weights, probs = c(.25, .50, .75)
    )

    # ---- weighted MEDIAN over all reconstructed results --------------------
    o      <- order(es_grid_ext)
    es_srt <- es_grid_ext[o]
    w_srt  <- w_grid_ext[o]
    keep   <- is.finite(es_srt) & is.finite(w_srt) & w_srt > 0
    es_srt <- es_srt[keep]
    w_srt  <- w_srt[keep]

    if (length(es_srt) >= 2 && sum(w_srt) > 0) {
      cum_w     <- cumsum(w_srt) / sum(w_srt)
      es_median <- approx(cum_w, es_srt, xout = 0.5, rule = 2,
                          ties = list("ordered", mean))$y
    } else {
      es_median <- NA_real_
    }

    # ---- local effect sizes by z-bin ----------------------------------------
    z_breaks <- seq(x_lim_min, x_lim_max, by = int_loc)
    if (tail(z_breaks, 1) < int_end) z_breaks <- c(z_breaks, int_end)

    local_es <- matrix(NA_real_, nrow = 6, ncol = length(z_breaks) - 1)
    rownames(local_es) <- c("z_min", "z_max", "weight", "mean_z", "es", "tau")
    colnames(local_es) <- paste0(head(z_breaks, -1), "-", tail(z_breaks, -1))

    for (j in seq_len(length(z_breaks) - 1)) {
      in_bin <- grid >= z_breaks[j] & grid < z_breaks[j + 1]
      if (j == length(z_breaks) - 1) {
        in_bin <- grid >= z_breaks[j] & grid <= z_breaks[j + 1]
      }
      if (!any(in_bin)) next

      bin_weight <- sum(w_grid[in_bin], na.rm = TRUE)
      if (!is.finite(bin_weight) || bin_weight <= 0) next

      w_bin <- w_grid[in_bin] / bin_weight

      local_es["z_min",  j] <- z_breaks[j]
      local_es["z_max",  j] <- z_breaks[j + 1]
      local_es["mean_z", j] <- sum(grid[in_bin] * w_bin, na.rm = TRUE)
      local_es["weight", j] <- bin_weight
      local_mean <- sum(es_grid[in_bin] * w_bin, na.rm = TRUE)
      local_second <- sum(es_second_grid[in_bin] * w_bin, na.rm = TRUE)
      local_es["es",  j] <- local_mean
      local_es["tau", j] <- sqrt(max(local_second - local_mean^2, 0))
    }

    list(
      es_mean_all   = es_mean_all,
      es_median_all = es_median,
      es_mean_sig   = es_mean_sig,
      es_q25_sig    = es_sig_quantiles[1],
      es_median_sig = es_sig_quantiles[2],
      es_q75_sig    = es_sig_quantiles[3],
      es_tau_sig    = es_tau_sig,
      es_tau        = es_tau,
      local_es      = local_es[5, ],
      local_tau     = local_es[6, ],
      local_w       = local_es[3, ]
    )

  } ### EOF get_es_estimates


  ################################
  ### Power Computation
  ################################

  compute_power <- function(est_inp, int_beg, int_end, grid,
        test_stat, crit, directional, alpha,
        x_lim_min, x_lim_max, int_loc, folded,
        z_sd = NULL, df = NULL) {

    # All downstream calculations operate on latent-NCP points. A discrete
    # component maps to one point; a continuous component is integrated with
    # normal quadrature. The fitted w_inp values describe selected components,
    # so quadrature weights must first be conditioned on selection.
    latent_components <- expand_ncp_components(
      ncp = est_inp$ncp,
      weight = est_inp$w_inp,
      z_sd = est_inp$z_sd
    )
    ncp   <- latent_components$ncp
    z_sd  <- latent_components$z_sd

    sel_crit <- max(0, int_beg)
    if (directional) {
      latent_selection_probability <-
        pnorm(ncp, mean = sel_crit) + pnorm(-sel_crit, mean = ncp)
    } else {
      latent_selection_probability <-
        pnorm(abs(ncp), mean = sel_crit) +
        pnorm(-sel_crit, mean = abs(ncp))
    }
    component_selection_probability <- as.numeric(tapply(
      latent_components$quadrature_weight * latent_selection_probability,
      latent_components$component_id,
      sum
    ))
    component_selection_probability <- pmax(
      component_selection_probability, .Machine$double.xmin
    )

    w_inp <- latent_components$component_weight *
      latent_components$quadrature_weight * latent_selection_probability /
      component_selection_probability[latent_components$component_id]

    # Preserve component-level population weights for fitted-density output.
    # The latent version below is used for power and effect-size integration.
    w_all_component <- est_inp$w_inp / component_selection_probability
    w_all_component <- w_all_component / sum(w_all_component)

    xz      <- grid[2] - grid[1]
    k_tests <- length(test_stat)

    # extreme-value input PROPORTIONS (of total)
    ext_neg_inp <- sum(test_stat < -int_end) / k_tests
    ext_pos_inp <- sum(test_stat >  int_end) / k_tests

    # component power (two-tailed, with sign error)
    pow_neg <- pnorm(-crit, ncp)
    pow_pos <- 1 - pnorm(crit, ncp)
    pow_dir <- ifelse(ncp < 0, pow_neg, pow_pos)
    pow     <- pow_neg + pow_pos

    # power conditional on selection interval; sel_crit floors window bound at 0
    pow_sel <- latent_selection_probability

    # selection correction: input -> population weights (grid only)
    w_all <- w_inp / pow_sel
    w_all <- w_all / sum(w_all)

    # ---- extend vectors with extreme components (used only for EDR/ERR) ----
    add_neg <- (ext_neg_inp > 0 & directional == FALSE)
    add_pos <- (ext_pos_inp > 0)

    # grid weight after removing whatever extreme mass is present (computed once)
    grid_share <- 1 - (if (add_neg) ext_neg_inp else 0) -
                      (if (add_pos) ext_pos_inp else 0)

    ncp_ext     <- ncp
    pow_ext     <- pow
    pow_dir_ext <- pow_dir
    pow_sel_ext <- pow_sel
    w_inp_ext   <- w_inp * grid_share

    if (add_neg) {
      ncp_ext     <- c(-(int_end + 1), ncp_ext)
      pow_ext     <- c(1, pow_ext)
      pow_dir_ext <- c(1, pow_dir_ext)
      pow_sel_ext <- c(1, pow_sel_ext)
      w_inp_ext   <- c(ext_neg_inp, w_inp_ext)
    }

    if (add_pos) {
      ncp_ext     <- c(ncp_ext, int_end + 1)
      pow_ext     <- c(pow_ext, 1)
      pow_dir_ext <- c(pow_dir_ext, 1)
      pow_sel_ext <- c(pow_sel_ext, 1)
      w_inp_ext   <- c(w_inp_ext, ext_pos_inp)
    }

    stopifnot(length(pow_ext)     == length(w_inp_ext),
              length(pow_dir_ext) == length(w_inp_ext),
              length(pow_sel_ext) == length(w_inp_ext))

    w_all_ext <- w_inp_ext / pow_sel_ext
    w_all_ext <- w_all_ext / sum(w_all_ext)

    # EDR: average power over all studies (pre-selection)
    EDR <- sum(w_all_ext * pow_ext)

    # ERR: average directional power over significant studies
    w_sig_ext <- w_all_ext * pow_ext
    w_sig_ext <- w_sig_ext / sum(w_sig_ext)
    ERR <- sum(w_sig_ext * pow_dir_ext)

    ERR <- min(max(ERR, alpha / 2), 1)
    EDR <- min(max(EDR, alpha),     1)

    # ---- mixture density at each grid point (grid-only ncp and w_all) ----
    # test_stat is folded/absolute, so the density of |Z| (or |T|) for a
    # component at mu needs the mirror term when folded = TRUE. build_dens
    # handles the fold, the component SDs (z_sd), and the z-vs-t choice (df).
    dens <- build_dens(d_x = grid, ncp = ncp, z_sd = z_sd, df = df,
                          curve_type = if (is.null(df)) "z" else "t",
                          folded = folded)

    lik_mat <- t(dens)                       # rows = grid, cols = components
    w_grid  <- as.vector(lik_mat %*% w_all)

    # local power: posterior-weighted average of component power
    # (numerator and denominator share the same w_all-based lik_mat, so this
    #  ratio is scale-invariant -- no normalization of w_grid needed or wanted)
    pow_grid <- as.vector((lik_mat %*% (w_all * pow)) / w_grid)

    # bin into intervals, density-weighted average within each bin
    int  <- seq(x_lim_min, x_lim_max, by = int_loc)
    bins <- cut(grid, breaks = int, include.lowest = TRUE)

    local_power <- as.vector(tapply(seq_along(grid), bins, function(idx) {
      denom <- sum(w_grid[idx])
      if (denom <= 0 || !is.finite(denom)) return(NA_real_)
      sum(pow_grid[idx] * w_grid[idx]) / denom
    }))

    midpoints          <- (int[-length(int)] + int[-1]) / 2
    names(local_power) <- paste0("lp.", midpoints)


    # ---- heterogeneity of the non-centrality distribution (z-scale) --------
    # Weighted SD of the fitted noncentrality mixture, including the existing
    # extreme component for observations beyond int_end. Previously tau used
    # only the interior grid even though EDR and ERR already used this tail.
    # This matches the folded simulation target sd(abs(ncp)) more closely.
    # A previous method-of-moments calculation subtracted 1 from the variance
    # of the folded/truncated z density. That variance need not exceed 1, so the
    # calculation frequently collapsed to zero and forced ncp_tau to zero.
    ncp_tau_values <- if (folded) abs(ncp_ext) else ncp_ext
    ncp_mean  <- sum(w_all_ext * ncp_tau_values)
    ncp_var   <- sum(w_all_ext * ncp_tau_values^2) - ncp_mean^2
    ncp_tau  <- sqrt(max(ncp_var, 0))

    # shape diagnostic: predicted vs observed nonsignificant z
    upper <- qnorm(1 - alpha)
    lower <- 0
    idx   <- grid > lower & grid < upper

    d_med  <- NA_real_
    d_mean <- NA_real_
    if (sum(idx) > 1) {
      x_ns <- grid[idx]
      pred <- w_grid[idx]
      pred <- pred / (sum(pred) * xz)

      pred_cdf  <- cumsum(pred) * xz
      pred_cdf  <- pred_cdf / pred_cdf[length(pred_cdf)]
      pred_med  <- approx(pred_cdf, x_ns, xout = 0.5, rule = 2)$y
      pred_mean <- sum(x_ns * pred) * xz

      obs <- test_stat[test_stat > lower & test_stat < upper]
      if (length(obs) > 0) {
        d_med  <- median(obs) - pred_med
        d_mean <- mean(obs)   - pred_mean
      }
    }

    list(
      EDR            = EDR,
      ERR            = ERR,
      w_all          = w_all_component,
      w_all_latent   = w_all,
      w_all_ext      = w_all_ext,
      pow_ext        = pow_ext,
      pow_dir_ext    = pow_dir_ext,
      pow_sel_ext    = pow_sel_ext,
      ncp_ext        = ncp_ext,
      latent_components = latent_components,
      grid           = grid,
      w_grid         = w_grid,
      pow_grid       = pow_grid,
      local_power    = local_power,
      ncp_tau        = ncp_tau,
      shape_d_median = d_med,
      shape_d_mean   = d_mean
    )

  } ### EOF compute_power


  ################################################################
  ### Main Power # ppp2
  ################################################################

  xz     <- 0.01
  grid <- seq(x_lim_min, x_lim_max, by = xz)

  est_out <- compute_power(
    est_inp     = est_inp,
    int_beg     = int_beg,
    int_end     = int_end,
    grid        = grid,
    test_stat   = test_stat,
    crit        = crit,
    directional = directional,
    folded      = folded,
    alpha       = alpha,
    x_lim_min   = x_lim_min,
    x_lim_max   = x_lim_max,
    int_loc     = int_loc
   )

  est_out$local_es      <- NULL
  est_out$local_tau     <- NULL
  est_out$es_mean_all   <- NULL
  est_out$es_median_all <- NULL
  est_out$es_mean_sig   <- NULL
  est_out$es_q25_sig    <- NULL
  est_out$es_median_sig <- NULL
  est_out$es_q75_sig    <- NULL
  est_out$es_tau_sig    <- NULL
  est_out$es_tau        <- NULL

  if (!is.null(yi)) {

    sei_est <- estimate_sei_function(
      test_stat = test_stat,
      sei       = sei,
      grid      = grid,
      int_beg   = int_beg,
      int_end   = int_end
    )

   es_est <- get_es_estimates(
      sei_grid    = sei_est[,2],
      grid        = est_out$grid,
      w_grid      = est_out$w_grid,
      ncp_ext     = est_out$ncp_ext,
      w_all_ext   = est_out$w_all_ext,
      test_stat   = test_stat,  
      yi          = yi,         
      directional = directional,
      folded      = folded,
      crit        = crit,
      int_beg     = int_beg,
      int_end     = int_end,
      int_loc     = int_loc,
      x_lim_min   = x_lim_min,
      x_lim_max   = x_lim_max
    )
 
    est_out$local_es      <- es_est$local_es
    est_out$local_tau     <- es_est$local_tau
    est_out$es_mean_all   <- es_est$es_mean_all
    est_out$es_median_all <- es_est$es_median_all   # FIX: now propagated
    est_out$es_mean_sig   <- es_est$es_mean_sig
    est_out$es_q25_sig    <- es_est$es_q25_sig
    est_out$es_median_sig <- es_est$es_median_sig
    est_out$es_q75_sig    <- es_est$es_q75_sig
    est_out$es_tau_sig    <- es_est$es_tau_sig
    est_out$es_tau        <- es_est$es_tau

  } # EOF es 

  est_out
 
  return(est_out)
 
} ### EOF Compute Estimates
 



############################################################
### SLOPE DIAGNOSTICS
############################################################

get_slope <- function(int, from, to, bw, spline_width = .5) {

  d <- get_densities(int, bw = bw, from = from, to = to,
       width = spline_width, augment = TRUE)
  fit <- lm(d$zy ~ d$zx)
  slope <- unname(coef(fit)[2])
  se <- unname(summary(fit)$coefficients[2, 2])
  ci <- confint(fit)[2, ]
  
  out <- c(slope = slope, se = se, ci_lb = ci[1], ci_ub = ci[2])
  out
} # EOF get_slope


######################################################
### END OF SLOPE DIAGNOSTICS
######################################################



##################################################
### ROOT OF FUNCTION
##################################################


#ooo

run_OF <- function(test_stat,yi=NULL,sei=NULL,ODR,int_beg,int_end,ncp,z_sd, 
                    cola = "springgreen4",bw_est = .2,
                    dens_pre = NULL, d_x_pre = NULL, width = 1, 
                    ncp_fixed, z_sd_fixed, w_fixed,
                    x_lim_min,x_lim_max,augment,density_step,
                    max_iter,tol_criterion) {


  testing = FALSE
  #testing = TRUE
  #width = 1

  int = test_stat[test_stat >= int_beg & test_stat <= int_end]
  #hist(int)

  ODR = mean(test_stat > 1.96)

  k_ncp <- length(ncp)

  #from = int_beg; to = int_end; width = 1
  ## ---- Observed density ----
  dens <- get_densities(int, bw = bw_est, from = int_beg, 
                         to = int_end, width = width, augment = augment,
                         density_step = density_step)

  d_x       <- dens[, 1]
  o_d_y     <- dens[, 2]
  n_bars    <- length(d_x) 
  bar_width <- d_x[2] - d_x[1]

  use_pre_dens <- k_ncp > 1 &&
                     ncp_fixed == TRUE &&
                     z_sd_fixed == TRUE


  ## ---- theta = [weights, ncp, z_sd] — always ----

  startval <- c(rep(1/k_ncp, k_ncp), ncp, z_sd)

  lowlim   <- c(if (k_ncp == 1) 1 else rep(0, k_ncp),
                if (ncp_fixed) ncp else rep(0,k_ncp),
                if (z_sd_fixed) z_sd else rep(1,k_ncp)
			)
  highlim  <- c(rep(1, k_ncp),
                if (ncp_fixed) ncp else rep(max(ncp),k_ncp),
                if (z_sd_fixed) z_sd else rep(10,k_ncp)
			)

  ## ---- Positions (always the same) ----
  idx_wt  <- 1:k_ncp
  idx_ncp <- (k_ncp + 1):(2 * k_ncp)
  idx_z_sd <- (2 * k_ncp + 1):(3 * k_ncp)


  curve_fitting <- function(theta) {

    wt <- theta[idx_wt]
    wt <- wt / sum(wt)
 
    if (use_pre_dens && ncp_fixed && z_sd_fixed) {
      dens_now <- dens_pre
    } else {
      dens_now <- build_dens(d_x, theta[idx_ncp], theta[idx_z_sd], df, 
           curve_type, folded=folded)
    }
    dens_now <- dens_now / (rowSums(dens_now) * bar_width)

    e_d_y <- colSums(dens_now * wt)

    rmse = sqrt(mean((e_d_y - o_d_y)^2))

    scale = 1
    ## ---- Diagnostic plot ----
    if (testing && runif(1) > 0.6) {
      plot(d_x, o_d_y * scale, type = 'l',
           xlim = c(x_lim_min, x_lim_max), ylim = c(0, ymax),
           xlab = "", ylab = "", axes = FALSE)
      lines(d_x, e_d_y * scale, col = "red3")
    }


    rmse

    } #EOF curve_fitting

  ## ---- Optimize ----
  auto <- nlminb(startval, curve_fitting,
                 lower = lowlim, upper = highlim,
                 control = list(eval.max = max_iter, abs.tol = tol_criterion))

  of_fit =  list(
    w_inp   = auto$par[idx_wt]/sum(auto$par[idx_wt]),
    ncp     = auto$par[idx_ncp],
    z_sd    = auto$par[idx_z_sd],
    fit     = auto$objective
  )

  of_fit$w_inp

  est_inp = of_fit
 
  ## ---- Compute power for each bootstrap result ----
  res_est = get_estimates(
        est_inp = est_inp, 
        test_stat = test_stat,
        int_beg = int_beg,
	        int_end = int_end, 
        yi = yi, 
        sei = sei,
        x_lim_min = x_lim_min,
        x_lim_max = x_lim_max)

  ret = list(
    w_inp           = est_inp$w_inp,
    ncp             = est_inp$ncp,
    z_sd            = est_inp$z_sd,
    model_fit       = est_inp$fit,
    ODR             = ODR,   
    EDR             = res_est$EDR,
    ERR             = res_est$ERR,
    w_all           = res_est$w_all, 
    local_power     = res_est$local_power,
    ncp_tau         = res_est$ncp_tau,
    shape_d_median  = res_est$shape_d_median,
    shape_d_mean    = res_est$shape_d_mean,
    local_es        = res_est$local_es,
    local_tau       = res_est$local_tau,
    es_mean_all     = res_est$es_mean_all,
    es_median_all   = res_est$es_median_all,
    es_mean_sig     = res_est$es_mean_sig,
    es_q25_sig      = res_est$es_q25_sig,
    es_median_sig   = res_est$es_median_sig,
    es_q75_sig      = res_est$es_q75_sig,
    es_tau_sig      = res_est$es_tau_sig,
    es_tau          = res_est$es_tau
  )

  ret
  
  return(ret)


} ### EOF run_OF



##################################################
### ROOT EM FUNCTION
##################################################


#eee

run_EM <- function(
  test_stat,
  yi, 
  sei, 
  int_beg,
  int_end,
  ncp,
  z_sd,
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  x_lim_min,
  x_lim_max,
  max_iter,
  tol_criterion) {


n_starts = 1

components = length(ncp)
w_inp = rep(1/components,components)

int = test_stat[test_stat >= int_beg & test_stat <= int_end]
ODR = mean(test_stat > 1.96);ODR

EM_ext_signed <- function(
  int,
  ncp,
  z_sd,
  int_beg = -Inf,        # signed lower window bound (e.g. -8, or -Inf for none)
  int_end =  Inf,        # signed upper window bound
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  max_iter,
  tol_criterion
) {

  x <- int
  x <- x[!is.na(x)]
  k_int <- length(x)
  components <- length(ncp)

  loglik_old <- -Inf
  loglik_trace <- numeric(max_iter)

  # ---- Newton update for a single component's mean (signed, truncated) ----
  update_ncp_newton_signed <- function(
    mu_init, sigma, tau_k, x, int_beg, int_end,
    max_iter = 20, tol_criterion, lower = -10, upper = 10
  ) {
    mu  <- mu_init
    n_k <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: sufficient statistic is just tau-weighted sum of x
      s1 <- sum(tau_k * x)

      # single truncation interval [int_beg, int_end] on the signed axis
      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      cmu <- pnorm(b) - pnorm(a)
      if (cmu <= 0 || !is.finite(cmu)) return(mu)

      phi_a <- dnorm(a); phi_b <- dnorm(b)
      delta <- (phi_a - phi_b) / cmu
      kappa <- (a * phi_a - b * phi_b) / cmu

      u <- (s1 - n_k * mu) / sigma^2 - n_k * delta / sigma
      h <- -n_k / sigma^2 - n_k / sigma^2 * (delta^2 + kappa)
      if (h >= 0) h <- -n_k / sigma^2

      step   <- u / h
      mu_new <- max(lower, min(upper, mu - step))
      if (abs(mu_new - mu) < tol_criterion) break
      mu <- mu_new
    }
    return(mu)
  }

  # ---- Newton update for a single component's sd (signed, truncated) ----
  update_z_sd_newton_signed <- function(
    sigma_init, mu, tau_k, x, int_beg, int_end,
    max_iter = 20, tol_criteiron, lower = 0.5, upper = 10
  ) {
    sigma <- sigma_init
    n_k   <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: second moment about mu directly
      s2 <- sum(tau_k * (x - mu)^2)

      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      cmu <- pnorm(b) - pnorm(a)
      if (cmu <= 0 || !is.finite(cmu)) return(sigma)

      phi_a <- dnorm(a); phi_b <- dnorm(b)
      kappa <- (a * phi_a - b * phi_b) / cmu

      u <- s2 / sigma^3 - n_k / sigma - n_k * kappa / sigma
      h <- -2 * n_k / sigma^2

      step      <- u / h
      sigma_new <- max(lower, min(upper, sigma - step))
      if (abs(sigma_new - sigma) < tol_criterion) break
      sigma <- sigma_new
    }
    return(sigma)
  }

  ######## EM start #####
  for (iter in 1:max_iter) {

    # ---- E-STEP (signed, no fold) ----
    log_g <- matrix(NA_real_, k_int, components)
    for (k in 1:components) {
      mu    <- ncp[k]
      sigma <- z_sd[k]

      # single signed density, no + dnorm(x, -mu)
      log_dens <- dnorm(x, mu, sigma, log = TRUE)

      # single truncation interval on the signed axis
      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      norm_const <- pnorm(b) - pnorm(a)
      if (norm_const <= 0 || !is.finite(norm_const)) return(NULL)

      log_g[, k] <- log_dens - log(norm_const)
    }

    log_num <- sweep(log_g, 2, log(w_inp), "+")
    log_denom <- apply(log_num, 1, function(v) {
      m <- max(v); m + log(sum(exp(v - m)))
    })
    tau <- exp(log_num - log_denom)

    # ---- M-STEP ----
    n_k <- colSums(tau)
    if (!w_fixed) {
      w_inp <- n_k / k_int
      w_inp <- pmax(w_inp, 1e-12)
      w_inp <- w_inp / sum(w_inp)
    }
    if (!ncp_fixed) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        ncp[k] <- update_ncp_newton_signed(
          mu_init = ncp[k], sigma = z_sd[k], tau_k = tau[, k],
          x = x, int_beg = int_beg, int_end = int_end
        )
      }
    }
    if (!z_sd_fixed) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        z_sd[k] <- max(1, update_z_sd_newton_signed(
          sigma_init = z_sd[k], mu = ncp[k], tau_k = tau[, k],
          x = x, int_beg = int_beg, int_end = int_end
        ))
      }
    }

    # ---- LOG-LIKELIHOOD ----
    loglik <- sum(log_denom)
    loglik_trace[iter] <- loglik
    if (!is.finite(loglik)) return(NULL)
    if (abs(loglik - loglik_old) < tol_criterion) break
    loglik_old <- loglik
  }

  list(w_inp = w_inp, ncp = ncp, z_sd = z_sd, fit = loglik)
}




EM_ext_folded <- function(
  int,
  ncp,
  z_sd,
  int_beg = 1.96,
  int_end = max(int),
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  max_iter,
  tol_criterion
) {


x <- int
x <- x[!is.na(x)]
k_int <- length(x)
components <- length(ncp)

loglik_old <- -Inf
loglik_trace <- numeric(max_iter)



update_ncp_newton_folded <- function(
  mu_init,
  sigma,
  tau_k,
  x,
  int_beg,
  int_end,
  max_iter,
  tol,criterion,
  lower    = -10,
  upper    = 10
) {
  mu  <- mu_init
  n_k <- sum(tau_k)

  for (iter in 1:max_iter) {
    # --- posterior fold probabilities (depend on current mu) ---
    f_pos <- dnorm(x,  mu, sigma)
    f_neg <- dnorm(x, -mu, sigma)
    p_pos <- f_pos / (f_pos + f_neg)

    # signed sufficient statistic
    r         <- 2 * p_pos - 1
    s1_signed <- sum(tau_k * x * r)

    # --- truncation correction ---
    a1 <- ( int_beg - mu) / sigma
    b1 <- ( int_end - mu) / sigma
    a2 <- (-int_end - mu) / sigma
    b2 <- (-int_beg - mu) / sigma

    cmu <- (pnorm(b1) - pnorm(a1)) +
           (pnorm(b2) - pnorm(a2))
    if (cmu <= 0 || !is.finite(cmu)) return(mu)

    phi_a1 <- dnorm(a1); phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2); phi_b2 <- dnorm(b2)

    delta <- (phi_a1 - phi_b1 + phi_a2 - phi_b2) / cmu
    kappa <- (a1*phi_a1 - b1*phi_b1 +
              a2*phi_a2 - b2*phi_b2) / cmu

    # --- score ---
    u <- (s1_signed - n_k * mu) / sigma^2 -
         n_k * delta / sigma

    # --- Hessian (with folding term) ---
    p_neg  <- 1 - p_pos
    h_fold <- 4 / sigma^4 * sum(tau_k * x^2 * p_pos * p_neg)
    h      <- h_fold - n_k / sigma^2 -
              n_k / sigma^2 * (delta^2 + kappa)

    if (h >= 0) h <- -n_k / sigma^2   # safeguard

    step   <- u / h
    mu_new <- max(lower, min(upper, mu - step))

    if (abs(mu_new - mu) < tol_criterion) break
    mu <- mu_new
  }
  return(mu)
}


update_z_sd_newton_folded <- function(
  sigma_init,
  mu,
  tau_k,
  x,
  int_beg,
  int_end,
  max_iter,
  tol_criterion,
  lower    = 0.5,
  upper    = 10
) {
  sigma <- sigma_init
  n_k   <- sum(tau_k)

  for (iter in 1:max_iter) {
    # --- posterior fold probabilities (depend on current sigma) ---
    f_pos <- dnorm(x,  mu, sigma)
    f_neg <- dnorm(x, -mu, sigma)
    p_pos <- f_pos / (f_pos + f_neg)
    p_neg <- 1 - p_pos

    # signed second moment
    s2_signed <- sum(tau_k * (p_pos * (x - mu)^2 +
                              p_neg * (x + mu)^2))

    # --- truncation correction ---
    a1 <- ( int_beg - mu) / sigma
    b1 <- ( int_end - mu) / sigma
    a2 <- (-int_end - mu) / sigma
    b2 <- (-int_beg - mu) / sigma

    cmu <- (pnorm(b1) - pnorm(a1)) +
           (pnorm(b2) - pnorm(a2))
    if (cmu <= 0 || !is.finite(cmu)) return(sigma)

    phi_a1 <- dnorm(a1); phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2); phi_b2 <- dnorm(b2)

    kappa <- (a1*phi_a1 - b1*phi_b1 +
              a2*phi_a2 - b2*phi_b2) / cmu

    # --- score ---
    u <- s2_signed / sigma^3 - n_k / sigma -
         n_k * kappa / sigma

    # --- Hessian (conservative approximation) ---
    h <- -2 * n_k / sigma^2

    step      <- u / h
    sigma_new <- max(lower, min(upper, sigma - step))

    if (abs(sigma_new - sigma) < tol_criterion) break
    sigma <- sigma_new
  }
  return(sigma)
}

    ######## EM start ##### 

    for (iter in 1:max_iter) {

    # ============================================================
    # E-STEP
    # ============================================================
    log_g <- matrix(NA_real_, k_int, components)
    for (k in 1:components) {
      mu    <- ncp[k]
      sigma <- z_sd[k]
      l1 <- dnorm(x,  mu, sigma, log = TRUE)
      l2 <- dnorm(x, -mu, sigma, log = TRUE)
      m  <- pmax(l1, l2)
      log_fold <- m + log(exp(l1 - m) + exp(l2 - m))
      a1 <- ( int_beg - mu) / sigma
      b1 <- ( int_end - mu) / sigma
      a2 <- (-int_end - mu) / sigma
      b2 <- (-int_beg - mu) / sigma
      norm_const <-
        (pnorm(b1) - pnorm(a1)) +
        (pnorm(b2) - pnorm(a2))
      if (norm_const <= 0 || !is.finite(norm_const))
        return(NULL)
      log_g[, k] <- log_fold - log(norm_const)
    }
    log_num <- sweep(log_g, 2, log(w_inp), "+")
    log_denom <- apply(log_num, 1, function(v) {
      m <- max(v)
      m + log(sum(exp(v - m)))
    })
    tau <- exp(log_num - log_denom)
    # ============================================================
    # M-STEP
    # ============================================================
    n_k <- colSums(tau)
    if (!w_fixed) {
      w_inp <- n_k / k_int
      w_inp <- pmax(w_inp, 1e-12)
      w_inp <- w_inp / sum(w_inp)
    }
    if (!ncp_fixed) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        ncp[k] <- update_ncp_newton_folded(
          mu_init = ncp[k],
          sigma   = z_sd[k],
          tau_k   = tau[, k],
          x       = x,
          int_beg = int_beg,
          int_end = int_end
        )
      }
    }
    if (!z_sd_fixed) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        z_sd[k] <- max(1, update_z_sd_newton_folded(
          sigma_init = z_sd[k],
          mu         = ncp[k],
          tau_k      = tau[, k],
          x          = x,
          int_beg    = int_beg,
          int_end    = int_end
        ))
      }
    }

    # ============================================================
    # LOG-LIKELIHOOD
    # ============================================================

    loglik <- sum(log_denom)
    loglik_trace[iter] <- loglik

    if (!is.finite(loglik))
      return(NULL)

    if (abs(loglik - loglik_old) < tol_criterion)
      break

    loglik_old <- loglik
  }

  list(
    w_inp        = w_inp,
    ncp          = ncp,
    z_sd         = z_sd,
    fit          = loglik
  )

  
}



### Main
### eee2

  components <- length(ncp)
  ncp_init   <- ncp
  z_sd_init  <- z_sd

  best_fit <- NULL
  best_ll  <- -Inf

  for (s in 1:n_starts) {

      raw <- rgamma(components, shape = 1)   # random Dirichlet
      w_s <- raw / sum(raw)

      if(folded) {

        fit <- EM_ext_folded(
          int        = int,
          ncp        = ncp,
          z_sd       = z_sd,
          int_beg    = int_beg,
          int_end    = int_end,
          ncp_fixed  = ncp_fixed,
          z_sd_fixed = z_sd_fixed,
          w_fixed    = w_fixed,
          max_iter   = max_iter,
          tol_criterion = tol_criterion
        ) 

      } else {

        fit <- EM_ext_signed(
          int        = int,
          ncp        = ncp,
          z_sd       = z_sd,
          int_beg    = int_beg,
          int_end    = int_end,
          ncp_fixed  = ncp_fixed,
          z_sd_fixed = z_sd_fixed,
          w_fixed    = w_fixed,
          max_iter   = max_iter,
          tol_criterion = tol_criterion
         ) 

      }

    if (!is.null(fit) && is.finite(fit$fit) && fit$fit > best_ll) {
      best_ll  <- fit$fit
      best_fit <- fit
      best_fit$start <- s   # track which start won
    }
  }

  em_est = best_fit

  #print(em.est)

  cp_inp = em_est

  ## ---- Compute power
  cp_res = get_estimates(
        est_inp = cp_inp, 
        test_stat = test_stat,
        int_beg = int_beg,
        int_end = int_end, 
        yi = yi, 
        sei = sei,
        x_lim_min = x_lim_min,
        x_lim_max = x_lim_max)

  ret = list(
    w_inp           = em_est$w_inp,
    ncp             = em_est$ncp,
    z_sd            = em_est$z_sd,
    model_fit       = em_est$fit,
    ODR             = ODR,   
    EDR             = cp_res$EDR,
    ERR             = cp_res$ERR,
    w_all           = cp_res$w_all, 
    local_power     = cp_res$local_power,
    ncp_tau         = cp_res$ncp_tau,
    shape_d_median  = cp_res$shape_d_median,
    shape_d_mean    = cp_res$shape_d_mean,
    local_es        = cp_res$local_es,
    local_tau       = cp_res$local_tau,
    es_mean_all     = cp_res$es_mean_all,
    es_median_all   = cp_res$es_median_all,
    es_mean_sig     = cp_res$es_mean_sig,
    es_q25_sig      = cp_res$es_q25_sig,
    es_median_sig   = cp_res$es_median_sig,
    es_q75_sig      = cp_res$es_q75_sig,
    es_tau_sig      = cp_res$es_tau_sig,
    es_tau          = cp_res$es_tau
  )

  ret

  return(ret)

}




##################################################
### Bootstrap
##################################################

### BOOTSTRAP FUNCTION FOR CI

#xxx

run_bootstrap <- function(
                   test_stat,
                   yi,
                   sei,
                   crit,   
                   ncp,
                   z_sd,
                   int_beg,
                   int_end,
                   directional,
                   folded,
                   dens_pre,
                   d_x_pre,
                   width = 1,
                   boot_iter,
                   ncp_fixed,
                   z_sd_fixed,
                   w_fixed,
                   cluster_id,
                   est_method,
                   x_lim_min,
                   x_lim_max,
                   bw_est,
                   augment,
                   density_step,
                   max_iter_boot,
                   tol_criterion,
                   parallel_boot,
                   cores,
                   seed
                  ) {

    

      k_ncp = length(ncp)

      if (!is.null(cluster_id)) {
   
          ## ---- Cluster Bootstrap ----

          stopifnot(length(test_stat) == length(cluster_id))	
          z_cluster_rows  <- split(seq_along(test_stat), cluster_id)
          z_cluster_names <- names(z_cluster_rows)
          n_z_clusters    <- length(z_cluster_names)

      } else {

          z_cluster_rows = NULL; z_cluster_names = NULL; n_z_clusters = NULL

      }

 
      var_list = c("est_method","test_stat","crit",
           "cluster_id","z_cluster_rows","z_cluster_names","n_z_clusters",
           "ncp", "z_sd","int_beg", "int_end", 
           "curve_type","directional","folded","x_lim_min", "x_lim_max", 
           "ncp_fixed", "z_sd_fixed","w_fixed",
            "get_estimates","expand_ncp_components",
            "fit_fixed_weight_selection","fit_folded_effect_model","fit_tau_sig",
            "use_fit_tau_sig",
            "alpha","int_loc",
           "density_step","tol_criterion","max_iter_boot",
           "use_folded_es_model"
           )


      if (!is.null(yi)) var_list = c(var_list,"yi","sei")

      bootstrap_lapply <- function(FUN, var_list, packages = character()) {
        detected_cores <- parallel::detectCores(logical = FALSE)
        if (!is.finite(detected_cores) || detected_cores < 1) detected_cores <- 1L
        requested_cores <- if (is.null(cores)) {
          max(1L, min(8L, floor(detected_cores * .8)))
        } else {
          as.integer(cores)
        }
        worker_count <- min(as.integer(boot_iter), requested_cores)

        safe_FUN <- function(i) {
          tryCatch(
            FUN(i),
            error = function(error) structure(
              list(message = conditionMessage(error)),
              class = "zcurve_bootstrap_error"
            )
          )
        }

        if (isTRUE(parallel_boot) && worker_count > 1L) {
          cl <- parallel::makeCluster(worker_count)
          on.exit(parallel::stopCluster(cl), add = TRUE)
          parallel::clusterExport(cl, varlist = unique(var_list), envir = environment())
          for (package in packages) {
            parallel::clusterCall(
              cl,
              function(package) suppressPackageStartupMessages(
                library(package, character.only = TRUE)
              ),
              package
            )
          }
          rng_seed <- if (is.null(seed)) {
            sample.int(.Machine$integer.max, 1L)
          } else {
            as.integer(seed)
          }
          parallel::clusterSetRNGStream(cl, iseed = rng_seed)
          result <- parallel::parLapply(cl, seq_len(boot_iter), safe_FUN)
        } else {
          if (!is.null(seed)) set.seed(seed)
          result <- lapply(seq_len(boot_iter), safe_FUN)
        }

        failed <- vapply(
          result,
          function(x) is.null(x) || inherits(x, "zcurve_bootstrap_error"),
          logical(1)
        )
        if (any(failed)) {
          warning(
            sum(failed), " of ", length(result),
            " bootstrap replicates failed and were omitted.",
            call. = FALSE
          )
        }
        successful <- result[!failed]
        attr(successful, "failed_replicates") <- sum(failed)
        attr(successful, "worker_count") <- if (isTRUE(parallel_boot)) worker_count else 1L
        successful
      }

      boot_res = NULL

      if (est_method == "OF") {

         var_list = c(var_list,"run_OF","build_dens", "get_densities",
                  "bw_est","width","dens_pre", "d_x_pre","augment")

         boot_res <- bootstrap_lapply(function(i) {

         if (!is.null(cluster_id)) {
              sampled  <- sample(z_cluster_names, n_z_clusters, replace = TRUE)
              rows     <- unlist(z_cluster_rows[sampled], use.names = FALSE)
              boot_all <- test_stat[rows]
          } else {
              rows     <- sample(seq_along(test_stat), length(test_stat), replace = TRUE)
              boot_all <- test_stat[rows]
          }

          boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit_res <- #tryCatch(
                 run_OF(
                   test_stat     = boot_all,
                   yi            = yi_boot,
                   sei           = sei_boot,
                   int_beg       = int_beg,
                   int_end       = int_end,
                   ncp           = ncp,
                   z_sd          = z_sd,
                   dens_pre      = dens_pre,
                   d_x_pre       = d_x_pre,
                   bw_est        = bw_est, 
                   width         = width,
                   ncp_fixed     = ncp_fixed,
                   z_sd_fixed    = z_sd_fixed,
                   w_fixed       = w_fixed,
                   x_lim_min     = x_lim_min,
                   x_lim_max     = x_lim_max,
                   augment       = augment,
                   density_step  = density_step,
                   max_iter      = max_iter_boot,
                   tol_criterion = tol_criterion
                 )
	        #  ,error = function(e) NULL )
          if (!is.null(yi_boot)) {
            selection_fit <- fit_fixed_weight_selection(
              yi = yi_boot,
              sei = sei_boot,
              ODR = fit_res$ODR,
              EDR = fit_res$EDR,
              crit = crit,
              selection = "one_sided"
            )
            fit_res$es_tau <- unname(selection_fit["tau"])
            if (use_folded_es_model && any(fit_res$z_sd > 1 + 1e-10, na.rm = TRUE)) {
              effect_fit <- fit_folded_effect_model(
                yi_boot, sei_boot, fit_res$ODR, fit_res$EDR, crit
              )
              if (is.finite(effect_fit["mean_sig"])) {
                fit_res$es_mean_sig <- unname(effect_fit["mean_sig"])
              }
              if (is.finite(effect_fit["tau_sig"])) {
                fit_res$es_tau_sig <- unname(effect_fit["tau_sig"])
              }
            }

            if (use_fit_tau_sig) {
              tau_sig_fit <- fit_tau_sig(yi = yi_boot, sei = sei_boot,
                                         crit = crit, folded = folded)
              if (is.finite(tau_sig_fit["tau_sig"])) {
                fit_res$es_tau_sig <- unname(tau_sig_fit["tau_sig"])
              }
            }
          }

          fit_res
  
      }, var_list = var_list, packages = "KernSmooth") # EOF boot_res

      if(is.null(boot_res)) stop("boostrap failure.")

    } # EOF OF

     # boot_res

      if (est_method == "EM") {

         var_list = c(var_list,"run_EM")

         boot_res <- bootstrap_lapply(function(i) {


         if (!is.null(cluster_id)) {
              sampled  <- sample(z_cluster_names, n_z_clusters, replace = TRUE)
              rows     <- unlist(z_cluster_rows[sampled], use.names = FALSE)
              boot_all <- test_stat[rows]
          } else {
             rows     <- sample(seq_along(test_stat), length(test_stat), replace = TRUE)
             boot_all <- test_stat[rows]
          }

          boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit_EM <- #tryCatch(
                 run_EM(
                 test_stat  = boot_all,
                 yi         = yi_boot,
                 sei        = sei_boot,
                 int_beg    = int_beg,
                 int_end    = int_end,
                 ncp        = ncp,
                 z_sd       = z_sd,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed    = w_fixed,
                 x_lim_min  = x_lim_min,
                 x_lim_max  = x_lim_max,
                 max_iter      = max_iter_boot,
                 tol_criterion = tol_criterion
                 )
	         #,error = function(e) NULL )
          if (!is.null(yi_boot)) {
            selection_fit <- fit_fixed_weight_selection(
              yi = yi_boot,
              sei = sei_boot,
              ODR = fit_EM$ODR,
              EDR = fit_EM$EDR,
              crit = crit,
              selection = "one_sided"
            )
            fit_EM$es_tau <- unname(selection_fit["tau"])
            if (use_folded_es_model && any(fit_EM$z_sd > 1 + 1e-10, na.rm = TRUE)) {
              effect_fit <- fit_folded_effect_model(
                yi_boot, sei_boot, fit_EM$ODR, fit_EM$EDR, crit
              )
              if (is.finite(effect_fit["mean_sig"])) {
                fit_EM$es_mean_sig <- unname(effect_fit["mean_sig"])
              }
              if (is.finite(effect_fit["tau_sig"])) {
                fit_EM$es_tau_sig <- unname(effect_fit["tau_sig"])
              }
            }

            if (use_fit_tau_sig) {
              tau_sig_fit <- fit_tau_sig(yi = yi_boot, sei = sei_boot,
                                         crit = crit, folded = folded)
              if (is.finite(tau_sig_fit["tau_sig"])) {
                fit_EM$es_tau_sig <- unname(tau_sig_fit["tau_sig"])
              }
            }
          }

          fit_EM

         }, var_list = var_list, packages = "stats") # close boot_res

    } # EOF EM 

    #boot_res 

    #est_method = "ML_gamma"; boot.iter = 500
    if (est_method == "ML_gamma") {

         n_starts = 2
         if (est_method == "ML_gamma") var_list = c(var_list,"run_ML_gamma","n_starts")

         boot_res <- bootstrap_lapply(function(i) {

           if (!is.null(cluster_id)) {
               sampled <- sample(z_cluster_names, n_z_clusters, replace = TRUE)
               boot_all <- unlist(z_clusters[sampled], use.names = FALSE)
               ODR_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]
           } else {
               boot_all <- sample(test_stat, length(test_stat), replace = TRUE)
               ODR_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]
           }


           gamma_fit <- tryCatch(
                 run_ml_gamma(boot_all, int_beg = int_beg, int_end = int_end, crit = crit,
                   n_starts = n_starts)
               ,error = function(e) NULL )

           fit_res = list(
             ODR    = ODR_boot,
             gamma_edr = if (!is.null(gamma_fit)) gamma_fit$gamma_edr[1] else NA,
             gamma_err = if (!is.null(gamma_fit)) gamma_fit$gamma_err[1] else NA
             #fit = gamma_fit$negloglik
             )

          fit_res
  
      }, var_list = var_list, packages = "KernSmooth") # EOF boot_res

   } # EOF ML_gamma 


   #summary(sapply(boot_res, function(x) x$gamma_edr))
   #summary(sapply(boot_res, function(x) x$gamma_err))

   #length(boot_res)
   if (length(boot_res) == 0)
      stop("All bootstrap runs failed.")

  return(boot_res)

}


##################################################
### End of New Extended Zcurve
##################################################

#zzz

run_zcurve = function(est_method,test_stat,crit,int_beg,int_end,
           ncp,z_sd,ncp_fixed,z_sd_fixed,w_fixed,
           cluster_id,boot_iter=0,edr_ci_adj=0,err_ci_adj=0,
           yi=NULL,sei=NULL,folded,directional,alpha,int_loc,df,
           x_lim_min,x_lim_max,spline_width,bw_est,augment,
           density_step,bootstrap_density_step,
           tol_criterion,max_iter,max_iter_boot,
           ci_alpha,parallel_boot,cores,seed) {


  zcurve_time = system.time({

  ci_probs <- c(ci_alpha / 2, 1 - ci_alpha / 2)

  int = test_stat[test_stat >= int_beg & test_stat <= int_end]

  ODR = mean(test_stat > crit)

  width = spline_width

  # ---- es_tau_sig floor correction (plug-in homogeneous null + blend) ----
  # estimate_es_tau_floor: es_tau_sig that a HOMOGENEOUS population (constant
  # effect = the fitted significant mean, observed sei, same fit) would report;
  # that residual is the floor. Distribution-free (no assumption on the effect
  # distribution -- it just measures what the mixture invents from se-spread).
  # The null is refit with the SAME engine as the main fit (run_EM vs run_OF):
  # EM spreads the ncp cloud more than OF at high-mean homogeneous data, so its
  # floor is larger, and an OF-estimated floor would under-correct EM.
  es_tau_floor <- NA_real_
  estimate_es_tau_floor <- function(mu, reps) {
    n <- length(sei); vals <- numeric(reps)
    for (r in seq_len(reps)) {
      y  <- mu + rnorm(n, 0, sei)              # one constant effect, no heterogeneity
      zz <- y / sei
      ft <- tryCatch(
        if (est_method == "EM") {
          run_EM(
            test_stat = abs(zz), yi = abs(y), sei = sei,
            int_beg = int_beg, int_end = int_end, ncp = ncp, z_sd = z_sd,
            ncp_fixed = ncp_fixed, z_sd_fixed = z_sd_fixed, w_fixed = w_fixed,
            x_lim_min = x_lim_min, x_lim_max = x_lim_max,
            max_iter = max_iter, tol_criterion = tol_criterion)
        } else {
          run_OF(
            test_stat = abs(zz), yi = abs(y), sei = sei,
            int_beg = int_beg, int_end = int_end, ncp = ncp, z_sd = z_sd,
            dens_pre = dens_pre, d_x_pre = d_x_pre, bw_est = bw_est, width = width,
            ncp_fixed = ncp_fixed, z_sd_fixed = z_sd_fixed, w_fixed = w_fixed,
            x_lim_min = x_lim_min, x_lim_max = x_lim_max, augment = augment,
            density_step = density_step, max_iter = max_iter,
            tol_criterion = tol_criterion)
        }, error = function(e) NULL)
      vals[r] <- if (is.null(ft)) NA_real_ else ft$es_tau_sig
    }
    mean(vals, na.rm = TRUE)
  }
  # variance-share blend: full floor removed at tau ~ floor, quartic taper
  # (reduction ~ floor^4 / tau^2) so genuine high-tau heterogeneity is intact.
  apply_floor_blend <- function(tau, floor) {
    if (!is.finite(tau) || !is.finite(floor) || floor <= 0 || tau <= 0) return(tau)
    if (tau <= floor) return(0)
    sqrt(max(0, tau^2 - floor^4 / tau^2))
  }

  d_x_pre  <- NULL
  dens_pre <- NULL

 if(est_method == "OF") { 

 ## ---- Precompute densities if both fixed ----
 if (ncp_fixed && z_sd_fixed) {
      densy     <- get_densities(int, bw = bw_est, from = int_beg,
                             to = int_end, width = width, augment = augment,
                             density_step = density_step)
      d_x_pre   <- densy[, 1]
      dens_pre  <- build_dens(d_x_pre, ncp, z_sd, df, 
           curve_type = curve_type, folded=folded)
      bar_width <- d_x_pre[2] - d_x_pre[1]
      dens_pre  <- dens_pre / (rowSums(dens_pre) * bar_width)
 }

         res_of <- run_OF(
                 test_stat  = test_stat, 
                 yi         = yi,
                 sei        = sei,
                 int_beg    = int_beg,
                 int_end    = int_end,   
                 ncp        = ncp,
                 z_sd       = z_sd,
                 dens_pre   = dens_pre,
                 d_x_pre    = d_x_pre,
                 bw_est     = bw_est,
                 width      = width,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed       = w_fixed,
                 x_lim_min     = x_lim_min,
                 x_lim_max     = x_lim_max,
                 augment       = augment,
                 density_step  = density_step,
                 max_iter      = max_iter,
                 tol_criterion = tol_criterion
           )

        res_of
        res_run = res_of 
  }  


   if (est_method == "EM") { 

         d_x_pre = NULL
         dens    = NULL

         res_EM <- run_EM(
                 test_stat  = test_stat, 
                 yi         = yi,
                 sei        = sei,
                 int_beg    = int_beg,
                 int_end    = int_end,   
                 ncp        = ncp,
                 z_sd       = z_sd,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed    = w_fixed,
                 x_lim_min  = x_lim_min,
                 x_lim_max  = x_lim_max,
                 max_iter      = max_iter,
                 tol_criterion = tol_criterion
        )


        res_EM
        res_run = res_EM
     
   }

    res_pe <- list(
      ODR              = ODR,
      EDR              = res_run$EDR,
      ERR              = res_run$ERR,
      ncp              = res_run$ncp,
      z_sd             = res_run$z_sd,
      w_inp            = res_run$w_inp,
      w_all            = res_run$w_all,
      local_power      = res_run$local_power,
      ncp_tau          = res_run$ncp_tau,
      model_fit        = res_run$model_fit,
      ODR_EDR_d        = ODR - res_run$EDR,
      shape_d_median   = res_run$shape_d_median,
      shape_d_mean     = res_run$shape_d_mean,
      local_es_pe      = res_run$local_es,
      local_tau_pe     = res_run$local_tau,
      es_mean_all_pe   = res_run$es_mean_all,
      es_median_all_pe = res_run$es_median_all,
      es_mean_sig_pe   = res_run$es_mean_sig,
      es_q25_sig_pe    = res_run$es_q25_sig,
      es_median_sig_pe = res_run$es_median_sig,
      es_q75_sig_pe    = res_run$es_q75_sig,
      es_tau_sig_pe    = res_run$es_tau_sig,
      es_tau_pe        = res_run$es_tau
     )

    if (!is.null(yi)) {
      selection_fit <- fit_fixed_weight_selection(
        yi = yi,
        sei = sei,
        ODR = res_pe$ODR,
        EDR = res_pe$EDR,
        crit = crit,
        selection = "one_sided"
      )
      res_pe$es_tau_pe <- unname(selection_fit["tau"])

      if (use_folded_es_model && any(res_pe$z_sd > 1 + 1e-10, na.rm = TRUE)) {
        effect_fit <- fit_folded_effect_model(
          yi = yi,
          sei = sei,
          ODR = res_pe$ODR,
          EDR = res_pe$EDR,
          crit = crit
        )
        if (is.finite(effect_fit["mean_sig"])) {
          res_pe$es_mean_sig_pe <- unname(effect_fit["mean_sig"])
        }
        if (is.finite(effect_fit["tau_sig"])) {
          res_pe$es_tau_sig_pe <- unname(effect_fit["tau_sig"])
        }
      }

      # Significant-results heterogeneity: study-level RE (off by default; the
      # distribution-free mixture es_tau_sig is used unless this is enabled).
      if (use_fit_tau_sig) {
        tau_sig_fit <- fit_tau_sig(yi = yi, sei = sei, crit = crit, folded = folded)
        if (is.finite(tau_sig_fit["tau_sig"])) {
          res_pe$es_tau_sig_pe <- unname(tau_sig_fit["tau_sig"])
        }
      }

      # Floor-correct es_tau_sig with the variance-share blend (mu = fitted
      # significant mean). Computed once here and reused for the bootstrap.
      if (floor_correct && isTRUE(folded) && is.finite(res_pe$es_mean_sig_pe) &&
          res_pe$es_mean_sig_pe > 0 && is.finite(res_pe$es_tau_sig_pe)) {
        es_tau_floor <- estimate_es_tau_floor(res_pe$es_mean_sig_pe, floor_reps)
        res_pe$es_tau_sig_pe <- apply_floor_blend(res_pe$es_tau_sig_pe, es_tau_floor)
      }
    }

    compute_power_h1_significant <- function(EDR, ERR) {
      FDR <- (1 / EDR - 1) * (alpha / (1 - alpha))
      power <- (ERR - FDR * (alpha / 2)) / (1 - FDR)
      pmin(1, pmax(0, power))
    }

    pow_h1_sig_pe <- compute_power_h1_significant(res_pe$EDR, res_pe$ERR)


  #################################

  ### finsh # res_pe

  #################################


  ### start boostrap

  local_es <- NULL
  local_tau <- NULL
  bootstrap_failures <- 0L
  bootstrap_workers <- 0L
  es_mean_all <- NULL
  es_median_all <- NULL
  es_mean_sig <- NULL
  es_q25_sig <- NULL
  es_median_sig <- NULL
  es_q75_sig <- NULL
  es_tau_sig <- NULL
  es_tau <- NULL
  pow_h1_sig <- NULL
  w_all_boot <- NULL
  ncp_boot <- NULL
  z_sd_boot <- NULL

  if (boot_iter > 0) {

       boot_dens_pre <- dens_pre
       boot_d_x_pre <- d_x_pre
       if (est_method == "OF" && ncp_fixed && z_sd_fixed) {
         boot_density <- get_densities(
           int,
           bw = bw_est,
           from = int_beg,
           to = int_end,
           width = width,
           augment = augment,
           density_step = bootstrap_density_step
         )
         boot_d_x_pre <- boot_density[, 1]
         boot_dens_pre <- build_dens(
           boot_d_x_pre, res_pe$ncp, res_pe$z_sd, df,
           curve_type = curve_type, folded = folded
         )
         boot_bar_width <- boot_d_x_pre[2] - boot_d_x_pre[1]
         boot_dens_pre <- boot_dens_pre /
           (rowSums(boot_dens_pre) * boot_bar_width)
       }

       boot_res <- run_bootstrap(
          test_stat     = test_stat, 
          yi            = yi,
          sei           = sei,
          crit          = crit,
          ncp           = res_pe$ncp,
          z_sd          = res_pe$z_sd,
          int_beg       = int_beg,
          int_end       = int_end,
          boot_iter     = boot_iter,
          directional   = directional,
          folded        = folded,
          dens_pre      = boot_dens_pre,
          d_x_pre       = boot_d_x_pre,
          ncp_fixed     = ncp_fixed,
          z_sd_fixed    = z_sd_fixed,
          w_fixed       = w_fixed,
          cluster_id    = cluster_id,
          est_method    = est_method,
          x_lim_min     = x_lim_min,
          x_lim_max     = x_lim_max,
          bw_est        = bw_est,
          augment       = augment,
          density_step  = bootstrap_density_step,
          tol_criterion = tol_criterion,
          max_iter_boot = max_iter_boot,
          parallel_boot = parallel_boot,
          cores         = cores,
          seed          = seed)

          bootstrap_failures <- attr(boot_res, "failed_replicates")
          bootstrap_workers <- attr(boot_res, "worker_count")
          if (is.null(bootstrap_failures)) bootstrap_failures <- 0L
          if (is.null(bootstrap_workers)) bootstrap_workers <- 1L

          # Per-replicate mixture, retained so individual effects can propagate
          # the population (w_all) uncertainty into each study's posterior.
          w_all_boot <- sapply(boot_res, function(x) x$w_all)
          ncp_boot   <- sapply(boot_res, function(x) x$ncp)
          z_sd_boot  <- sapply(boot_res, function(x) x$z_sd)

          ODR_boot   <- sapply(boot_res, function(x) x$ODR)
          ODR_boot[is.nan(ODR_boot)] = NA
          ODR_ci = c(
             quantile(ODR_boot, probs = ci_probs[1], na.rm = TRUE),
             quantile(ODR_boot, probs = ci_probs[2], na.rm = TRUE)
          )


          EDR_boot   <- sapply(boot_res, function(x) x$EDR)

          EDR_boot_lb <- pmax(alpha, EDR_boot - edr_ci_adj)
          EDR_boot_ub <- pmin(1,     EDR_boot + edr_ci_adj)

          EDR_ci = c(
             quantile(EDR_boot_lb, probs = ci_probs[1], na.rm = TRUE),
             quantile(EDR_boot_ub, probs = ci_probs[2], na.rm = TRUE)
          )


          ERR_boot   <- sapply(boot_res, function(x) x$ERR)

          ERR_boot_lb <- pmax(alpha, ERR_boot - err_ci_adj)
          ERR_boot_ub <- pmin(1,     ERR_boot + err_ci_adj)

          ERR_ci = c(
             quantile(ERR_boot_lb, probs = ci_probs[1], na.rm = TRUE),
             quantile(ERR_boot_ub, probs = ci_probs[2], na.rm = TRUE)
          )

          pow_h1_sig_boot <- compute_power_h1_significant(EDR_boot, ERR_boot)
          pow_h1_sig_ci <- quantile(
            pow_h1_sig_boot,
            probs = ci_probs,
            na.rm = TRUE
          )
          pow_h1_sig <- c(pow_h1_sig_pe, pow_h1_sig_ci)


          # Each endpoint is a one-sided (1 - ci_alpha) bound. The lower
          # bound is the publication-bias test and therefore subtracts the
          # adjusted upper EDR, not the lower EDR.
          ODR_EDR_d_ci  <- c(
               quantile(ODR_boot - EDR_boot_ub, ci_alpha, na.rm = TRUE),
               quantile(ODR_boot - EDR_boot_lb, 1 - ci_alpha, na.rm = TRUE)
           )
          ODR_EDR_d = c(res_pe$ODR_EDR_d,ODR_EDR_d_ci)

          if(length(res_pe$ncp) == 1) {
             ncp_ci = quantile(sapply(boot_res, function(x) x$ncp), ci_probs)
             ncp = c(res_pe$ncp,ncp_ci)
          } else {
             ncp_ci = apply(sapply(boot_res, function(x) x$ncp),1,function(x) quantile(x,ci_probs) )
             ncp  = rbind(res_pe$ncp, ncp_ci) 
          }

          ncp_tau_boot <- sapply(boot_res, function(x) x$ncp_tau)
          ncp_tau_ci <- if (any(is.finite(ncp_tau_boot))) {
            quantile(ncp_tau_boot, ci_probs, na.rm = TRUE)
          } else {
            c(NA_real_, NA_real_)
          }
          ncp_tau <- c(res_pe$ncp_tau, ncp_tau_ci)

          if(length(res_pe$z_sd) == 1) {
             z_sd_ci = quantile(sapply(boot_res, function(x) x$z_sd), ci_probs)
             z_sd = c(res_pe$z_sd,z_sd_ci)
          } else {
             z_sd_ci = apply(sapply(boot_res, function(x) x$z_sd),1,function(x) quantile(x,ci_probs) )
             z_sd  = rbind(res_pe$z_sd, z_sd_ci) 
          }

          if(length(res_pe$w_all) == 1) {
             w_all_ci = quantile(sapply(boot_res, function(x) x$w_all), ci_probs)
             w_all = c(res_pe$w_all,w_all_ci)
          } else {
             w_all_ci = apply(sapply(boot_res, function(x) x$w_all),1,function(x) quantile(x,ci_probs) )
             w_all  = rbind(res_pe$w_all, w_all_ci) 
          }

          if(length(res_pe$w_inp) == 1) {
             w_inp_ci = quantile(sapply(boot_res, function(x) x$w_inp), ci_probs)
             w_inp = c(res_pe$w_inp,w_inp_ci)
          } else {
             w_inp_ci = apply(sapply(boot_res, function(x) x$w_inp),1,function(x) quantile(x,ci_probs) )
             w_inp  = rbind(res_pe$w_inp, w_inp_ci) 
          }

          if(length(res_pe$local_power) == 1) {
             local_power_ci = quantile(sapply(boot_res, function(x) x$local_power), ci_probs)
             local_power = c(res_pe$local_power,local_power_ci)
          } else {
             local_power_ci = apply(sapply(boot_res, function(x) x$local_power),1,function(x) quantile(x,ci_probs) )
             local_power  = rbind(res_pe$local_power, local_power_ci) 
          }


          shape_d_median_boot <- sapply(boot_res, function(x) x$shape_d_median)
          shape_d_mean_boot <- sapply(boot_res, function(x) x$shape_d_mean)

          shape_d_median_ci <- if (any(is.finite(shape_d_median_boot))) {
            quantile(shape_d_median_boot, ci_probs, na.rm = TRUE)
          } else {
            c(NA_real_, NA_real_)
          }
          shape_d_mean_ci <- if (any(is.finite(shape_d_mean_boot))) {
            quantile(shape_d_mean_boot, ci_probs, na.rm = TRUE)
          } else {
            c(NA_real_, NA_real_)
          }

          shape_d_median <- c(res_pe$shape_d_median, shape_d_median_ci)
          shape_d_mean <- c(res_pe$shape_d_mean, shape_d_mean_ci)


          ODR = c(res_pe$ODR,ODR_ci)
          EDR = c(res_pe$EDR,EDR_ci)
          ERR = c(res_pe$ERR,ERR_ci)

          EDR[2] = max(alpha,EDR[2])
          EDR[3] = min(1,EDR[3])

          ERR[2] = max(alpha/2,ERR[2])
          ERR[3] = min(1,ERR[3])

          ###

          if(!is.null(yi)) {

            es_mean_all_boot = sapply(boot_res, function(x) x$es_mean_all)

            es_mean_all_ci = quantile(es_mean_all_boot, probs = ci_probs, na.rm=TRUE)
            es_mean_all = c(res_pe$es_mean_all_pe,es_mean_all_ci)

            es_median_all_boot = sapply(boot_res, function(x) x$es_median_all)
            es_median_all_ci = quantile(es_median_all_boot, probs = ci_probs, na.rm=TRUE)
            es_median_all = c(res_pe$es_median_all_pe,es_median_all_ci)

            es_tau_boot = sapply(boot_res, function(x) x$es_tau)
            es_tau_ci = quantile(es_tau_boot, probs = ci_probs, na.rm=TRUE)
            es_tau = c(res_pe$es_tau_pe,es_tau_ci)


            es_mean_sig_boot = sapply(boot_res, function(x) x$es_mean_sig)

            es_mean_sig_ci = quantile(es_mean_sig_boot, probs = ci_probs, na.rm=TRUE)
            es_mean_sig = c(res_pe$es_mean_sig_pe,es_mean_sig_ci)

            es_q25_sig_boot = sapply(boot_res, function(x) x$es_q25_sig)
            es_q25_sig_ci = quantile(es_q25_sig_boot, probs = ci_probs, na.rm=TRUE)
            es_q25_sig = c(res_pe$es_q25_sig_pe, es_q25_sig_ci)

            es_median_sig_boot = sapply(boot_res, function(x) x$es_median_sig)
            es_median_sig_ci = quantile(es_median_sig_boot, probs = ci_probs, na.rm=TRUE)
            es_median_sig = c(res_pe$es_median_sig_pe, es_median_sig_ci)

            es_q75_sig_boot = sapply(boot_res, function(x) x$es_q75_sig)
            es_q75_sig_ci = quantile(es_q75_sig_boot, probs = ci_probs, na.rm=TRUE)
            es_q75_sig = c(res_pe$es_q75_sig_pe, es_q75_sig_ci)

            es_tau_sig_boot = sapply(boot_res, function(x) x$es_tau_sig)
            if (floor_correct && is.finite(es_tau_floor)) {
              es_tau_sig_boot = vapply(es_tau_sig_boot,
                function(t) apply_floor_blend(t, es_tau_floor), numeric(1))
            }
            es_tau_sig_ci = quantile(es_tau_sig_boot, probs = ci_probs, na.rm=TRUE)
            es_tau_sig = c(res_pe$es_tau_sig_pe, es_tau_sig_ci)

            es_loc = sapply(boot_res, function(x) x$local_es)

            if(length(res_pe$local_es_pe) == 1) {
               local_es_ci = quantile(es_loc, ci_probs,na.rm=TRUE)
               local_es = c(res_pe$local_es_pe,local_es_ci)
            } else {
               local_es_ci = apply(es_loc,1,function(x) quantile(x,ci_probs,na.rm=TRUE) )
               local_es  = rbind(res_pe$local_es_pe, local_es_ci) 
            }

            local_tau_boot = sapply(boot_res, function(x) x$local_tau)
            if(length(res_pe$local_tau_pe) == 1) {
               local_tau_ci = quantile(local_tau_boot, ci_probs, na.rm=TRUE)
               local_tau = c(res_pe$local_tau_pe, local_tau_ci)
            } else {
               local_tau_ci = apply(
                 local_tau_boot, 1,
                 function(x) quantile(x, ci_probs, na.rm=TRUE)
               )
               local_tau = rbind(res_pe$local_tau_pe, local_tau_ci)
            }


          }# EOF if yi 


   } else {
        ODR = t(c(res_pe$ODR,rep(NA,2)))
        EDR = t(c(res_pe$EDR,rep(NA,2)))
        ERR = t(c(res_pe$ERR,rep(NA,2)))
        ncp = t(cbind(res_pe$ncp, matrix(NA,length(res_pe$w_all),2)))
        ncp_tau = c(res_pe$ncp_tau, NA, NA)
        z_sd = t(cbind(res_pe$z_sd, matrix(NA,length(res_pe$w_inp),2)))
        w_all = t(cbind(res_pe$w_all, matrix(NA,length(res_pe$w_all),2)))
        w_inp = t(cbind(res_pe$w_inp, matrix(NA,length(res_pe$w_inp),2)))
        local_power = t(cbind(res_pe$local_power,matrix(NA,length(res_pe$local_power),2)))
        ODR_EDR_d = c(res_pe$ODR_EDR_d,NA,NA)
        shape_d_median  = c(res_pe$shape_d_median,NA,NA)
        shape_d_mean    = c(res_pe$shape_d_mean,NA,NA)
        pow_h1_sig = c(pow_h1_sig_pe, NA, NA)

        if (!is.null(yi)) {
          local_es = t(cbind(res_pe$local_es_pe,matrix(NA,length(res_pe$local_es_pe),2)))
          local_tau = t(cbind(res_pe$local_tau_pe,matrix(NA,length(res_pe$local_tau_pe),2)))
          es_mean_all = c(res_pe$es_mean_all_pe,NA,NA)
          es_median_all = c(res_pe$es_median_all_pe,NA,NA)
          es_mean_sig = c(res_pe$es_mean_sig_pe,NA,NA)
          es_q25_sig = c(res_pe$es_q25_sig_pe,NA,NA)
          es_median_sig = c(res_pe$es_median_sig_pe,NA,NA)
          es_q75_sig = c(res_pe$es_q75_sig_pe,NA,NA)
          es_tau_sig = c(res_pe$es_tau_sig_pe,NA,NA)
          es_tau = c(res_pe$es_tau_pe,NA,NA)
        }

  } # EOF no bootstrap

}) # EOF Timer

  zcurve_res = list(
            ODR                = ODR,
            EDR                = EDR,
            ERR                = ERR,
            pow_h1_sig         = pow_h1_sig,
            ncp                = ncp,
            ncp_tau            = ncp_tau,
            z_sd               = z_sd,
            w_inp              = w_inp,
            w_all              = w_all,
            w_all_boot         = w_all_boot,
            ncp_boot           = ncp_boot,
            z_sd_boot          = z_sd_boot,
            local_power        = local_power,
            local_es           = local_es,
            local_tau          = local_tau,
            ODR_EDR_d          = ODR_EDR_d,
            shape_d_median     = shape_d_median,
            shape_d_mean       = shape_d_mean,
            es_mean_all        = es_mean_all,
            es_median_all      = es_median_all,  
            es_mean_sig        = es_mean_sig,
            es_q25_sig         = es_q25_sig,
            es_median_sig      = es_median_sig,
            es_q75_sig         = es_q75_sig,
            es_tau_sig         = es_tau_sig,
            es_tau             = es_tau,
            bootstrap_failures = bootstrap_failures,
            bootstrap_workers  = bootstrap_workers,
            time               = zcurve_time
          )



#zcurve.res

#print(zcurve.res)

return(zcurve_res)

}


#################################################
#################################################
#################################################

#For quick testing
#test_stat = abs(c(rnorm(500,0,1),rnorm(100,2,1)));hist(test_stat)
#cluster.id = NULL;yi=NULL;sei=NULL


#################################################
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
#################################################


  ### --- curve type / crit derivation ---
  curve_type <- "z"
  # crit already set via the formal-parameter default (qnorm(1 - alpha/2)) or user override

  ### --- determine test_stat from whichever input was supplied ---
  if (xor(is.null(yi), is.null(sei))) {
    stop("yi and sei must be supplied together.")
  }
  if (!is.null(yi) && length(yi) != length(sei)) {
    stop("yi and sei must have the same length.")
  }

  n_test_inputs <- sum(!is.null(zval), !is.null(tval), !is.null(pval))
  if (n_test_inputs == 0 && is.null(yi)) {
    stop("Supply one of zval, tval (+ df), pval, or yi (+ sei).")
  }
  if (n_test_inputs > 1) {
    stop("Supply only one of zval, tval (+ df), or pval.")
  }

  if (!is.null(zval)) test_stat <- zval
  if (!is.null(pval)) test_stat <- qnorm(1 - pval / 2)
  if (n_test_inputs == 0 && !is.null(yi)) test_stat <- yi / sei

  input_type <- if (!is.null(zval)) {
    if (is.null(yi)) "z_values" else "z_values_with_effect_sizes"
  } else if (!is.null(tval)) {
    if (is.null(yi)) "t_values" else "t_values_with_effect_sizes"
  } else if (!is.null(pval)) {
    if (is.null(yi)) "p_values" else "p_values_with_effect_sizes"
  } else {
    "effect_sizes"
  }

  if (!is.null(tval)) {
    if (is.null(df)) stop("df must be supplied when tval is used.")
    curve_type <- "t"
    if (missing(crit)) crit <- qt(alpha / 2, df, lower.tail = FALSE)
    test_stat <- tval
  }

  # Nondirectional (folded) analyses model |z|: fold the input so BOTH
  # significant tails contribute to the fit and the plotted density (which is
  # symmetric about 0) sits over the |z| data instead of signed data. Already-
  # folded input (|z|) is unchanged. Directional analyses keep the sign.
  if (isTRUE(folded)) {
    test_stat <- abs(test_stat)
    if (!is.null(yi)) yi <- abs(yi)
  }

  if (!is.null(yi) && length(test_stat) != length(yi)) {
    stop("The test-statistic input, yi, and sei must have the same length.")
  }
  if (!is.null(study_id) && length(study_id) != length(test_stat)) {
    stop("study_id must have the same length as the test-statistic input.")
  }
  if (!is.null(cluster_id) && length(cluster_id) != length(test_stat)) {
    stop("cluster_id must have the same length as the test-statistic input.")
  }
  if (is.null(study_id)) study_id <- seq_along(test_stat)

  if (curve_type == "t" && est_method != "OF") {
    stop("Genuine t-curve fitting (tval + df) is only supported for est_method = \"OF\".")
  }

  ### --- remove NA (one mask, built from test_stat, not from zval/tval specifically) ---
  sel <- !is.na(test_stat)
  if (!is.null(cluster_id)) sel <- sel & !is.na(cluster_id)
  if (!is.null(yi))         sel <- sel & !is.na(yi)
  if (!is.null(sei))        sel <- sel & !is.na(sei)

  test_stat  <- test_stat[sel]
  if (!is.null(cluster_id)) cluster_id <- cluster_id[sel]
  study_id  <- study_id[sel]
  if (!is.null(yi))         yi         <- yi[sel]
  if (!is.null(sei))        sei        <- sei[sel]

  yi_observed <- yi

  ### --- fold to |test_stat| when sign doesn't matter ---
  if (!directional) {
    test_stat <- abs(test_stat)
    if (!is.null(yi)) yi <- abs(yi)
  }

  n_tests <- length(test_stat)
  n_sig   <- sum(test_stat > crit)
  if (n_sig < 20) stop("Insufficient data: need at least 20 significant test statistics.")

  ### --- extreme-value diagnostic ---
  extreme <- c(
    ext_neg = sum(test_stat < -int_end) / max(1, sum(test_stat < -crit)),
    ext_pos = sum(test_stat >  int_end) / max(1, sum(test_stat >  crit))
  )

  n_clusters <- NULL
  if (!is.null(cluster_id)) n_clusters <- length(unique(cluster_id))

  ### --- density-fit diagnostic (informational only -- not fed into any adjustment) ---
  slope_width <- 1
  slope_int   <- test_stat[test_stat >= int_beg + 2 * control$bw_est &
                            test_stat <  int_beg + 2 * control$bw_est + slope_width]
  slope_k <- length(slope_int)
  slope   <- rep(NA, 6)
  if (slope_k > 10) {
    s <- get_slope(int = slope_int, bw = control$bw_est, 
         from = int_beg, to = int_beg + slope_width)
    reliability    <- s["slope"]^2 / (s["slope"]^2 + s["se"]^2)
    adjusted_slope <- reliability * s["slope"]
    slope <- c(slope_k, adjusted_slope, s)
  }
  names(slope) <- c("slope_k", "slope_adj", "slope_pe", "slope_se", "slope_lb", "slope_ub")

  ### --- fit ---
  res_main <- run_zcurve(
    est_method             = est_method,
    test_stat              = test_stat,
    crit                   = crit,
    int_beg                = int_beg,
    int_end                = int_end,
    cluster_id             = cluster_id,
    boot_iter              = control$boot_iter,
    ncp                    = ncp,
    z_sd                   = z_sd,
    ncp_fixed              = ncp_fixed,
    z_sd_fixed             = z_sd_fixed,
    w_fixed                = w_fixed,
    augment                = control$augment,
    bw_est                 = control$bw_est,
    directional            = directional,   # NOTE: run_zcurve() doesn't accept this yet -- see file header
    folded                 = folded,
    alpha                  = alpha,         # NOTE: same
    int_loc                = int_loc,       # NOTE: same
    df                     = df,            # NOTE: same
    yi                     = yi,
    sei                    = sei,
    x_lim_min              = x_lim_min,
    x_lim_max              = x_lim_max,
    spline_width           = control$spline_width,
    density_step           = control$density_step,
    bootstrap_density_step = control$bootstrap_density_step,
    tol_criterion          = control$tol_criterion,
    max_iter               = control$max_iter,
    max_iter_boot          = control$max_iter_boot,
    ci_alpha               = control$ci_alpha,
    parallel_boot          = control$parallel,
    cores                  = control$cores,
    seed                   = control$seed
  )

  ### --- assemble output ---
  FDR <- as.numeric((1 / res_main$EDR - 1) * (alpha / (1 - alpha)))
  FDR <- FDR[c(1, 3, 2)]
  names(FDR) <- c("FDR_pe", "FDR_lb", "FDR_ub")

  pow_h1_sig <- as.numeric(res_main$pow_h1_sig)
  names(pow_h1_sig) <- c("pow_h1_sig_pe", "pow_h1_sig_lb", "pow_h1_sig_ub")

  ODR_EDR_d <- as.numeric(res_main$ODR_EDR_d)
  names(ODR_EDR_d) <- c("ODR_EDR_d_pe", "ODR_EDR_d_lb", "ODR_EDR_d_ub")

  res_main$ODR <- as.numeric(res_main$ODR); names(res_main$ODR) <- c("ODR_pe", "ODR_lb", "ODR_ub")
  res_main$EDR <- as.numeric(res_main$EDR); names(res_main$EDR) <- c("EDR_pe", "EDR_lb", "EDR_ub")
  res_main$ERR <- as.numeric(res_main$ERR); names(res_main$ERR) <- c("ERR_pe", "ERR_lb", "ERR_ub")
  res_main$ncp_tau <- as.numeric(res_main$ncp_tau)
  names(res_main$ncp_tau) <- c("ncp_tau_pe", "ncp_tau_lb", "ncp_tau_ub")

  names(res_main$shape_d_mean) <- c(
    "shape_d_mean_pe", "shape_d_mean_lb", "shape_d_mean_ub"
  )
  names(res_main$shape_d_median) <- c(
    "shape_d_median_pe", "shape_d_median_lb", "shape_d_median_ub"
  )

  local_interval_labels <- NULL
  if (is.finite(int_loc) && int_loc > 0) {
    local_breaks <- seq(x_lim_min, x_lim_max, by = int_loc)
    local_interval_labels <- paste0(
      head(local_breaks, -1), "-", tail(local_breaks, -1)
    )
  }

  if (is.matrix(res_main$local_power)) {
    if (length(local_interval_labels) == ncol(res_main$local_power)) {
      colnames(res_main$local_power) <- local_interval_labels
    } else if (is.null(colnames(res_main$local_power))) {
      colnames(res_main$local_power) <- paste0("local_power", seq_len(ncol(res_main$local_power)))
    }
    rownames(res_main$local_power) <- c("pe", "lb", "ub")
  }

  if (k_ncp > 1) {
    colnames(res_main$ncp)   <- paste0("ncp",   seq_len(k_ncp)); rownames(res_main$ncp)   <- c("pe", "lb", "ub")
    colnames(res_main$z_sd)  <- paste0("z_sd",  seq_len(k_ncp)); rownames(res_main$z_sd)  <- c("pe", "lb", "ub")
    colnames(res_main$w_inp) <- paste0("w_inp", seq_len(k_ncp)); rownames(res_main$w_inp) <- c("pe", "lb", "ub")
    colnames(res_main$w_all) <- paste0("w_all", seq_len(k_ncp)); rownames(res_main$w_all) <- c("pe", "lb", "ub")
  }

  point_values <- function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      return(as.numeric(x[1, , drop = TRUE]))
    }
    as.numeric(x)
  }

  plot_values <- test_stat[test_stat >= x_lim_min & test_stat <= x_lim_max]
  observed_density <- get_densities(
    plot_values,
    bw = control$bw_est,
    from = x_lim_min,
    to = x_lim_max,
    width = control$spline_width,
    augment = control$augment,
    density_step = control$density_step
  )
  plot_ncp <- point_values(res_main$ncp)
  plot_z_sd <- point_values(res_main$z_sd)
  plot_weights <- point_values(res_main$w_all)
  component_density <- build_dens(
    d_x = observed_density$zx,
    ncp = plot_ncp,
    z_sd = plot_z_sd,
    df = df,
    curve_type = curve_type,
    folded = folded
  )
  plot_bin_width <- observed_density$zx[2] - observed_density$zx[1]
  component_density <- component_density /
    (rowSums(component_density) * plot_bin_width)
  fitted_density <- colSums(component_density * plot_weights)

  plot_data <- list(
    values = plot_values,
    grid = observed_density$zx,
    observed_density = observed_density$zy,
    fitted_density = fitted_density,
    interval = c(lower = x_lim_min, upper = x_lim_max),
    selection_interval = c(lower = int_beg, upper = int_end),
    critical_value = crit,
    folded = folded,
    directional = directional
  )

  individual_effects <- NULL
  individual_posterior <- NULL
  if (!is.null(yi)) {
    component_order <- order(plot_ncp)
    component_ncp <- plot_ncp[component_order]
    component_weight <- plot_weights[component_order]
    component_z_sd <- plot_z_sd[component_order]

    # A fitted continuous component already supplies its latent-NCP
    # distribution through quadrature. Discrete mixtures retain the existing
    # interpolation used to avoid artificial zero-width posterior intervals.
    if (any(component_z_sd > 1 + 1e-10)) {
      individual_components <- expand_ncp_components(
        ncp = component_ncp,
        weight = component_weight,
        z_sd = component_z_sd
      )
      individual_order <- order(individual_components$ncp)
      individual_ncp <- individual_components$ncp[individual_order]
      individual_prior <- individual_components$weight[individual_order]
      individual_z_sd <- individual_components$z_sd[individual_order]
    } else if (length(component_ncp) > 1L && diff(range(component_ncp)) > 0) {
      individual_ncp <- seq(
        min(component_ncp), max(component_ncp), by = .01
      )
      if (tail(individual_ncp, 1) < max(component_ncp)) {
        individual_ncp <- c(individual_ncp, max(component_ncp))
      }
      individual_prior <- approx(
        component_ncp, component_weight,
        xout = individual_ncp,
        method = "linear", rule = 2,
        ties = "ordered"
      )$y
      individual_z_sd <- approx(
        component_ncp, component_z_sd,
        xout = individual_ncp,
        method = "linear", rule = 2,
        ties = "ordered"
      )$y
    } else {
      individual_ncp <- component_ncp
      individual_prior <- component_weight
      individual_z_sd <- component_z_sd
    }
    individual_prior <- pmax(individual_prior, 0)
    individual_prior <- individual_prior / sum(individual_prior)

    # Rebuild the prior on the fixed latent-NCP grid from ANY set of component
    # weights, so the per-replicate bootstrap w_all can be pushed through the
    # same posterior machinery (grid and likelihood do not depend on weights).
    if (any(component_z_sd > 1 + 1e-10)) {
      .expanded  <- expand_ncp_components(ncp = component_ncp,
                                          weight = component_weight, z_sd = component_z_sd)
      .eo        <- order(.expanded$ncp)
      .node_comp <- .expanded$component_id[.eo]
      .node_qw   <- .expanded$quadrature_weight[.eo]
      prior_builder <- function(cw) {
        p <- pmax(cw[.node_comp] * .node_qw, 0); s <- sum(p); if (s <= 0) p else p / s
      }
    } else if (length(component_ncp) > 1L && diff(range(component_ncp)) > 0) {
      prior_builder <- function(cw) {
        p <- pmax(approx(component_ncp, cw, xout = individual_ncp,
                         method = "linear", rule = 2, ties = "ordered")$y, 0)
        s <- sum(p); if (s <= 0) p else p / s
      }
    } else {
      prior_builder <- function(cw) {
        p <- pmax(cw, 0); s <- sum(p); if (s <= 0) p else p / s
      }
    }

    individual_likelihood <- outer(
      test_stat, seq_along(individual_ncp),
      function(observed_z, component) {
        mu <- individual_ncp[component]
        component_sd <- individual_z_sd[component]
        if (folded) {
          dnorm(observed_z, mean = mu, sd = component_sd) +
            dnorm(-observed_z, mean = mu, sd = component_sd)
        } else {
          dnorm(observed_z, mean = mu, sd = component_sd)
        }
      }
    )
    latent_ncp <- if (folded || !directional) abs(individual_ncp) else individual_ncp
    individual_theta <- outer(sei, latent_ncp, "*")

    # Null-share smoothing of the INDIVIDUAL-effects prior only. Where the low-
    # power region (ncp <= 1) carries weight, give the null (ncp 0) at least
    # `es_null_share` of it, taken proportionally from the other low-power rungs,
    # so a genuinely-null study can shrink to 0. Self-limiting: with no low-power
    # mass (a high-power fit) it does nothing. Applied to w_all here; the fit's
    # w_inp / EDR / ERR are computed upstream and never see this.
    es_null_share <- if (is.null(control$es_null_share)) 0 else control$es_null_share
    lp_idx   <- which(component_ncp <= 1 + 1e-9)
    null_idx <- which.min(component_ncp)
    smooth_cw <- function(cw) {
      if (es_null_share <= 0 || length(lp_idx) < 2L) return(cw)
      W_lp <- sum(cw[lp_idx])
      target <- es_null_share * W_lp
      if (W_lp <= 0 || cw[null_idx] >= target) return(cw)
      add <- target - cw[null_idx]
      donors <- setdiff(lp_idx, null_idx)
      dW <- sum(cw[donors])
      if (dW <= 0) return(cw)
      cw[donors] <- cw[donors] * (1 - add / dW)
      cw[null_idx] <- cw[null_idx] + add
      cw
    }

    posterior_from_weights <- function(cw) {
      pp <- sweep(individual_likelihood, 2, prior_builder(smooth_cw(cw)), "*")
      pp <- sweep(pp, 1, rowSums(pp), "/")
      pp[!is.finite(pp)] <- NA_real_
      pp
    }

    posterior_quantile <- function(values, weights, probability) {
      keep <- is.finite(values) & is.finite(weights) & weights > 0
      if (!any(keep)) return(NA_real_)
      values <- values[keep]
      weights <- weights[keep]
      ord <- order(values)
      values <- values[ord]
      cumulative <- cumsum(weights[ord]) / sum(weights)
      values[which(cumulative >= probability)[1]]
    }

    # Adjusted individual effects REQUIRE the bootstrap. An honest interval must
    # reflect uncertainty in the fitted population (w_all) -- not just the study's
    # own se -- so each study's posterior is POOLED over the per-replicate w_all.
    # For an imprecise study this widens the interval by the population/mean
    # uncertainty (the (1 - B)^2 Var(mu) term). Without a bootstrap the adjusted
    # columns are left NA and the forest plot refuses to run.
    w_all_boot <- res_main$w_all_boot
    have_boot  <- is.matrix(w_all_boot) && ncol(w_all_boot) >= 2L &&
      nrow(w_all_boot) == length(component_order)

    if (have_boot) {
      pooled <- matrix(0, nrow = n_tests, ncol = length(individual_ncp))
      n_ok <- 0L
      for (b in seq_len(ncol(w_all_boot))) {
        cw_b <- w_all_boot[component_order, b]
        if (any(!is.finite(cw_b)) || sum(pmax(cw_b, 0)) <= 0) next
        pb <- posterior_from_weights(cw_b)
        pb[!is.finite(pb)] <- 0
        pooled <- pooled + pb
        n_ok <- n_ok + 1L
      }
      if (n_ok > 0L) pooled <- pooled / n_ok
      individual_posterior_prob <- pooled

      adjusted_effect <- rowSums(individual_posterior_prob * individual_theta, na.rm = TRUE)
      adjusted_second <- rowSums(individual_posterior_prob * individual_theta^2, na.rm = TRUE)
      adjusted_sd <- sqrt(pmax(adjusted_second - adjusted_effect^2, 0))
      adjusted_lb <- vapply(seq_len(n_tests), function(i) posterior_quantile(
        individual_theta[i, ], individual_posterior_prob[i, ], control$ci_alpha / 2), numeric(1))
      adjusted_median <- vapply(seq_len(n_tests), function(i) posterior_quantile(
        individual_theta[i, ], individual_posterior_prob[i, ], .5), numeric(1))
      adjusted_ub <- vapply(seq_len(n_tests), function(i) posterior_quantile(
        individual_theta[i, ], individual_posterior_prob[i, ], 1 - control$ci_alpha / 2), numeric(1))
      prob_meaningful <- rowSums(
        individual_posterior_prob * (abs(individual_theta) > min_effect), na.rm = TRUE)
      # Adjusted power = the study's replication probability: the power at each
      # latent ncp, averaged over the study's pooled posterior.
      comp_power <- pnorm((latent_ncp - crit) / individual_z_sd) +
        pnorm((-crit - latent_ncp) / individual_z_sd)
      power_adjusted <- as.numeric(individual_posterior_prob %*% comp_power)
    } else {
      individual_posterior_prob <- posterior_from_weights(component_weight)
      adjusted_effect <- rep(NA_real_, n_tests)
      adjusted_sd     <- rep(NA_real_, n_tests)
      adjusted_lb     <- rep(NA_real_, n_tests)
      adjusted_median <- rep(NA_real_, n_tests)
      adjusted_ub     <- rep(NA_real_, n_tests)
      prob_meaningful <- rep(NA_real_, n_tests)
      power_adjusted  <- rep(NA_real_, n_tests)
    }

    observed_effect <- if (folded || !directional) abs(yi_observed) else yi_observed
    sampling_critical <- qnorm(1 - control$ci_alpha / 2)
    if (folded || !directional) {
      observed_lb <- pmax(0, observed_effect - sampling_critical * sei)
      observed_ub <- observed_effect + sampling_critical * sei
    } else {
      observed_lb <- observed_effect - sampling_critical * sei
      observed_ub <- observed_effect + sampling_critical * sei
    }

    individual_effects <- data.frame(
      study_id = study_id,
      cluster_id = if (is.null(cluster_id)) rep(NA, length(study_id)) else cluster_id,
      z = test_stat,
      significant = test_stat > crit,
      effect_observed = observed_effect,
      effect_observed_lb = observed_lb,
      effect_observed_ub = observed_ub,
      effect_adjusted = adjusted_median,
      effect_adjusted_mean = adjusted_effect,
      effect_adjusted_median = adjusted_median,
      effect_adjusted_sd = adjusted_sd,
      effect_adjusted_lb = adjusted_lb,
      effect_adjusted_ub = adjusted_ub,
      min_effect = min_effect,
      prob_meaningful = prob_meaningful,
      power_adjusted = power_adjusted,
      sei = sei,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    individual_posterior <- list(
      probabilities = individual_posterior_prob,
      component_ncp = latent_ncp,
      sei = sei
    )
  }

  results <- list(
    call           = zcurve_call,
    metadata       = list(
      method = est_method,
      input = input_type,
      n_tests = n_tests,
      n_significant = n_sig,
      n_clusters = n_clusters,
      alpha = alpha,
      critical_value = crit,
      interval = c(lower = int_beg, upper = int_end),
      bootstrap_iterations = control$boot_iter,
      bootstrap_failures = res_main$bootstrap_failures,
      bootstrap_parallel = res_main$bootstrap_workers > 1L,
      bootstrap_workers = res_main$bootstrap_workers,
      bootstrap_seed = control$seed,
      confidence_level = 1 - control$ci_alpha,
      min_effect = min_effect,
      effect_scale = if (folded || !directional) "absolute" else "signed",
      version = version,
      date = date,
      elapsed_seconds = unname(res_main$time[["elapsed"]])
    ),
    clusters       = n_clusters,
    slope          = slope,
    extreme        = extreme,
    ODR            = res_main$ODR,
    EDR            = res_main$EDR,
    ERR            = res_main$ERR,
    FDR            = FDR,
    pow_h1_sig     = pow_h1_sig,
    ncp            = res_main$ncp,
    ncp_tau        = res_main$ncp_tau,
    z_sd           = res_main$z_sd,
    w_inp          = res_main$w_inp,
    w_all          = res_main$w_all,
    w_all_boot     = res_main$w_all_boot,
    ncp_boot       = res_main$ncp_boot,
    z_sd_boot      = res_main$z_sd_boot,
    local_power    = res_main$local_power,
    individual_effects = individual_effects,
    individual_posterior = individual_posterior,
    val_input      = test_stat,
    plot_data      = plot_data,
    ODR_EDR_d      = ODR_EDR_d,
    shape_d_median = res_main$shape_d_median,
    shape_d_mean   = res_main$shape_d_mean,
    fit            = res_main$model_fit

#   version        = if (is.null(version) || version == "") NULL else version,
#   date           = if (is.null(date)    || date    == "") NULL else date
  )

  if (!is.null(res_main$gamma_shape)) {
    results$gamma <- list(
      shape = res_main$gamma_shape,
      rate = res_main$gamma_rate,
      EDR = res_main$gamma_edr,
      ERR = res_main$gamma_err
    )
  }

  if (!is.null(yi)) {
    names(res_main$es_mean_all) <- c(
      "es_mean_all_pe", "es_mean_all_lb", "es_mean_all_ub"
    )
    names(res_main$es_median_all) <- c(
      "es_median_all_pe", "es_median_all_lb", "es_median_all_ub"
    )
    names(res_main$es_mean_sig) <- c(
      "es_mean_sig_pe", "es_mean_sig_lb", "es_mean_sig_ub"
    )
    names(res_main$es_q25_sig) <- c(
      "es_q25_sig_pe", "es_q25_sig_lb", "es_q25_sig_ub"
    )
    names(res_main$es_median_sig) <- c(
      "es_median_sig_pe", "es_median_sig_lb", "es_median_sig_ub"
    )
    names(res_main$es_q75_sig) <- c(
      "es_q75_sig_pe", "es_q75_sig_lb", "es_q75_sig_ub"
    )
    names(res_main$es_tau_sig) <- c(
      "es_tau_sig_pe", "es_tau_sig_lb", "es_tau_sig_ub"
    )
    names(res_main$es_tau) <- c("es_tau_pe", "es_tau_lb", "es_tau_ub")

    if (is.matrix(res_main$local_es)) {
      if (length(local_interval_labels) == ncol(res_main$local_es)) {
        colnames(res_main$local_es) <- local_interval_labels
      } else if (is.null(colnames(res_main$local_es))) {
        colnames(res_main$local_es) <- paste0("local_es_", seq_len(ncol(res_main$local_es)))
      }
      rownames(res_main$local_es) <- c("pe", "lb", "ub")
    }

    if (is.matrix(res_main$local_tau)) {
      if (length(local_interval_labels) == ncol(res_main$local_tau)) {
        colnames(res_main$local_tau) <- local_interval_labels
      } else if (is.null(colnames(res_main$local_tau))) {
        colnames(res_main$local_tau) <- paste0("local_tau_", seq_len(ncol(res_main$local_tau)))
      }
      rownames(res_main$local_tau) <- c("pe", "lb", "ub")
    }

    results$effect_size <- list(
      mean_all = res_main$es_mean_all,
      median_all = res_main$es_median_all,
      mean_significant = res_main$es_mean_sig,
      q25_significant = res_main$es_q25_sig,
      median_significant = res_main$es_median_sig,
      q75_significant = res_main$es_q75_sig,
      heterogeneity_significant = res_main$es_tau_sig,
      heterogeneity = res_main$es_tau,
      local = res_main$local_es,
      local_heterogeneity = res_main$local_tau
    )
  }

  class(results) <- "zcurve"

  if (isTRUE(show_plot)) {
    if (!exists("plot.zcurve", mode = "function")) {
      stop("plot.zcurve() is unavailable. Load the package plotting file before calling zcurve().")
    }
    plot(results, xlim = c(x_lim_min, x_lim_max))
  }

  if (isTRUE(show_summary)) {
    print(summary(results))
    return(invisible(results))
  }

  results

} ### EOF zcurve


zcurve_interval_values <- function(x) {
  values <- as.numeric(x)
  c(
    estimate = if (length(values) >= 1) values[1] else NA_real_,
    ci_lower = if (length(values) >= 2) values[2] else NA_real_,
    ci_upper = if (length(values) >= 3) values[3] else NA_real_
  )
}


zcurve_estimate_table <- function(object) {
  estimates <- list(
    ODR = zcurve_interval_values(object$ODR),
    EDR = zcurve_interval_values(object$EDR),
    ERR = zcurve_interval_values(object$ERR),
    FDR = zcurve_interval_values(object$FDR),
    power_h1_significant = zcurve_interval_values(object$pow_h1_sig),
    ncp_tau = zcurve_interval_values(object$ncp_tau)
  )

  result <- data.frame(
    measure = c(
      "ODR", "EDR", "ERR", "FDR", "Power H1 | significant",
      "NCP heterogeneity (tau)"
    ),
    estimate = vapply(estimates, `[[`, numeric(1), "estimate"),
    row.names = NULL,
    check.names = FALSE
  )

  ci_lower <- vapply(estimates, `[[`, numeric(1), "ci_lower")
  ci_upper <- vapply(estimates, `[[`, numeric(1), "ci_upper")
  if (any(is.finite(ci_lower) & is.finite(ci_upper))) {
    result$ci_lower <- ci_lower
    result$ci_upper <- ci_upper
  }

  result
}


zcurve_effect_size_table <- function(effect_size) {
  estimates <- list(
    mean_significant = zcurve_interval_values(effect_size$mean_significant),
    median_significant = zcurve_interval_values(effect_size$median_significant),
    q25_significant = zcurve_interval_values(effect_size$q25_significant),
    q75_significant = zcurve_interval_values(effect_size$q75_significant),
    heterogeneity_significant = zcurve_interval_values(effect_size$heterogeneity_significant)
  )

  result <- data.frame(
    measure = c(
      "Mean", "Median", "25th percentile", "75th percentile",
      "Heterogeneity (tau)"
    ),
    estimate = vapply(estimates, `[[`, numeric(1), "estimate"),
    row.names = NULL,
    check.names = FALSE
  )

  ci_lower <- vapply(estimates, `[[`, numeric(1), "ci_lower")
  ci_upper <- vapply(estimates, `[[`, numeric(1), "ci_upper")
  if (any(is.finite(ci_lower) & is.finite(ci_upper))) {
    result$ci_lower <- ci_lower
    result$ci_upper <- ci_upper
  }

  result
}


zcurve_shape_diagnostic_table <- function(object) {
  diagnostics <- list(
    median = zcurve_interval_values(object$shape_d_median),
    mean = zcurve_interval_values(object$shape_d_mean)
  )

  result <- data.frame(
    statistic = c("Median difference", "Mean difference"),
    estimate = vapply(diagnostics, `[[`, numeric(1), "estimate"),
    row.names = NULL,
    check.names = FALSE
  )

  if (!any(is.finite(result$estimate))) {
    return(NULL)
  }

  ci_lower <- vapply(diagnostics, `[[`, numeric(1), "ci_lower")
  ci_upper <- vapply(diagnostics, `[[`, numeric(1), "ci_upper")
  if (any(is.finite(ci_lower) & is.finite(ci_upper))) {
    result$ci_lower <- ci_lower
    result$ci_upper <- ci_upper
  }

  result
}


zcurve_point_estimates <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(as.numeric(x[1, , drop = TRUE]))
  }
  as.numeric(x)
}


zcurve_print_header <- function(metadata) {
  version_text <- if (!is.null(metadata$version) && nzchar(metadata$version)) {
    paste0(" ", metadata$version)
  } else {
    ""
  }

  cat("zcurve", version_text, " fit\n", sep = "")
  cat("Method: ", metadata$method,
      " | Input: ", metadata$input,
      " | Tests: ", metadata$n_tests,
      " | Significant: ", metadata$n_significant,
      "\n", sep = "")

  if (!is.null(metadata$n_clusters)) {
    cat("Clusters: ", metadata$n_clusters, "\n", sep = "")
  }

  cat("Alpha: ", format(metadata$alpha),
      " | Bootstrap iterations: ", metadata$bootstrap_iterations,
      "\n", sep = "")
}


print.zcurve <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  zcurve_print_header(x$metadata)
  cat("\nEstimates:\n")
  print(zcurve_estimate_table(x), digits = digits, row.names = FALSE)

  if (!is.null(x$effect_size)) {
    cat("\nEffect sizes among significant studies:\n")
    print(zcurve_effect_size_table(x$effect_size), digits = digits, row.names = FALSE)
  }

  invisible(x)
}


summary.zcurve <- function(object, ...) {
  component_count <- min(
    length(zcurve_point_estimates(object$ncp)),
    length(zcurve_point_estimates(object$z_sd)),
    length(zcurve_point_estimates(object$w_inp)),
    length(zcurve_point_estimates(object$w_all))
  )

  components <- data.frame(
    component = seq_len(component_count),
    ncp = zcurve_point_estimates(object$ncp)[seq_len(component_count)],
    z_sd = zcurve_point_estimates(object$z_sd)[seq_len(component_count)],
    weight_selected = zcurve_point_estimates(object$w_inp)[seq_len(component_count)],
    weight_all = zcurve_point_estimates(object$w_all)[seq_len(component_count)],
    row.names = NULL
  )

  local_power_matrix <- as.matrix(object$local_power)
  local_power <- data.frame(
    interval = colnames(local_power_matrix),
    power = as.numeric(local_power_matrix[1, , drop = TRUE]),
    row.names = NULL
  )
  if (nrow(local_power_matrix) >= 3L &&
      any(is.finite(local_power_matrix[2, ]) & is.finite(local_power_matrix[3, ]))) {
    local_power$ci_lower <- as.numeric(local_power_matrix[2, ])
    local_power$ci_upper <- as.numeric(local_power_matrix[3, ])
  }

  local_effect_size <- NULL
  if (!is.null(object$effect_size$local) && ncol(object$effect_size$local) > 0) {
    local_es_matrix <- as.matrix(object$effect_size$local)
    local_effect_size <- data.frame(
      interval = colnames(local_es_matrix),
      effect_size = as.numeric(local_es_matrix[1, , drop = TRUE]),
      row.names = NULL
    )
    if (nrow(local_es_matrix) >= 3L &&
        any(is.finite(local_es_matrix[2, ]) & is.finite(local_es_matrix[3, ]))) {
      local_effect_size$ci_lower <- as.numeric(local_es_matrix[2, ])
      local_effect_size$ci_upper <- as.numeric(local_es_matrix[3, ])
    }
  }

  local_heterogeneity <- NULL
  if (!is.null(object$effect_size$local_heterogeneity) &&
      ncol(object$effect_size$local_heterogeneity) > 0) {
    local_tau_matrix <- as.matrix(object$effect_size$local_heterogeneity)
    local_heterogeneity <- data.frame(
      interval = colnames(local_tau_matrix),
      tau = as.numeric(local_tau_matrix[1, , drop = TRUE]),
      row.names = NULL
    )
    if (nrow(local_tau_matrix) >= 3L &&
        any(is.finite(local_tau_matrix[2, ]) & is.finite(local_tau_matrix[3, ]))) {
      local_heterogeneity$ci_lower <- as.numeric(local_tau_matrix[2, ])
      local_heterogeneity$ci_upper <- as.numeric(local_tau_matrix[3, ])
    }
  }

  result <- list(
    call = object$call,
    metadata = object$metadata,
    estimates = zcurve_estimate_table(object),
    components = components,
    local_power = local_power,
    diagnostics = list(
      slope = object$slope,
      extreme = object$extreme,
      shape = zcurve_shape_diagnostic_table(object),
      fit = object$fit,
      elapsed_seconds = object$metadata$elapsed_seconds
    ),
    effect_size = if (is.null(object$effect_size)) NULL else zcurve_effect_size_table(object$effect_size),
    local_effect_size = local_effect_size,
    local_heterogeneity = local_heterogeneity
  )

  class(result) <- "summary_zcurve"
  result
}


print.summary_zcurve <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  zcurve_print_header(x$metadata)

  cat("\nEstimates:\n")
  print(x$estimates, digits = digits, row.names = FALSE)

  if (!is.null(x$effect_size)) {
    cat("\nEffect sizes among significant studies:\n")
    print(x$effect_size, digits = digits, row.names = FALSE)
  }

  cat("\nMixture components:\n")
  print(x$components, digits = digits, row.names = FALSE)

  cat("\nLocal power:\n")
  print(x$local_power, digits = digits, row.names = FALSE)

  if (!is.null(x$local_effect_size)) {
    cat("\nLocal effect sizes:\n")
    print(x$local_effect_size, digits = digits, row.names = FALSE)
  }

  if (!is.null(x$local_heterogeneity)) {
    cat("\nLocal effect-size heterogeneity:\n")
    print(x$local_heterogeneity, digits = digits, row.names = FALSE)
  }

  if (!is.null(x$diagnostics$shape)) {
    cat("\nShape diagnostics:\n")
    print(x$diagnostics$shape, digits = digits, row.names = FALSE)
  }

  cat("\nElapsed time: ", format(round(x$diagnostics$elapsed_seconds, 3)), " seconds\n", sep = "")
  invisible(x)
}
