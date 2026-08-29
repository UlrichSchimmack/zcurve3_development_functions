

### Load weightr in main program. 

if (!"parallel" %in% loadedNamespaces()) {
  library(parallel)
}


if (!"weightr" %in% loadedNamespaces()) {
  library(weightr)
}


print.myresult <- function(x, ...) {
  cat(as.character(x), sep = "")
  invisible(x)
}

### Clustered Step-Function Selection Model (Weight-R with parallel bootstrap)

#nboot = 500;steps=c(.5,.025);ncores=16;seed = 2026
#yi = tm$yi; vi = tm$vi;cluster.id = tm$cluster.id


run_boot_cluster_weightr <- function(
  yi, vi, 
  cluster.id,
  nboot = 500,
  steps = c(.5, .05, .025),
  weights = NULL,
  fe = FALSE,
  ncores = parallel::detectCores() - 1,
  seed = 2026
) {

  print("Version 26.06.20")

  # ---------------------------------------------------------
  # Original estimate
  # ---------------------------------------------------------


yy <- yi[sel]
vv <- vi[sel]
cc <- cluster_id[sel]

sort(vv[is.finite(vv) & vv > 0])[1:10]

  if (is.null(weights)) {
    orig <- weightr::weightfunct(
      yi, vi,
      steps = steps,
      fe = fe
    )
  } else {
    orig <- weightfunct(
      yi, vi,
      steps = steps,
      weights = weights,
      fe = fe
    )
  }

  adj_par   <- orig[[2]]$par
  unadj_par <- orig[[1]]$par

  # ---------------------------------------------------------
  # Extract parameters: FE and RE have different structures
  # ---------------------------------------------------------

  if (fe) {

    # FE: mean + estimated weights
    est_mean_orig <- adj_par[1]
    tau_orig <- 0

    if (is.null(weights)) {
      weights_orig <- adj_par[-1]
      n_weights <- length(weights_orig)
    } else {
      weights_orig <- numeric(0)
      n_weights <- 0
    }

    boot_names <- c(
      "adj_mean",
      "unadj_mean",
      paste0("weight_", seq_len(n_weights))
    )

  } else {

    # RE: tau2 + mean + estimated weights
    tau_orig      <- sqrt(adj_par[1])
    est_mean_orig <- adj_par[2]

    if (is.null(weights)) {
      weights_orig <- adj_par[-c(1, 2)]
      n_weights <- length(weights_orig)
    } else {
      weights_orig <- numeric(0)
      n_weights <- 0
    }

    boot_names <- c(
      "adj_tau2",
      "unadj_tau2",
      "adj_mean",
      paste0("weight_", seq_len(n_weights))
    )
  }


  # ---------------------------------------------------------
  # Interval labels
  # ---------------------------------------------------------

  steps_full <- sort(unique(c(steps, 1)))
  boundaries <- c(0, steps_full)

  interval_labels <- paste0(
    "p in (",
    boundaries[-length(boundaries)],
    ", ",
    boundaries[-1],
    "]"
  )


  # ---------------------------------------------------------
  # Cluster bootstrap
  # ---------------------------------------------------------

  clusters <- unique(cluster.id)
  k_clust  <- length(clusters)

  cat(
    "Clusters:", k_clust,
    " | Effects:", length(yi),
    " | Bootstrap iterations:", nboot, "\n"
  )

  set.seed(seed)

  boot_samples <- lapply(
    seq_len(nboot),
    function(i) sample(clusters, k_clust, replace = TRUE)
  )


  boot_one <- function(b) {

    idx <- unlist(
      lapply(
        boot_samples[[b]],
        function(cid) which(cluster.id == cid)
      ),
      use.names = FALSE
    )

    yi_b <- yi[idx]
    vi_b <- vi[idx]

    res <- tryCatch({

      if (is.null(weights)) {

        fit <- weightfunct(
          yi_b, vi_b,
          steps = steps,
          fe = fe
        )

      } else {

        fit <- weightfunct(
          yi_b, vi_b,
          steps = steps,
          weights = weights,
          fe = fe
        )
      }

      if (fe) {

        # FE: adjusted mean, unadjusted mean, weights
        if (is.null(weights)) {

          c(
            adj_mean   = fit[[2]]$par[1],
            unadj_mean = fit[[1]]$par[1],
            fit[[2]]$par[-1]
          )

        } else {

          c(
            adj_mean   = fit[[2]]$par[1],
            unadj_mean = fit[[1]]$par[1]
          )
        }

      } else {

        # RE: adjusted tau2, unadjusted tau2, mean, weights
        if (is.null(weights)) {

          c(
            adj_tau2   = fit[[2]]$par[1],
            unadj_tau2 = fit[[1]]$par[1],
            adj_mean   = fit[[2]]$par[2],
            fit[[2]]$par[-c(1, 2)]
          )

        } else {

          c(
            adj_tau2   = fit[[2]]$par[1],
            unadj_tau2 = fit[[1]]$par[1],
            adj_mean   = fit[[2]]$par[2]
          )
        }
      }

    }, error = function(e) {

      rep(NA_real_, length(boot_names))

    })

    return(res)
  }


  # ---------------------------------------------------------
  # Parallel execution
  # ---------------------------------------------------------

  ncores <- max(1, ncores)

  cl <- parallel::makeCluster(ncores)

  clusterEvalQ(cl, library(weightr))

  clusterExport(
    cl,
    varlist = c(
      "yi", "vi", "cluster.id",
      "boot_samples",
      "steps", "weights", "fe",
      "boot_names"
    ),
    envir = environment()
  )

  cat(
    "Running", nboot,
    "bootstrap iterations on",
    ncores, "cores...\n"
  )

  t0 <- proc.time()

  boot_results <- parLapply(
    cl,
    seq_len(nboot),
    boot_one
  )

  stopCluster(cl)

  elapsed <- (proc.time() - t0)[3]

  cat(
    "Done in",
    round(elapsed, 1),
    "seconds\n\n"
  )


  # ---------------------------------------------------------
  # Combine results
  # ---------------------------------------------------------

  boot_mat <- do.call(rbind, boot_results)

  # Safety check -- useful if something ever changes in weightr
  if (ncol(boot_mat) != length(boot_names)) {

    stop(
      paste0(
        "Bootstrap returned ",
        ncol(boot_mat),
        " columns but ",
        length(boot_names),
        " column names were expected. ",
        "Lengths returned: ",
        paste(unique(lengths(boot_results)), collapse = ", ")
      )
    )
  }

  colnames(boot_mat) <- boot_names


  # ---------------------------------------------------------
  # Remove failed fits
  # ---------------------------------------------------------

  n_fail <- sum(!complete.cases(boot_mat))

  cat(
    "Convergence failures:",
    n_fail, "of", nboot,
    "(",
    round(100 * n_fail / nboot, 1),
    "%)\n"
  )

  boot_clean <- boot_mat[
    complete.cases(boot_mat),
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------
  # Mean CI
  # ---------------------------------------------------------

  ci_mean_adj <- quantile(
    boot_clean[, "adj_mean"],
    c(.025, .975),
    na.rm = TRUE
  )

  med_adj <- median(
    boot_clean[, "adj_mean"],
    na.rm = TRUE
  )


  # ---------------------------------------------------------
  # FE versus RE output
  # ---------------------------------------------------------

  if (fe) {

    output_mu_tau <- paste0(
      "\n--- Results ---\n",
      sprintf(
        "Adjusted mean:     %.4f    %.4f   [%.4f, %.4f]\n",
        est_mean_orig,
        med_adj,
        ci_mean_adj[1],
        ci_mean_adj[2]
      )
    )

    pred.interval <- c(NA, NA)
    output_pi <- ""

  } else {

    ci_tau_adj <- quantile(
      sqrt(boot_clean[, "adj_tau2"]),
      c(.025, .975),
      na.rm = TRUE
    )

    med_tau <- median(
      sqrt(boot_clean[, "adj_tau2"]),
      na.rm = TRUE
    )

    pred_draws <- unlist(
      lapply(seq_len(nrow(boot_clean)), function(bb) {

        rnorm(
          500,
          mean = boot_clean[bb, "adj_mean"],
          sd = sqrt(max(boot_clean[bb, "adj_tau2"], 0))    
        )
      })
    )

    pred.interval <- quantile(
      pred_draws,
      c(.025, .975),
      na.rm = TRUE
    )

    output_mu_tau <- paste0(
      "\n--- Results ---\n",
      sprintf(
        "Adjusted mean:     %.4f    %.4f   [%.4f, %.4f]\n",
        est_mean_orig,
        med_adj,
        ci_mean_adj[1],
        ci_mean_adj[2]
      ),
      sprintf(
        "Tau:               %.4f    %.4f   [%.4f, %.4f]\n",
        tau_orig,
        med_tau,
        ci_tau_adj[1],
        ci_tau_adj[2]
      )
    )

    output_pi <- sprintf(
      "Prediction Interval ranges from %.2f to %.2f\n",
      pred.interval[1],
      pred.interval[2]
    )
  }


  # ---------------------------------------------------------
  # Selection weights
  # ---------------------------------------------------------

  output_weight <- ""

  if (n_weights > 0) {

    for (j in seq_len(n_weights)) {

      wname <- paste0("weight_", j)

      ci_w <- quantile(
        boot_clean[, wname],
        c(.025, .975),
        na.rm = TRUE
      )

      med_w <- median(
        boot_clean[, wname],
        na.rm = TRUE
      )

      output_weight[j] <- sprintf(
        "Weight %d (%s):\n                   %.4f    %.4f   [%.4f, %.4f]\n",
        j,
        interval_labels[j + 1],
        weights_orig[j],
        med_w,
        ci_w[1],
        ci_w[2]
      )
    }

    output_weight <- paste0(
      output_weight,
      collapse = ""
    )
  }


  # ---------------------------------------------------------
  # Final output
  # ---------------------------------------------------------

  output <- paste0(
    output_mu_tau,
    output_weight,
    output_pi
  )

  cluster_results <- structure(
    output,
    class = c("myresult", "character")
  )

  res <- list(
    original_results = orig,
    cluster_results  = cluster_results,
    boot_matrix      = boot_mat,
    n_fail           = n_fail,
    pred.interval    = pred.interval
  )

  return(res)
}

##############################################

##############################################

##############################################


run_boot_cluster_petpeese <- function(
  yi,
  vi,
  cluster.id,
  nboot = 500,
  ncores = parallel::detectCores() - 1,
  seed = 2026,
  alpha = .05
) {

  # =========================================================
  # 1. CHECK AND CLEAN DATA
  # =========================================================

  if (length(yi) != length(vi) ||
      length(yi) != length(cluster.id)) {
    stop("yi, vi, and cluster.id must have the same length.")
  }

  ok <- is.finite(yi) &
        is.finite(vi) &
        vi > 0 &
        !is.na(cluster.id)

  yi         <- yi[ok]
  vi         <- vi[ok]
  cluster.id <- cluster.id[ok]

  if (length(yi) < 3) {
    stop("Too few valid effect sizes.")
  }

  clusters <- unique(cluster.id)
  k_clust  <- length(clusters)

  if (k_clust < 3) {
    stop("Too few clusters for cluster-robust inference.")
  }

  ncores <- max(1, ncores)
  ncores <- min(ncores, nboot)


  # =========================================================
  # 2. FIT PET AND PEESE
  #
  # PET:
  #   yi = beta0 + beta1 * SE
  #
  # PEESE:
  #   yi = beta0 + beta1 * SE^2
  #
  # In both cases:
  #   weights = 1 / vi
  #
  # Cluster-robust inference is used for coefficient SEs,
  # p-values, and the PET decision rule.
  # =========================================================

  fit_one <- function(yi, vi, cluster.id, alpha = .05) {

    sei <- sqrt(vi)

    # -------------------------------------------------------
    # PET: effect size ~ SE
    # -------------------------------------------------------

    pet_fit <- metafor::rma(
      yi = yi,
      vi = vi,
      mods = ~ sei,
      method = "FE"
    )

    pet_rob <- metafor::robust(
      pet_fit,
      cluster = cluster.id,
      adjust = TRUE,
      clubSandwich = FALSE
    )


    # -------------------------------------------------------
    # PEESE: effect size ~ variance
    # -------------------------------------------------------

    peese_fit <- metafor::rma(
      yi = yi,
      vi = vi,
      mods = ~ vi,
      method = "FE"
    )

    peese_rob <- metafor::robust(
      peese_fit,
      cluster = cluster.id,
      adjust = TRUE,
      clubSandwich = FALSE
    )


    # -------------------------------------------------------
    # PET intercept
    # -------------------------------------------------------

    pet_est <- as.numeric(coef(pet_fit)[1])
    pet_se  <- as.numeric(pet_rob$se[1])
    pet_p   <- as.numeric(pet_rob$pval[1])

    pet_lb <- as.numeric(pet_rob$ci.lb[1])
    pet_ub <- as.numeric(pet_rob$ci.ub[1])


    # -------------------------------------------------------
    # PET slope = funnel asymmetry / small-study association
    # -------------------------------------------------------

    pet_slope    <- as.numeric(coef(pet_fit)[2])
    pet_slope_se <- as.numeric(pet_rob$se[2])
    pet_slope_p  <- as.numeric(pet_rob$pval[2])

    pet_slope_lb <- as.numeric(pet_rob$ci.lb[2])
    pet_slope_ub <- as.numeric(pet_rob$ci.ub[2])


    # -------------------------------------------------------
    # PEESE intercept
    # -------------------------------------------------------

    peese_est <- as.numeric(coef(peese_fit)[1])
    peese_se  <- as.numeric(peese_rob$se[1])
    peese_p   <- as.numeric(peese_rob$pval[1])

    peese_lb <- as.numeric(peese_rob$ci.lb[1])
    peese_ub <- as.numeric(peese_rob$ci.ub[1])


    # -------------------------------------------------------
    # PEESE slope
    # -------------------------------------------------------

    peese_slope    <- as.numeric(coef(peese_fit)[2])
    peese_slope_se <- as.numeric(peese_rob$se[2])
    peese_slope_p  <- as.numeric(peese_rob$pval[2])

    peese_slope_lb <- as.numeric(peese_rob$ci.lb[2])
    peese_slope_ub <- as.numeric(peese_rob$ci.ub[2])


    # -------------------------------------------------------
    # Conditional PET-PEESE rule
    #
    # PET significant -> use PEESE
    # PET non-significant -> use PET
    # -------------------------------------------------------

    if (is.finite(pet_p) && pet_p < alpha) {

      pp_est <- peese_est
      method_used <- 1    # PEESE

    } else {

      pp_est <- pet_est
      method_used <- 0    # PET
    }


    # -------------------------------------------------------
    # Return all relevant quantities
    # -------------------------------------------------------

    c(
      PET_ES       = pet_est,
      PET_SE       = pet_se,
      PET_P        = pet_p,
      PET_LB       = pet_lb,
      PET_UB       = pet_ub,

      PET_SLOPE    = pet_slope,
      PET_SLOPE_SE = pet_slope_se,
      PET_SLOPE_P  = pet_slope_p,
      PET_SLOPE_LB = pet_slope_lb,
      PET_SLOPE_UB = pet_slope_ub,

      PEESE_ES       = peese_est,
      PEESE_SE       = peese_se,
      PEESE_P        = peese_p,
      PEESE_LB       = peese_lb,
      PEESE_UB       = peese_ub,

      PEESE_SLOPE    = peese_slope,
      PEESE_SLOPE_SE = peese_slope_se,
      PEESE_SLOPE_P  = peese_slope_p,
      PEESE_SLOPE_LB = peese_slope_lb,
      PEESE_SLOPE_UB = peese_slope_ub,

      PP_ES      = pp_est,
      PP_METHOD  = method_used
    )
  }


  # =========================================================
  # 3. ORIGINAL FIT
  # =========================================================

  orig <- fit_one(
    yi = yi,
    vi = vi,
    cluster.id = cluster.id,
    alpha = alpha
  )


  # =========================================================
  # 4. GENERATE CLUSTER BOOTSTRAP SAMPLES
  # =========================================================

  set.seed(seed)

  boot_samples <- lapply(
    seq_len(nboot),
    function(b) {
      sample(
        clusters,
        size = k_clust,
        replace = TRUE
      )
    }
  )

  cat(
    "Clusters:", k_clust,
    "| Effects:", length(yi),
    "| Bootstrap iterations:", nboot,
    "\n"
  )


  # =========================================================
  # 5. ONE CLUSTER BOOTSTRAP REPLICATION
  # =========================================================

  boot_one <- function(b) {

    sampled_clusters <- boot_samples[[b]]

    # Each sampled occurrence becomes a separate bootstrap
    # cluster. This matters when the same original cluster
    # is sampled more than once.

    idx_list <- lapply(
      sampled_clusters,
      function(cid) {
        which(cluster.id == cid)
      }
    )

    idx <- unlist(
      idx_list,
      use.names = FALSE
    )

    cluster_b <- unlist(
      lapply(
        seq_along(idx_list),
        function(j) {
          rep(j, length(idx_list[[j]]))
        }
      ),
      use.names = FALSE
    )

    tryCatch(

      fit_one(
        yi = yi[idx],
        vi = vi[idx],
        cluster.id = cluster_b,
        alpha = alpha
      ),

      error = function(e) {

        out <- rep(
          NA_real_,
          length(orig)
        )

        names(out) <- names(orig)

        out
      }
    )
  }


  # =========================================================
  # 6. PARALLEL BOOTSTRAP
  # =========================================================

  cl <- parallel::makeCluster(ncores)

  on.exit(
    {
      if (!is.null(cl)) {
        try(
          parallel::stopCluster(cl),
          silent = TRUE
        )
      }
    },
    add = TRUE
  )

  parallel::clusterExport(
    cl,
    varlist = c(
      "yi",
      "vi",
      "cluster.id",
      "boot_samples",
      "fit_one",
      "orig",
      "alpha"
    ),
    envir = environment()
  )

  cat(
    "Running",
    nboot,
    "bootstrap iterations on",
    ncores,
    "cores...\n"
  )

  t0 <- proc.time()

  boot_results <- parallel::parLapply(
    cl,
    seq_len(nboot),
    boot_one
  )

  parallel::stopCluster(cl)
  cl <- NULL

  elapsed <- (proc.time() - t0)[3]

  cat(
    "Done in",
    round(elapsed, 1),
    "seconds\n"
  )


  # =========================================================
  # 7. COMBINE BOOTSTRAP RESULTS
  # =========================================================

  boot_mat <- do.call(
    rbind,
    boot_results
  )

  if (ncol(boot_mat) != length(orig)) {

    stop(
      "Unexpected number of columns in bootstrap output."
    )
  }

  colnames(boot_mat) <- names(orig)

  n_fail <- sum(
    !complete.cases(boot_mat)
  )

  cat(
    "Bootstrap failures:",
    n_fail,
    "of",
    nboot,
    "(",
    round(100 * n_fail / nboot, 1),
    "%)\n"
  )

  boot_clean <- boot_mat[
    complete.cases(boot_mat),
    ,
    drop = FALSE
  ]

  if (nrow(boot_clean) < 10) {
    stop("Too few successful bootstrap replications.")
  }


  # =========================================================
  # 8. BOOTSTRAP SUMMARY FUNCTION
  # =========================================================

  boot_summary <- function(var) {

    x <- boot_clean[, var]

    q <- stats::quantile(
      x,
      probs = c(.025, .50, .975),
      na.rm = TRUE,
      names = FALSE
    )

    c(
      median = q[2],
      lb = q[1],
      ub = q[3]
    )
  }


  # =========================================================
  # 9. BOOTSTRAP SUMMARIES
  # =========================================================

  PET_boot <- boot_summary("PET_ES")

  PET_slope_boot <- boot_summary(
    "PET_SLOPE"
  )

  PEESE_boot <- boot_summary(
    "PEESE_ES"
  )

  PEESE_slope_boot <- boot_summary(
    "PEESE_SLOPE"
  )

  PP_boot <- boot_summary(
    "PP_ES"
  )


  # ---------------------------------------------------------
  # Proportion of bootstrap samples selecting PEESE
  # ---------------------------------------------------------

  prop_boot_peese <- mean(
    boot_clean[, "PP_METHOD"] == 1,
    na.rm = TRUE
  )


  # =========================================================
  # 10. FINAL RETURN OBJECT
  # =========================================================

  list(

    k = length(yi),

    k_clusters = k_clust,


    # -------------------------------------------------------
    # Original complete result vector
    # -------------------------------------------------------

    original = orig,


    # -------------------------------------------------------
    # PET corrected effect
    # -------------------------------------------------------

    pet = c(
      pe = unname(orig["PET_ES"]),
      se = unname(orig["PET_SE"]),
      p  = unname(orig["PET_P"]),

      median = unname(PET_boot["median"]),
      lb = unname(PET_boot["lb"]),
      ub = unname(PET_boot["ub"])
    ),


    # -------------------------------------------------------
    # PET slope
    # -------------------------------------------------------

    pet_slope = c(
      pe = unname(orig["PET_SLOPE"]),
      se = unname(orig["PET_SLOPE_SE"]),
      p  = unname(orig["PET_SLOPE_P"]),

      model_lb = unname(orig["PET_SLOPE_LB"]),
      model_ub = unname(orig["PET_SLOPE_UB"]),

      median = unname(PET_slope_boot["median"]),
      boot_lb = unname(PET_slope_boot["lb"]),
      boot_ub = unname(PET_slope_boot["ub"])
    ),


    # -------------------------------------------------------
    # PEESE corrected effect
    # -------------------------------------------------------

    peese = c(
      pe = unname(orig["PEESE_ES"]),
      se = unname(orig["PEESE_SE"]),
      p  = unname(orig["PEESE_P"]),

      median = unname(PEESE_boot["median"]),
      lb = unname(PEESE_boot["lb"]),
      ub = unname(PEESE_boot["ub"])
    ),


    # -------------------------------------------------------
    # PEESE slope
    # -------------------------------------------------------

    peese_slope = c(
      pe = unname(orig["PEESE_SLOPE"]),
      se = unname(orig["PEESE_SLOPE_SE"]),
      p  = unname(orig["PEESE_SLOPE_P"]),

      model_lb = unname(orig["PEESE_SLOPE_LB"]),
      model_ub = unname(orig["PEESE_SLOPE_UB"]),

      median = unname(PEESE_slope_boot["median"]),
      boot_lb = unname(PEESE_slope_boot["lb"]),
      boot_ub = unname(PEESE_slope_boot["ub"])
    ),


    # -------------------------------------------------------
    # Conditional PET-PEESE estimate
    # -------------------------------------------------------

    petpeese = c(
      pe = unname(orig["PP_ES"]),

      median = unname(PP_boot["median"]),
      lb = unname(PP_boot["lb"]),
      ub = unname(PP_boot["ub"])
    ),


    # -------------------------------------------------------
    # Which estimator was selected in original sample?
    # -------------------------------------------------------

    method_original = ifelse(
      orig["PP_METHOD"] == 1,
      "PEESE",
      "PET"
    ),


    # -------------------------------------------------------
    # Stability of PET versus PEESE decision
    # -------------------------------------------------------

    prop_boot_peese = prop_boot_peese,


    # -------------------------------------------------------
    # Full bootstrap matrix
    # -------------------------------------------------------

    boot_matrix = boot_mat,

    n_fail = n_fail
  )
}


##############################################

##############################################

##############################################





extract_petpeese <- function(x) {

  c(
    k = x$k,
    k_clusters = x$k_clusters,

    # PET
    PET_ES = unname(x$pet["pe"]),
    PET_SE = unname(x$pet["se"]),
    PET_P  = unname(x$pet["p"]),
    PET_LB = unname(x$pet["lb"]),
    PET_UB = unname(x$pet["ub"]),

    # PET slope
    PET_SLOPE    = unname(x$pet_slope["pe"]),
    PET_SLOPE_SE = unname(x$pet_slope["se"]),
    PET_SLOPE_P  = unname(x$pet_slope["p"]),
    PET_SLOPE_LB = unname(x$pet_slope["boot_lb"]),
    PET_SLOPE_UB = unname(x$pet_slope["boot_ub"]),

    # PEESE
    PEESE_ES = unname(x$peese["pe"]),
    PEESE_SE = unname(x$peese["se"]),
    PEESE_P  = unname(x$peese["p"]),
    PEESE_LB = unname(x$peese["lb"]),
    PEESE_UB = unname(x$peese["ub"]),

    # PEESE slope
    PEESE_SLOPE    = unname(x$peese_slope["pe"]),
    PEESE_SLOPE_SE = unname(x$peese_slope["se"]),
    PEESE_SLOPE_P  = unname(x$peese_slope["p"]),
    PEESE_SLOPE_LB = unname(x$peese_slope["boot_lb"]),
    PEESE_SLOPE_UB = unname(x$peese_slope["boot_ub"]),

    # Conditional PET-PEESE
    PP_ES = unname(x$petpeese["pe"]),
    PP_LB = unname(x$petpeese["lb"]),
    PP_UB = unname(x$petpeese["ub"]),

    # 0 = PET, 1 = PEESE
    PP_METHOD = ifelse(
      x$method_original == "PEESE",
      1,
      0
    ),

    # proportion of bootstrap samples choosing PEESE
    PP_PEESE_PROP = x$prop_boot_peese,

    n_fail = x$n_fail
  )
}

##############################################

##############################################

##############################################



extract_weightr <- function(x) {

  r <- x$original_results

  ## -----------------------------
  ## Basic model information
  ## -----------------------------

  k     <- r$k
  p     <- r$p
  steps <- r$steps

  ## adjusted parameters
  ## par[1] = tau^2
  ## par[2] = adjusted mean
  ## par[3:] = selection weights
  par <- as.numeric(r$output_adj$par)

  tau2_adj <- par[1]
  tau_adj  <- sqrt(tau2_adj)
  mean_adj <- par[2]

  weights <- par[3:length(par)]

  ## unadjusted parameters
  par_unadj <- as.numeric(r$output_unadj$par)

  tau2_unadj <- par_unadj[1]
  tau_unadj  <- sqrt(tau2_unadj)
  mean_unadj <- par_unadj[2]


  ## -----------------------------
  ## Count observations in p bins
  ## intervals:
  ## (-Inf,.025], (.025,.5],
  ## (.5,.975], (.975,1]
  ## -----------------------------

  bins <- cut(
    p,
    breaks = c(-Inf, steps),
    right = TRUE,
    include.lowest = TRUE
  )

  bin_counts <- as.numeric(table(bins))

  names(bin_counts) <- c(
    "p_le_.025",
    "p_.025_.5",
    "p_.5_.975",
    "p_gt_.975"
  )


  ## -----------------------------
  ## Full weight vector
  ## reference significant bin = 1
  ## -----------------------------

  weights_full <- c(1, weights)

  names(weights_full) <- names(bin_counts)


  ## -----------------------------
  ## Estimated number before selection
  ## observed_j = latent_j * weight_j
  ## -----------------------------

  latent_counts <- bin_counts / weights_full


  ## Overall relative publication weight
  ## count-weighted harmonic mean
  w_overall <- sum(bin_counts) /
               sum(latent_counts)


  ## -----------------------------
  ## Two-sided ODR and implied EDR
  ##
  ## p <= .025 and p > .975
  ## correspond to two-sided p < .05
  ## -----------------------------

  ODR <- (bin_counts[1] + bin_counts[4]) /
         sum(bin_counts)

  EDR <- (latent_counts[1] + latent_counts[4]) /
         sum(latent_counts)

  bias <- ODR - EDR


  ## -----------------------------
  ## Bootstrap quantities
  ## -----------------------------

  boot <- x$boot_matrix

  mean_adj_boot <- mean(boot[, "adj_mean"], na.rm = TRUE)

  mean_adj_ci <- quantile(
    boot[, "adj_mean"],
    c(.025, .975),
    na.rm = TRUE
  )

  tau_adj_boot <- sqrt(
    mean(boot[, "adj_tau2"], na.rm = TRUE)
  )

  tau_adj_ci <- sqrt(
    quantile(
      boot[, "adj_tau2"],
      c(.025, .975),
      na.rm = TRUE
    )
  )


  ## return one named vector
  c(
    k = k,

    k_sig_pos = bin_counts[1],
    k_ns_pos  = bin_counts[2],
    k_ns_neg  = bin_counts[3],
    k_sig_neg = bin_counts[4],

    mean_unadj = mean_unadj,
    tau_unadj  = tau_unadj,

    mean_adj = mean_adj,
    tau_adj  = tau_adj,

    weight_1 = weights[1],
    weight_2 = weights[2],
    weight_3 = weights[3],

    weight_overall = w_overall,

    ODR_weightr = ODR,
    EDR_weightr = EDR,
    bias_weightr = bias,

    mean_adj_boot = mean_adj_boot,
    mean_adj_lb = mean_adj_ci[1],
    mean_adj_ub = mean_adj_ci[2],

    tau_adj_boot = tau_adj_boot,
    tau_adj_lb = tau_adj_ci[1],
    tau_adj_ub = tau_adj_ci[2],

    n_fail = x$n_fail
  )
}



