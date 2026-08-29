### Analysis of the cluster-bootstrap validation.
### Reads Cluster_Bootstrap_Validation_results.csv and reports, per cell and
### estimand, the coverage and mean CI width of the naive vs cluster bootstrap,
### with the empirical SD of the point estimate as the gold-standard yardstick.

options(scipen = 999)
project_dir <- "C:/Users/ulric/Documents/zcurve3"
res <- read.csv(file.path(project_dir, "Cluster_Bootstrap_Validation_results.csv"),
                stringsAsFactors = FALSE)

est_labels <- c("EDR", "ERR", "es_mean_sig", "es_tau_sig")
z975 <- qnorm(.975)

cover <- function(lb, ub, tr) mean(lb <= tr & tr <= ub, na.rm = TRUE)
width <- function(lb, ub)     mean(ub - lb, na.rm = TRUE)

cat("Replications per cell:\n"); print(table(res$cell))
cat("\nClusters per dataset (should be ~33, range 20-100):\n")
cl_summ <- do.call(rbind, tapply(res$n_clusters, res$cell,
                                 function(x) c(mean = mean(x), min = min(x), max = max(x))))
print(round(cl_summ, 1))

for (ci in sort(unique(res$cell))) {
  D <- res[res$cell == ci, ]
  cat(sprintf("\n================ cell %d  (d=%.1f, SD=%.1f, p_h0=%.2f, n=%d) ================\n",
              ci, D$d[1], D$SD[1], D$p_h0[1], nrow(D)))
  tab <- data.frame()
  for (e in est_labels) {
    tr <- D[[paste0("true_", e)]]
    n_lb <- D[[paste0("naive_",   e, "_lb")]]; n_ub <- D[[paste0("naive_",   e, "_ub")]]
    c_lb <- D[[paste0("cluster_", e, "_lb")]]; c_ub <- D[[paste0("cluster_", e, "_ub")]]
    pe   <- D[[paste0("cluster_", e, "_pe")]]          # point est identical to naive
    tab <- rbind(tab, data.frame(
      estimand   = e,
      true       = round(tr[1], 3),
      cov_naive  = round(cover(n_lb, n_ub, tr), 3),
      cov_clustr = round(cover(c_lb, c_ub, tr), 3),
      w_naive    = round(width(n_lb, n_ub), 3),
      w_clustr   = round(width(c_lb, c_ub), 3),
      emp_SD     = round(sd(pe, na.rm = TRUE), 3),
      SE_naive   = round(width(n_lb, n_ub) / (2 * z975), 3),
      SE_clustr  = round(width(c_lb, c_ub) / (2 * z975), 3)
    ))
  }
  print(tab, row.names = FALSE)
}

cat("\nReading: cov_* is CI coverage of the fixed truth (target .95).\n")
cat("SE_* = mean CI width / 3.92; compare to emp_SD (the true sampling SD).\n")
cat("Cluster bootstrap validated if cov_clustr ~ .95 >= cov_naive and SE_clustr ~ emp_SD > SE_naive.\n")
