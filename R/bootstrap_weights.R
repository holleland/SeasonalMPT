# =====================================================================
#  bootstrap_weights.R
#
#  Quantifies the effect of covariance-estimation uncertainty on the
#  optimal portfolio weights and Seasonal Risk Score, via a moving-block
#  (circular) residual bootstrap, K = 4 fixed:
#    1. Fit the Fourier regression  X_t = mu + A c_t + Z_t.
#    2. Resample the seasonally-adjusted residual vectors Z_t in
#       contiguous blocks of length L (blocks preserve the short-lag
#       autocorrelation documented in Appendix F; the whole residual
#       vector is resampled so cross-asset dependence is kept).
#    3. Rebuild X*_t = mu_hat + A_hat c_t + Z*_t on the original time
#       grid, re-estimate (mu, A, Sigma_SEASON, Sigma_Z, Sigma), and
#       re-solve the three benchmark portfolios (psi = 0, 0.5, 1).
#    4. Collect the optimal solar share and SRS over B replications and
#       report percentile 95% confidence intervals + bootstrap SE.
#
#  The optimiser is the exact quadratic-programming solution used throughout
#  Section 4 (Section 4.1 was updated from the coarse target-grid search to
#  this exact solver; the difference is <= 1 percentage point):
#    Case 4.1 - minimise w' M(psi) w over {w'1 = 1, w >= 0};
#    Case 4.2 - solve.QP with the demand constraint of
#               norwegian_mix_of_solar_and_wind_with_demand.R at p_target = 100%.
#
#  Also reports sensitivity of the intervals to the block length L.
#
#  Run with working directory = repository root.
#  Outputs: output/bootstrap_*.csv  and  figures/bootstrap_weights.png
# =====================================================================

source("R/robustness_helpers.R")
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

set.seed(2024)
K       <- 4
B       <- 2000          # replications for the main results
L_main  <- 30            # block length for the main results (~1 month)
L_grid  <- c(15, 30, 60) # block lengths for the sensitivity check
B_Lsens <- 1000

strategies <- c("Season-adjusted" = 0, "Balanced" = 0.5, "Seasonal" = 1)

# ---------------------------------------------------------------------
# Circular block-bootstrap index vector of length n
# ---------------------------------------------------------------------
cbb_index <- function(n, L) {
  nb     <- ceiling(n / L)
  starts <- sample.int(n, nb, replace = TRUE)
  idx    <- as.vector(vapply(starts,
                             function(s) as.integer((s + 0:(L - 1) - 1) %% n) + 1L,
                             integer(L)))
  idx[seq_len(n)]
}

# ---------------------------------------------------------------------
# One bootstrap pass for a case.
#   opt_fun(SigS, SigZ, Sig, mu, psi) -> weight vector
#   summ_fun(w, SigS, SigZ, Sig, mu)  -> named numeric (quantities of interest)
# ---------------------------------------------------------------------
run_bootstrap <- function(Y, K, tau, B, L, opt_fun, summ_fun, label) {
  n <- nrow(Y); m <- ncol(Y)
  D <- cbind(1, fourier_design(n, K, tau))          # fixed design
  P <- solve(crossprod(D), t(D))                    # (2K+1) x n  hat operator
  co0    <- P %*% Y
  fit0   <- D %*% co0
  resid0 <- Y - fit0

  one_rep <- function(resmat) {
    Xs   <- fit0 + resmat
    co   <- P %*% Xs
    mu   <- co[1, ]
    A    <- t(co[-1, , drop = FALSE])
    SigS <- n / (n - 1) * A %*% t(A) / 2
    SigZ <- cov(Xs - D %*% co)
    Sig  <- cov(Xs)
    dimnames(SigS) <- dimnames(SigZ) <- dimnames(Sig) <- list(colnames(Y), colnames(Y))
    out <- lapply(names(strategies), function(nm) {
      w <- tryCatch(opt_fun(SigS, SigZ, Sig, mu, strategies[[nm]]),
                    error = function(e) rep(NA_real_, m))
      c(strategy = nm, summ_fun(w, SigS, SigZ, Sig, mu))
    })
    do.call(rbind, out)
  }

  # point estimate (observed residuals, no resampling)
  point <- as.data.frame(one_rep(resid0), stringsAsFactors = FALSE)

  reps <- vector("list", B)
  for (b in seq_len(B)) reps[[b]] <- one_rep(resid0[cbb_index(n, L), , drop = FALSE])
  boot <- as.data.frame(do.call(rbind, reps), stringsAsFactors = FALSE)

  num <- setdiff(names(boot), "strategy")
  boot[num]  <- lapply(boot[num],  as.numeric)
  point[num] <- lapply(point[num], as.numeric)

  summ <- boot %>%
    pivot_longer(all_of(num), names_to = "quantity", values_to = "value") %>%
    group_by(strategy, quantity) %>%
    summarise(boot_mean = mean(value, na.rm = TRUE),
              boot_se   = sd(value,  na.rm = TRUE),
              lwr       = quantile(value, 0.025, na.rm = TRUE),
              med       = quantile(value, 0.5,   na.rm = TRUE),
              upr       = quantile(value, 0.975, na.rm = TRUE),
              n_ok      = sum(!is.na(value)),
              .groups = "drop") %>%
    left_join(point %>%
                pivot_longer(all_of(num), names_to = "quantity", values_to = "estimate"),
              by = c("strategy", "quantity")) %>%
    mutate(case = label, L = L, B = B, .before = 1) %>%
    mutate(strategy = factor(strategy, levels = names(strategies))) %>%
    arrange(quantity, strategy)

  list(summary = summ, draws = boot %>% mutate(case = label, L = L))
}

# =====================================================================
# CASE 4.1
# =====================================================================
cat("\n=== Case 4.1 bootstrap (B =", B, ", L =", L_main, ") ===\n")
d1 <- load_case1_data()

# Case 4.1 optimiser: exact QP minimiser of w' M(psi) w over the long-only
# simplex, M(psi) = psi * Sigma_SEASON + (1 - psi) * Sigma_Z.
opt1 <- function(SigS, SigZ, Sig, mu, psi)
  opt_case1(risk_matrix(SigS, SigZ, psi))
summ1 <- function(w, SigS, SigZ, Sig, mu) {
  if (anyNA(w)) return(c(solar_weight = NA, SRS = NA, port_var = NA, rho_season = NA))
  c(solar_weight = unname(w[1]),
    SRS          = srs(w, SigS, SigZ),
    port_var     = port_var(w, Sig),
    rho_season   = seasonal_cor(SigS))
}

bt1 <- run_bootstrap(d1$Y, K, d1$tau, B, L_main, opt1, summ1, "Case 4.1")
print(bt1$summary %>%
        filter(quantity %in% c("solar_weight", "SRS")) %>%
        transmute(quantity, strategy,
                  estimate = round(estimate, 3),
                  boot_se  = round(boot_se, 3),
                  CI95 = sprintf("[%.3f, %.3f]", lwr, upr)) %>%
        as.data.frame(), row.names = FALSE)

# =====================================================================
# CASE 4.2   (p_target = 100%)
# =====================================================================
cat("\n=== Case 4.2 bootstrap (B =", B, ", L =", L_main, ", p_target = 100%) ===\n")
d2 <- load_case2_data()

opt2 <- function(SigS, SigZ, Sig, mu, psi)
  opt_case2(risk_matrix(SigS, SigZ, psi), mu, p_target = 1)
summ2 <- function(w, SigS, SigZ, Sig, mu) {
  if (anyNA(w)) return(c(solar_prop = NA, SRS = NA))
  c(solar_prop = unname(w[1] / (w[1] + w[2])),
    SRS        = srs(w, SigS, SigZ))
}

bt2 <- run_bootstrap(d2$Y, K, d2$tau, B, L_main, opt2, summ2, "Case 4.2")
print(bt2$summary %>%
        filter(quantity %in% c("solar_prop", "SRS")) %>%
        transmute(quantity, strategy,
                  estimate = round(estimate, 3),
                  boot_se  = round(boot_se, 3),
                  CI95 = sprintf("[%.3f, %.3f]", lwr, upr)) %>%
        as.data.frame(), row.names = FALSE)

main_summary <- bind_rows(bt1$summary, bt2$summary)
write_csv(main_summary, "output/bootstrap_main.csv")

# =====================================================================
# Block-length sensitivity
# =====================================================================
cat("\n=== Block-length sensitivity (B =", B_Lsens, ") ===\n")
Lsens <- map_dfr(L_grid, function(L) {
  a <- run_bootstrap(d1$Y, K, d1$tau, B_Lsens, L, opt1, summ1, "Case 4.1")$summary
  b <- run_bootstrap(d2$Y, K, d2$tau, B_Lsens, L, opt2, summ2, "Case 4.2")$summary
  bind_rows(a, b)
})
Lsens_tab <- Lsens %>%
  filter(quantity %in% c("solar_weight", "solar_prop", "SRS")) %>%
  transmute(case, quantity, strategy, L,
            CI95 = sprintf("[%.3f, %.3f]", lwr, upr), boot_se = round(boot_se, 3))
print(as.data.frame(Lsens_tab), row.names = FALSE)
write_csv(Lsens, "output/bootstrap_blocklength_sensitivity.csv")

# =====================================================================
# Figure: bootstrap distribution of solar share and SRS
# =====================================================================
draws <- bind_rows(
  bt1$draws %>% transmute(case, strategy, `Solar share (%)` = 100 * solar_weight, SRS),
  bt2$draws %>% transmute(case, strategy, `Solar share (%)` = 100 * solar_prop,   SRS)
) %>%
  pivot_longer(c(`Solar share (%)`, SRS), names_to = "quantity", values_to = "value") %>%
  mutate(strategy = factor(strategy, levels = names(strategies)))

pts <- bind_rows(bt1$summary, bt2$summary) %>%
  filter(quantity %in% c("solar_weight", "solar_prop", "SRS")) %>%
  mutate(quantity = ifelse(quantity == "SRS", "SRS", "Solar share (%)"),
         estimate = ifelse(quantity == "SRS", estimate, 100 * estimate),
         lwr      = ifelse(quantity == "SRS", lwr, 100 * lwr),
         upr      = ifelse(quantity == "SRS", upr, 100 * upr),
         strategy = factor(strategy, levels = names(strategies)))

p <- ggplot(draws, aes(strategy, value, fill = strategy)) +
  geom_violin(colour = NA, alpha = 0.35, scale = "width") +
  geom_linerange(data = pts, aes(strategy, ymin = lwr, ymax = upr),
                 inherit.aes = FALSE, linewidth = 0.6) +
  geom_point(data = pts, aes(strategy, estimate), inherit.aes = FALSE, size = 2) +
  facet_grid(quantity ~ case, scales = "free_y", switch = "y") +
  labs(x = NULL, y = NULL,
       title = "Residual moving-block bootstrap: optimal solar share and SRS",
       subtitle = sprintf("B = %d, block length L = %d; dot = point estimate, bar = percentile 95%% CI",
                          B, L_main)) +
  theme_bw() +
  theme(legend.position = "none", strip.placement = "outside",
        strip.background = element_blank())

ggsave("figures/bootstrap_weights.png", p, width = 9, height = 6, dpi = 150)
cat("\nSaved figures/bootstrap_weights.png and output/bootstrap_*.csv\n")
