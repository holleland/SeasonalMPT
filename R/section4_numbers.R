# =====================================================================
#  section4_numbers.R
#
#  Canonical Section 4.1 quantities under the *exact* QP solver
#  (minimise w' M(psi) w over {w'1 = 1, w >= 0}), replacing the coarse
#  target-grid search of the original norwegian_mix_of_solar_and_wind.R.
#  Also recomputes the Appendix F subperiod robustness check (Table 3)
#  with the same solver.
#
#  The subperiods match subperiod_robustness.R exactly: first 7 years
#  (2005-2011) and last 7 years (2013-2019); 2012 is in neither. The
#  Fourier design is built on the full sample and then row-subset, so the
#  seasonal phase stays anchored to the full series.
#
#  Run with working directory = repository root.
#  Output: output/section4_1_numbers.csv, output/table3_subperiods*.csv
# =====================================================================

source("R/robustness_helpers.R")
dir.create("output", showWarnings = FALSE)
K <- 4

strat <- c("Season-adjusted" = 0, "Balanced" = 0.5, "Seasonal" = 1)

strategy_block <- function(Y, tau, label) {
  dec <- estimate_decomposition(Y, K, tau)
  rows <- map_dfr(names(strat), function(nm) {
    psi <- strat[[nm]]
    w   <- opt_case1(risk_matrix(dec$SigS, dec$SigZ, psi))
    tibble(sample = label, strategy = nm, psi = psi,
           solar_w  = unname(w["Solar PV"]),
           wind_w   = unname(w["Wind offshore"]),
           cap_fac  = as.numeric(dec$mu %*% w),         # portfolio mean mu_Y
           SRS      = srs(w, dec$SigS, dec$SigZ),
           port_var = port_var(w, dec$Sig))
  })
  attr(rows, "dec") <- dec
  rows
}

# ---- full sample -----------------------------------------------------
d1  <- load_case1_data()
full <- strategy_block(d1$Y, d1$tau, "2005-2019")
dec  <- attr(full, "dec")

cat("=== Section 4.1: estimated inputs (K = 4, unchanged by solver) ===\n")
cat("\nmu_hat:\n"); print(round(dec$mu, 3))
cat("\nSigma_hat x100:\n");        print(round(dec$Sig  * 100, 2))
cat("\nSigma_SEASON_hat x100:\n"); print(round(dec$SigS * 100, 2))
cat("\nSigma_Z_hat x100:\n");      print(round(dec$SigZ * 100, 2))
cat(sprintf("\nrho_season = %.2f   rho_emp = %.2f   rho_Z = %.2f\n",
            seasonal_cor(dec$SigS), cov2cor(dec$Sig)[1, 2], cov2cor(dec$SigZ)[1, 2]))

cat("\n=== Section 4.1: optimal portfolios (EXACT QP) ===\n")
print(full %>%
        transmute(strategy,
                  `solar %`    = round(100 * solar_w, 1),
                  `wind %`     = round(100 * wind_w, 1),
                  `cap.factor` = round(cap_fac, 2),
                  SRS          = round(SRS, 3),
                  `SRS %`      = round(100 * SRS, 1),
                  `variance`   = round(port_var, 3)) %>%
        as.data.frame(), row.names = FALSE)

write_csv(full, "output/section4_1_numbers.csv")

# ---- subperiods (Table 3) ------------------------------------------
# Decomposition on a row-subset of the full-sample Fourier design,
# reproducing subperiod_robustness.R.
decompose_subperiod <- function(Y_full, tau, mask) {
  Cmat <- fourier_design(nrow(Y_full), K, tau)[mask, , drop = FALSE]
  Ysub <- Y_full[mask, , drop = FALSE]
  fit  <- lm(Ysub ~ 1 + Cmat)
  co   <- coef(fit); mu <- co[1, ]
  A    <- t(co[-1, , drop = FALSE])
  res  <- residuals(fit); if (is.null(dim(res))) res <- matrix(res, ncol = ncol(Ysub))
  nn   <- nrow(Ysub)
  SigS <- nn / (nn - 1) * A %*% t(A) / 2
  SigZ <- cov(res); Sig <- cov(Ysub)
  dimnames(SigS) <- dimnames(SigZ) <- dimnames(Sig) <-
    list(colnames(Y_full), colnames(Y_full))
  list(mu = mu, SigS = SigS, SigZ = SigZ, Sig = Sig)
}

strat_rows <- function(dec, label) {
  map_dfr(names(strat), function(nm) {
    psi <- strat[[nm]]
    w   <- opt_case1(risk_matrix(dec$SigS, dec$SigZ, psi))
    tibble(sample = label, strategy = nm, psi = psi,
           solar_w = unname(w["Solar PV"]), wind_w = unname(w["Wind offshore"]),
           cap_fac = as.numeric(dec$mu %*% w),
           SRS = srs(w, dec$SigS, dec$SigZ), port_var = port_var(w, dec$Sig))
  })
}

mask1 <- lubridate::year(d1$date) %in% 2005:2011
mask2 <- lubridate::year(d1$date) %in% 2013:2019
dec1  <- decompose_subperiod(d1$Y, d1$tau, mask1)
dec2  <- decompose_subperiod(d1$Y, d1$tau, mask2)

tab3 <- bind_rows(
  strat_rows(dec1, "2005-2011"),
  strat_rows(dec2, "2013-2019"),
  full
)

cov_row <- function(dec, label) tibble(
  sample = label,
  S_season = dec$SigS[1, 1], SW_season = dec$SigS[1, 2], W_season = dec$SigS[2, 2],
  S_Z      = dec$SigZ[1, 1], SW_Z      = dec$SigZ[1, 2], W_Z      = dec$SigZ[2, 2],
  S_emp    = dec$Sig[1, 1],  SW_emp    = dec$Sig[1, 2],  W_emp    = dec$Sig[2, 2])
tab3_cov <- bind_rows(
  cov_row(dec1, "2005-2011"),
  cov_row(dec2, "2013-2019"),
  cov_row(list(SigS = dec$SigS, SigZ = dec$SigZ, Sig = dec$Sig), "2005-2019"))

cat("\n=== Table 3 (Appendix F) subperiod check, EXACT QP ===\n")
cat("\ncovariance entries (x1, i.e. capacity-factor^2):\n")
print(as.data.frame(tab3_cov %>% mutate(across(-sample, ~round(., 4)))), row.names = FALSE)
cat("\nsolar weight (%) by strategy and sample:\n")
print(tab3 %>%
        transmute(sample, strategy, solar = round(100 * solar_w, 1)) %>%
        pivot_wider(names_from = sample, values_from = solar) %>%
        as.data.frame(), row.names = FALSE)

write_csv(tab3, "output/table3_subperiods.csv")
write_csv(tab3_cov, "output/table3_subperiods_cov.csv")
