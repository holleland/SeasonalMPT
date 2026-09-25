# =====================================================================
#  robustness_K.R
#
#  Checks whether the optimal energy mix and SRS remain stable under
#  alternative Fourier orders, to support the choice of K.
#
#  Re-runs the full estimation + optimisation pipeline of Sections 4.1
#  and 4.2 for a range of Fourier orders K and reports, for each of the
#  three benchmark strategies (season-adjusted psi=0, balanced psi=0.5,
#  seasonal psi=1):
#     * the optimal energy mix (solar share)
#     * the Seasonal Risk Score (SRS)
#     * the portfolio variance (case 4.1) and the p_target at which the
#       30 GW offshore-wind target is reached (case 4.2)
#     * the seasonal / empirical / residual correlations
#     * the per-asset AIC difference from its minimum (extends Table 2)
#
#  Run with working directory = repository root.
#  Outputs: output/robustness_K_*.csv  and  figures/robustness_K.png
# =====================================================================

source("R/robustness_helpers.R")
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

K_grid <- 2:8

# ---------------------------------------------------------------------
# AIC per asset for each K  (extends Table 2 of the manuscript)
# ---------------------------------------------------------------------
aic_table <- function(Y, K_grid, tau, asset_names = colnames(Y)) {
  n <- nrow(Y)
  out <- map_dfr(K_grid, function(K) {
    Cmat <- fourier_design(n, K, tau)
    tibble(K = K,
           asset = asset_names,
           AIC = sapply(seq_len(ncol(Y)), function(j) AIC(lm(Y[, j] ~ 1 + Cmat))))
  })
  out %>% group_by(asset) %>% mutate(dAIC = AIC - min(AIC)) %>% ungroup()
}

# =====================================================================
# CASE 4.1 : solar / wind capacity-factor mix
# =====================================================================
cat("\n=== Case 4.1: capacity-factor mix ===\n")
d1 <- load_case1_data()

res1 <- map_dfr(K_grid, function(K) {
  dec <- estimate_decomposition(d1$Y, K, d1$tau)
  summarise_case1(dec) %>%
    mutate(K = K,
           rho_season = seasonal_cor(dec$SigS),
           rho_emp    = cov2cor(dec$Sig)[1, 2],
           rho_Z      = cov2cor(dec$SigZ)[1, 2],
           seas_share_solar = dec$SigS[1, 1] / dec$Sig[1, 1],
           seas_share_wind  = dec$SigS[2, 2] / dec$Sig[2, 2],
           .before = 1)
})

aic1 <- aic_table(d1$Y, K_grid, d1$tau)

cat("\n-- optimal solar weight (%) by K and strategy --\n")
print(res1 %>%
        transmute(K, strategy, solar_weight = round(100 * solar_weight, 1)) %>%
        pivot_wider(names_from = strategy, values_from = solar_weight) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- Seasonal Risk Score by K and strategy --\n")
print(res1 %>%
        transmute(K, strategy, SRS = round(SRS, 3)) %>%
        pivot_wider(names_from = strategy, values_from = SRS) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- portfolio variance (x100) by K and strategy --\n")
print(res1 %>%
        transmute(K, strategy, var = round(100 * port_var, 3)) %>%
        pivot_wider(names_from = strategy, values_from = var) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- correlations by K --\n")
print(res1 %>% distinct(K, rho_season, rho_emp, rho_Z) %>%
        mutate(across(starts_with("rho"), ~round(., 2))) %>% as.data.frame(),
      row.names = FALSE)

cat("\n-- AIC difference from minimum, per asset (extends Table 2) --\n")
print(aic1 %>% transmute(K, asset, dAIC = round(dAIC, 1)) %>%
        pivot_wider(names_from = asset, values_from = dAIC) %>% as.data.frame(),
      row.names = FALSE)

write_csv(res1, "output/robustness_K_case1.csv")
write_csv(aic1, "output/robustness_K_case1_aic.csv")

# =====================================================================
# CASE 4.2 : net-balance mix with demand  (p_target = 100%)
# =====================================================================
cat("\n\n=== Case 4.2: net-balance mix with demand (p_target = 100%) ===\n")
d2 <- load_case2_data()

# p_target at which offshore wind reaches 2000 turbines (~30 GW).
# Matches norwegian_mix_of_solar_and_wind_with_demand.R: first value on the
# 1%-spaced p_target grid at which the wind-unit count is >= 2000.
p_at_30GW <- function(dec, psi, grid = seq(0.01, 6, by = 0.01)) {
  M <- risk_matrix(dec$SigS, dec$SigZ, psi)
  w2 <- vapply(grid, function(p)
    tryCatch(opt_case2(M, dec$mu, p)[2], error = function(e) NA_real_), numeric(1))
  hit <- which(w2 >= 2000)
  if (length(hit) == 0) NA_real_ else grid[hit[1]]
}

res2 <- map_dfr(K_grid, function(K) {
  dec <- estimate_decomposition(d2$Y, K, d2$tau)
  base <- summarise_case2(dec, p_target = 1)
  psis <- c("Season-adjusted" = 0, "Balanced" = 0.5, "Seasonal" = 1)
  base %>%
    mutate(K = K,
           p_target_30GW = map_dbl(strategy, ~p_at_30GW(dec, psis[[.x]])),
           rho_cons_solar_season = cov2cor(dec$SigS)[1, 3],
           rho_cons_wind_season  = cov2cor(dec$SigS)[2, 3],
           .before = 1)
})

cat("\n-- optimal solar-PV share (%) by K and strategy --\n")
print(res2 %>%
        transmute(K, strategy, solar_prop = round(100 * solar_prop, 1)) %>%
        pivot_wider(names_from = strategy, values_from = solar_prop) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- Seasonal Risk Score by K and strategy --\n")
print(res2 %>%
        transmute(K, strategy, SRS = round(SRS, 3)) %>%
        pivot_wider(names_from = strategy, values_from = SRS) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- p_target (%) at which 30 GW offshore wind is reached --\n")
print(res2 %>%
        transmute(K, strategy, p30 = round(100 * p_target_30GW)) %>%
        pivot_wider(names_from = strategy, values_from = p30) %>%
        as.data.frame(), row.names = FALSE)

cat("\n-- seasonal correlations with consumption by K --\n")
print(res2 %>% distinct(K, rho_cons_solar_season, rho_cons_wind_season) %>%
        mutate(across(starts_with("rho"), ~round(., 2))) %>% as.data.frame(),
      row.names = FALSE)

write_csv(res2, "output/robustness_K_case2.csv")

# =====================================================================
# Figure: solar share and SRS vs K for both cases
# =====================================================================
plot_df <- bind_rows(
  res1 %>% transmute(case = "Case 4.1: capacity factor", K, strategy,
                     `Solar share` = 100 * solar_weight, SRS),
  res2 %>% transmute(case = "Case 4.2: net balance (p=100%)", K, strategy,
                     `Solar share` = 100 * solar_prop, SRS)
) %>%
  pivot_longer(c(`Solar share`, SRS), names_to = "quantity", values_to = "value") %>%
  mutate(strategy = factor(strategy,
                           levels = c("Season-adjusted", "Balanced", "Seasonal")))

p <- ggplot(plot_df, aes(K, value, colour = strategy)) +
  geom_vline(xintercept = 4, linetype = 3, colour = "grey50") +
  geom_line() + geom_point() +
  facet_grid(quantity ~ case, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = K_grid) +
  labs(x = "Fourier order K", y = NULL,
       title = "Sensitivity of the optimal mix and SRS to the Fourier order K",
       subtitle = "Dotted line: K = 4 used in the manuscript") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        strip.placement = "outside", strip.background = element_blank())

ggsave("figures/robustness_K.png", p, width = 9, height = 6, dpi = 150)
cat("\nSaved figures/robustness_K.png and output/robustness_K_*.csv\n")
