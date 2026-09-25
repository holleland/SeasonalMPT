# =====================================================================
#  psi_sensitivity.R
#
#  Sensitivity analysis to support practical guidance on choosing psi.
#  Sweeps the preference parameter psi over [0, 1] in small steps and
#  records the main outputs along the psi-indexed optimal portfolios,
#  for both empirical applications:
#     Section 4.1  capacity-factor mix   -> solar weight, SRS,
#                                           portfolio capacity factor,
#                                           seasonal / residual variance
#     Section 4.2  net-balance mix       -> solar share, SRS,
#                  (p_target = 100%)       net-balance standard deviation
#
#  Run with working directory = repository root.
#  Outputs: output/psi_sensitivity_*.csv, figures/psi_sensitivity.png
# =====================================================================

source("R/robustness_helpers.R")
dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

K        <- 4
psi_grid <- seq(0, 1, by = 0.01)

# ---------------------------------------------------------------------
# Case 4.1: capacity-factor mix
# ---------------------------------------------------------------------
d1  <- load_case1_data()
dec1 <- estimate_decomposition(d1$Y, K, d1$tau)

res1 <- map_dfr(psi_grid, function(psi) {
  w <- opt_case1(risk_matrix(dec1$SigS, dec1$SigZ, psi))
  tibble(psi          = psi,
         solar        = unname(w["Solar PV"]),
         SRS          = srs(w, dec1$SigS, dec1$SigZ),
         cap_factor   = as.numeric(dec1$mu %*% w),
         seasonal_var = as.numeric(t(w) %*% dec1$SigS %*% w),
         resid_var    = as.numeric(t(w) %*% dec1$SigZ %*% w),
         total_var    = port_var(w, dec1$Sig))
})
write_csv(res1, "output/psi_sensitivity_case1.csv")

# ---------------------------------------------------------------------
# Case 4.2: net-balance mix with demand, p_target = 100%
#   net-balance SD = sqrt(w_D' Sigma w_D) since E[net balance] = 0 there.
# ---------------------------------------------------------------------
d2  <- load_case2_data()
dec2 <- estimate_decomposition(d2$Y, K, d2$tau)

res2 <- map_dfr(psi_grid, function(psi) {
  w <- tryCatch(opt_case2(risk_matrix(dec2$SigS, dec2$SigZ, psi), dec2$mu, 1),
                error = function(e) rep(NA_real_, 3))
  if (anyNA(w)) return(tibble(psi = psi, solar = NA, SRS = NA, nb_sd = NA))
  tibble(psi   = psi,
         solar = unname(w[1] / (w[1] + w[2])),
         SRS   = srs(w, dec2$SigS, dec2$SigZ),
         nb_sd = sqrt(port_var(w, dec2$Sig)))
})
write_csv(res2, "output/psi_sensitivity_case2.csv")

# ---------------------------------------------------------------------
# Text summary at psi = 0, 0.25, 0.5, 0.75, 1
# ---------------------------------------------------------------------
keys <- c(0, 0.25, 0.5, 0.75, 1)
cat("\n=== Case 4.1 (capacity-factor mix) ===\n")
print(res1 %>% filter(psi %in% keys) %>%
        transmute(psi,
                  `solar %`   = round(100 * solar, 1),
                  SRS         = round(SRS, 3),
                  `cap.fac`   = round(cap_factor, 3),
                  `seas.var`  = round(seasonal_var, 5),
                  `resid.var` = round(resid_var, 5)) %>% as.data.frame(),
      row.names = FALSE)

cat("\n=== Case 4.2 (net-balance mix, p_target = 100%) ===\n")
print(res2 %>% filter(psi %in% keys) %>%
        transmute(psi,
                  `solar %` = round(100 * solar, 1),
                  SRS       = round(SRS, 3),
                  `nb.sd`   = round(nb_sd, 1)) %>% as.data.frame(),
      row.names = FALSE)

# elbow of the residual-vs-seasonal trade-off (case 4.1): psi where the
# marginal residual-variance cost of further seasonal-variance reduction
# starts to climb steeply (max curvature of the trade-off curve).
tc <- res1 %>% arrange(psi)
dz <- c(NA, diff(tc$resid_var)); ds <- c(NA, diff(tc$seasonal_var))
slope <- dz / ds                      # d(resid) / d(seasonal), <0
tc$dslope <- c(NA, diff(slope))
cat(sprintf("\nCase 4.1 trade-off: |d resid / d seasonal| rises past psi ~ %.2f\n",
            tc$psi[which.max(abs(tc$dslope))]))

# ---------------------------------------------------------------------
# Figure: three shared quantities x two applications.
#   Row 3 is total variability expressed relative to its psi = 0 value,
#   so the two applications are on a common (dimensionless) scale.
# ---------------------------------------------------------------------
plot_df <- bind_rows(
  res1 %>% transmute(psi, case = "Case 4.1: capacity-factor mix",
                     `Solar share (%)` = 100 * solar, SRS,
                     `Total SD (relative to psi = 0)` =
                       sqrt(total_var) / sqrt(total_var[psi == 0])),
  res2 %>% transmute(psi, case = "Case 4.2: net-balance mix (p = 100%)",
                     `Solar share (%)` = 100 * solar, SRS,
                     `Total SD (relative to psi = 0)` = nb_sd / nb_sd[psi == 0])
) %>%
  pivot_longer(c(`Solar share (%)`, SRS, `Total SD (relative to psi = 0)`),
               names_to = "quantity", values_to = "value") %>%
  mutate(quantity = factor(quantity,
           levels = c("Solar share (%)", "SRS",
                      "Total SD (relative to psi = 0)")))

p <- ggplot(plot_df, aes(psi, value)) +
  geom_vline(xintercept = c(0, 0.5, 1), linetype = 3, colour = "grey65") +
  geom_line(linewidth = 0.7) +
  facet_grid(quantity ~ case, scales = "free_y", switch = "y") +
  scale_x_continuous(expression(psi), breaks = seq(0, 1, 0.25)) +
  labs(y = NULL) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("figures/psi_sensitivity.pdf", p, width = 9, height = 5.8)
cat("\nSaved figures/psi_sensitivity.pdf and output/psi_sensitivity_*.csv\n")
