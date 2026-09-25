# =====================================================================
#  robustness_helpers.R
#  Shared data-loading + estimation + optimisation helpers used by
#    R/robustness_K.R
#    R/bootstrap_weights.R
#
#  The data pipeline mirrors exactly
#    R/norwegian_mix_of_solar_and_wind.R              (Section 4.1)
#    R/norwegian_mix_of_solar_and_wind_with_demand.R  (Section 4.2)
#  so that the K = 4 results reproduce Table/Figure 3 and Section 4.2.
#
#  Run with working directory = repository root.
# =====================================================================

suppressWarnings(suppressMessages({
  library(tidyverse)
  library(lubridate)
  library(quadprog)
  library(forecast)
}))

# ---------------------------------------------------------------------
# 1. Data loading
# ---------------------------------------------------------------------

load_case1_data <- function() {
  ## ---- Solar PV (mean CF across the five locations) ----------------
  solar_files <- list.files("data/solar/", full.names = TRUE)
  PV <- tibble()
  for (file in solar_files) {
    PV <- bind_rows(
      PV,
      suppressWarnings(readr::read_csv(file, skip = 10, show_col_types = FALSE)) %>%
        mutate(datetime = as.POSIXct(time, format = "%Y%m%d:%H%M")) %>%
        mutate(CF = P / 1000, locID = "Solar PV") %>%
        select(datetime, locID, CF) %>%
        filter(year(datetime) %in% 2005:2019)
    )
  }
  PV <- PV %>% group_by(datetime, locID) %>% summarize(CF = mean(CF, na.rm = TRUE), .groups = "drop")

  ## ---- Offshore wind portfolio (Table A1 of Hoelleland et al. 2025) -
  wind <- readRDS("data/NVE.rds") %>%
    select(locID, value, datetime) %>%
    filter(year(datetime) %in% 2005:2019)

  port <- tibble(
    locID   = 1:20,
    turbines = c(269, 146, 100, 208, 204, 41, 102, 194, 37, 0,
                 0, 63, 0, 25, 0, 0, 218, 32, 28, 333)
  ) %>%
    left_join(wind, by = "locID") %>%
    group_by(datetime) %>%
    summarize(CF = sum(turbines * value) / 2000, .groups = "drop") %>%
    mutate(locID = "Wind offshore")

  power <- bind_rows(port, PV %>% filter(!is.na(datetime))) %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(locID, date) %>%
    summarize(CF = mean(CF, na.rm = TRUE), .groups = "drop") %>%
    filter(year(date) %in% 2005:2019)

  wide <- power %>%
    pivot_wider(names_from = locID, values_from = CF) %>%
    arrange(date)

  list(
    date = wide$date,
    Y    = as.matrix(wide[, c("Solar PV", "Wind offshore")]),
    tau  = 365.25
  )
}

load_case2_data <- function(size_of_solar_park = 10000) {
  solar_files <- list.files("data/solar/", full.names = TRUE)
  PV <- tibble()
  for (file in solar_files) {
    PV <- bind_rows(
      PV,
      suppressWarnings(readr::read_csv(file, skip = 10, show_col_types = FALSE)) %>%
        mutate(datetime = as.POSIXct(time, format = "%Y%m%d:%H%M")) %>%
        mutate(CF = P / 1000, locID = "Solar PV") %>%
        select(datetime, locID, CF) %>%
        filter(year(datetime) %in% 2005:2019)
    )
  }
  PV <- PV %>% group_by(datetime, locID) %>% summarize(CF = mean(CF, na.rm = TRUE), .groups = "drop")

  wind <- readRDS("data/NVE.rds") %>%
    select(locID, value, datetime) %>%
    filter(year(datetime) %in% 2005:2019)

  port <- tibble(
    locID   = 1:20,
    turbines = c(269, 146, 100, 208, 204, 41, 102, 194, 37, 0,
                 0, 63, 0, 25, 0, 0, 218, 32, 28, 333)
  ) %>%
    left_join(wind, by = "locID") %>%
    group_by(datetime) %>%
    summarize(CF = sum(turbines * value) / 2000, .groups = "drop") %>%
    mutate(locID = "WP portfolio")

  power <- bind_rows(port, PV %>% filter(!is.na(datetime))) %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(locID, date) %>%
    summarize(CF = mean(CF, na.rm = TRUE), .groups = "drop") %>%
    filter(year(date) %in% 2005:2019)

  ## ---- Statnett consumption ---------------------------------------
  consumption <- lapply(list.files("data/statnett", full.names = TRUE), readr::read_csv2) %>%
    bind_rows() %>%
    mutate(datetime = as.POSIXct(`Time(Local)`, format = "%d.%m.%Y %H:%M:%S")) %>%
    select(datetime, Consumption)

  daily_consumption <- consumption %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(date) %>%
    summarize(Consumption = sum(Consumption, na.rm = TRUE), .groups = "drop")

  power <- power %>%
    pivot_wider(names_from = "locID", values_from = "CF") %>%
    mutate(`Solar PV`     = `Solar PV` * 2000 * 24 * size_of_solar_park / 1e6,
           `WP portfolio` = `WP portfolio` * 15e6 * 24 / 1e6) %>%
    left_join(daily_consumption, by = "date") %>%
    filter(year(date) > 2005) %>%
    as_tibble() %>%
    filter(Consumption > 2e5) %>%
    na.omit() %>%
    arrange(date)

  ## detrend consumption to the 2019 level (linear growth removed)
  con.lm <- lm(Consumption ~ t, data = power %>% mutate(t = 1:n()))
  power$Consumption <- power$Consumption +
    coef(con.lm)[2] * nrow(power) - coef(con.lm)[2] * (1:nrow(power))

  list(
    date              = power$date,
    Y                 = as.matrix(power[, c("Solar PV", "WP portfolio", "Consumption")]),
    daily_consumption = daily_consumption,
    tau               = 365.25
  )
}

# ---------------------------------------------------------------------
# 2. Fourier design + covariance decomposition
# ---------------------------------------------------------------------

# Fourier design matrix of order K for n observations, period tau.
# Delegates to forecast::fourier(msts(...)) so the basis is byte-identical to
# norwegian_mix_of_solar_and_wind{,_with_demand}.R, then reorders the columns
# to (C1..CK, S1..SK) as those scripts do to the fitted A matrix. The covariance
# estimates are invariant to that permutation, so mu/SigS/SigZ/Sig reproduce the
# manuscript exactly for the full sample and for any subperiod (time restarts
# at 1, matching msts() built from the subsetted data).
fourier_design <- function(n, K, tau = 365.25) {
  x    <- forecast::msts(matrix(0, n, 1), seasonal.periods = tau)
  cols <- forecast::fourier(x, K = K)                       # order S1,C1,S2,C2,...
  cols <- cols[, c(seq(2, 2 * K, 2), seq(1, 2 * K - 1, 2)), drop = FALSE]
  colnames(cols) <- paste0(rep(c("C", "S"), each = K), rep(1:K, 2))
  cols
}

# Estimate mu, A, residuals and the three covariance matrices for a given K.
#   Y     : n x m data matrix
#   K     : Fourier order
# Returns list(mu, A, resid, SigS, SigZ, Sig, Cmat)
estimate_decomposition <- function(Y, K, tau = 365.25) {
  n <- nrow(Y); m <- ncol(Y)
  Cmat <- fourier_design(n, K, tau)
  fit  <- lm(Y ~ 1 + Cmat)
  co   <- coef(fit)
  mu   <- co[1, ]
  A    <- t(co[-1, , drop = FALSE])          # m x 2K, columns already cos|sin ordered
  res  <- residuals(fit)
  if (is.null(dim(res))) res <- matrix(res, ncol = m)
  SigS <- n / (n - 1) * A %*% t(A) / 2
  SigZ <- cov(res)
  Sig  <- cov(Y)
  dimnames(SigS) <- dimnames(SigZ) <- dimnames(Sig) <- list(colnames(Y), colnames(Y))
  list(mu = mu, A = A, resid = res, SigS = SigS, SigZ = SigZ, Sig = Sig, Cmat = Cmat)
}

# ---------------------------------------------------------------------
# 3. Portfolio optimisation
# ---------------------------------------------------------------------

# psi-weighted risk matrix
risk_matrix <- function(SigS, SigZ, psi) psi * SigS + (1 - psi) * SigZ

# --- Case 4.1: two assets, long-only, weights sum to one, no mean target
#     (any target mean pins the 2-asset mix, exactly as in Section 4.1).
#     We minimise w' M w over the simplex -> this is the min-risk-measure
#     portfolio the paper's arrows point to.
opt_case1 <- function(M) {
  # minimise w'Mw s.t. sum w = 1, w >= 0
  Dmat <- (M + t(M)) / 2
  # guard against numerical non-PD
  ev <- min(eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 1e-10) Dmat <- Dmat + diag(nrow(Dmat)) * (1e-10 - ev)
  Amat <- cbind(rep(1, 2), diag(2))
  bvec <- c(1, 0, 0)
  sol  <- solve.QP(Dmat, rep(0, 2), Amat, bvec, meq = 1)
  w <- sol$solution
  names(w) <- colnames(M)
  w
}

# --- Case 4.1, manuscript version: reproduce norwegian_mix_of_solar_and_wind.R
#     exactly. MPT() solves, for a given target mean, the equality-constrained
#     QP (w'1 = 1, w'mu = target, w >= 0); the minimum-risk portfolio is then
#     the grid argmin over target in seq(min(mu)+step, max(mu)-step, by = step).
#     `covmat` is Sigma_Z (season-adjusted), Sigma (balanced) or
#     Sigma_SEASON (seasonal), matching the three `cov =` arguments in the script.
mpt_qp_case1 <- function(mu, covmat, target) {
  Amat <- cbind(1, mu, diag(length(mu)))
  bvec <- c(1, target, rep(0, length(mu)))
  solve.QP(Dmat = covmat, dvec = rep(0, length(mu)), Amat = Amat, bvec = bvec, meq = 2)
}

opt_case1_grid <- function(covmat, mu, step = 0.01) {
  targets <- seq(min(mu) + step, max(mu) - step, by = step)
  vals <- vapply(targets, function(tg)
    tryCatch(mpt_qp_case1(mu, covmat, tg)$value, error = function(e) NA_real_),
    numeric(1))
  if (all(is.na(vals))) return(setNames(rep(NA_real_, length(mu)), colnames(covmat)))
  w <- mpt_qp_case1(mu, covmat, targets[which.min(vals)])$solution
  names(w) <- colnames(covmat)
  w
}

# Pick the covariance matrix used by strategy `psi` in the manuscript's Case 4.1
# (psi = 0 -> Sigma_Z, psi = 1 -> Sigma_SEASON, else the empirical Sigma).
covmat_case1 <- function(SigS, SigZ, Sig, psi)
  if (psi == 0) SigZ else if (psi == 1) SigS else Sig

# --- Case 4.2: three assets (Solar, Wind, Demand), demand weight fixed = -1,
#     w_D' mu = (p_target - 1) * mu_demand, solar/wind >= 0.
#     Mirrors MPT() in norwegian_mix_of_solar_and_wind_with_demand.R.
opt_case2 <- function(M, mu, p_target = 1) {
  Dmat <- (M + t(M)) / 2
  ev <- min(eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 1e-10) Dmat <- Dmat + diag(nrow(Dmat)) * (1e-10 - ev)
  Amat <- cbind(c(0, 0, 1), mu, c(1, 0, 0), c(0, 1, 0))
  bvec <- c(-1, (p_target - 1) * mu[3], 0, 0)
  sol  <- solve.QP(Dmat, rep(0, 3), Amat, bvec, meq = 2)
  w <- sol$solution
  names(w) <- names(mu)
  w
}

# Seasonal Risk Score for a weight vector
srs <- function(w, SigS, SigZ)
  as.numeric(t(w) %*% SigS %*% w / (t(w) %*% (SigS + SigZ) %*% w))

port_var <- function(w, Sig) as.numeric(t(w) %*% Sig %*% w)

seasonal_cor <- function(SigS) {
  d <- sqrt(diag(SigS))
  SigS[1, 2] / (d[1] * d[2])
}

# ---------------------------------------------------------------------
# 4. One-shot summaries for the three benchmark strategies
# ---------------------------------------------------------------------

# Case 4.1: returns a tibble row per strategy
summarise_case1 <- function(dec) {
  strategies <- c("Season-adjusted" = 0, "Balanced" = 0.5, "Seasonal" = 1)
  map_dfr(names(strategies), function(nm) {
    psi <- strategies[[nm]]
    M   <- risk_matrix(dec$SigS, dec$SigZ, psi)
    w   <- opt_case1(M)
    tibble(
      strategy    = nm,
      psi         = psi,
      solar_weight = unname(w["Solar PV"]),
      wind_weight  = unname(w["Wind offshore"]),
      SRS          = srs(w, dec$SigS, dec$SigZ),
      port_var     = port_var(w, dec$Sig)
    )
  })
}

# Case 4.2: returns a tibble row per strategy at a given p_target
summarise_case2 <- function(dec, p_target = 1) {
  strategies <- c("Season-adjusted" = 0, "Balanced" = 0.5, "Seasonal" = 1)
  map_dfr(names(strategies), function(nm) {
    psi <- strategies[[nm]]
    M   <- risk_matrix(dec$SigS, dec$SigZ, psi)
    w   <- tryCatch(opt_case2(M, dec$mu, p_target), error = function(e) rep(NA_real_, 3))
    solar_prop <- if (anyNA(w)) NA_real_ else w[1] / (w[1] + w[2])
    tibble(
      strategy      = nm,
      psi           = psi,
      p_target      = p_target,
      solar_units   = w[1],
      wind_units    = w[2],
      solar_prop    = as.numeric(solar_prop),
      SRS           = if (anyNA(w)) NA_real_ else srs(w, dec$SigS, dec$SigZ)
    )
  })
}
