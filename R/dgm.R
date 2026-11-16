# EVITAR-TB trial-level simulation: data-generating mechanism.
# Follows EVITAR-TB_Simulation_Spec_2026_11_12.docx. Values flagged PLACEHOLDER
# stand in for site/cohort data not yet available.

DAY_YR  <- 365.25
VISITS  <- (1:8) * DAY_YR / 4            # M3 ... M24 (days)
T_END   <- 2 * DAY_YR
WINDOW  <- 30

default_par <- function() list(
  # sites (protocol Table 6); large sites enrol 9 months, small 3 months
  site_n     = c(IPCG = 1250, EPJFC = 2030, EPMC = 570, PED = 2230,
                 PMEM = 890, PRSCS = 480, PEJTAA = 680, PEMSC = 520),
  site_state = c("MS", "MS", "MS", "MS", "RS", "RS", "RO", "RO"),
  site_start = c(0, 1, 3, 0, 2, 4, 3, 5),           # months; PLACEHOLDER
  site_dur   = c(9, 9, 3, 9, 3, 3, 3, 3),           # months (protocol 4.2)
  block_size = 90,                                  # PLACEHOLDER (>=30 rule)
  # outcome
  risk       = 0.03,          # annual placebo risk of incident TB onset
  state_risk = NULL,          # optional c(MS=, RS=, RO=) relative multipliers
  block_sd   = 0.35,          # SD of block log-hazard effect; PLACEHOLDER
  b_bmi      = -log(1.21),    # log HR per kg/m2 (HR 1.21 per kg/m2 decrease; TB R01 cohort)
  b_lunit    = 0.15,          # log HR per 10 Lunit points; PLACEHOLDER
  b_prevtb   = log(2),        # PLACEHOLDER
  cov_model  = "cohort",      # Lunit and previous TB from the TB R01 cohort; "placeholder": b_lunit, b_prevtb below
  hr_bmms    = 0.656,
  hr_b       = 0.656,
  effect     = "constant",
  window_eff = 1,             # effect on incident onset in days 0-30 (1 = full, 0 = none)
  nonadh     = 0,             # share of active-arm participants with low intake
  nonadh_eff = 0.3,           # fraction of log-HR retained if low intake
  contam     = 0,             # share of placebo receiving shared capsules
  contam_eff = 0.5,           # fraction of B+MMS log-HR they receive
  # window-period (prevalent/incipient) cases, no treatment effect
  prev_frac  = 0.006,         # PLACEHOLDER
  prev_leak  = 0,             # share diagnosed at M3 instead of within 30 days
  # detection
  symp_mean  = 180,           # mean days from onset to symptomatic presentation
  visit_sens = 0.85,          # detection at a scheduled visit
  # loss to follow-up (release/transfer out of network)
  loss       = 0.33,          # annual
  loss_shape = "constant",    # "constant" or "front" (higher in first 6 months)
  loss_risk  = 0,             # log rate ratio per SD of prognostic score
  loss_diff  = 0,             # absolute increase in annual loss in B+MMS arm
  enrol_delay = 0             # months added to every site's enrolment duration
)

# --- calibration of the onset hazard so marginal placebo annual risk = par$risk
calibrate_h0 <- function(par, n = 200000, seed = 1) {
  set.seed(seed)
  cv <- gen_covariates(n, par)
  m  <- exp(cv$lp) * exp(rnorm(n, -par$block_sd^2 / 2, par$block_sd)) * cv$smult
  target <- 1 - (1 - par$risk)^2
  uniroot(function(h) mean(1 - exp(-h * m * 2)) - target, c(1e-5, 1))$root
}

# Aggregate cohort inputs (TB R01 prison cohort; see inputs/README): Lunit score
# quantiles by previous TB, and Cox coefficients for a natural spline of Lunit
# (knots 6.65, 27.68, 49.48, 92.63) interacted with previous TB.
COHORT <- if (file.exists("../inputs/cohort_covariates.rds"))
  readRDS("../inputs/cohort_covariates.rds") else NULL    # scripts run from R/

cohort_lp_lunit <- function(lunit, prevtb, co) {
  kn <- co$knots
  B <- splines::ns(pmin(lunit, 100), knots = kn[2:3], Boundary.knots = kn[c(1, 4)])
  b <- co$coef
  drop(B %*% b[1:3]) + b[4] * prevtb + drop((B * prevtb) %*% b[5:7])
}

gen_covariates <- function(n, par, prison = NULL) {
  bmi    <- pmin(pmax(rnorm(n, 23.5, 3.2), 16), 29.9)
  age    <- sample(1:3, n, TRUE, c(0.45, 0.35, 0.20))
  previm <- rbinom(n, 1, 0.55)
  smoke  <- sample(1:3, n, TRUE, c(0.55, 0.40, 0.05))
  hiv    <- sample(1:3, n, TRUE, c(0.015, 0.935, 0.05))
  if (identical(par$cov_model, "cohort")) {
    # previous TB and Lunit from the cohort; Lunit drawn from its distribution
    # within previous-TB group; risk from the cohort spline-by-previous-TB model
    co <- COHORT
    prevtb <- rbinom(n, 1, co$prev_tb)
    p <- seq(0, 1, 0.005); u <- runif(n)
    lunit <- ifelse(prevtb == 1, approx(p, co$lunit_q_prev, u)$y, approx(p, co$lunit_q_noprev, u)$y)
    lp <- par$b_bmi * (bmi - 23.5) + cohort_lp_lunit(lunit, prevtb, co)
  } else {
    # PLACEHOLDER joint distribution (original runs)
    hi     <- runif(n) < 0.15
    lunit  <- ifelse(hi, runif(n, 40, 100), 100 * rbeta(n, 1, 9))
    prevtb <- rbinom(n, 1, 0.08)
    lp <- par$b_bmi * (bmi - 23.5) + par$b_lunit * (lunit - 10) / 10 + par$b_prevtb * prevtb
  }
  lp <- lp - log(mean(exp(lp)))
  smult <- 1
  if (!is.null(prison) && !is.null(par$state_risk)) {
    smult <- par$state_risk[par$site_state[prison]]
  } else if (is.null(prison) && !is.null(par$state_risk)) {
    w <- par$site_n / sum(par$site_n)
    smult <- par$state_risk[par$site_state[sample(seq_along(w), n, TRUE, w)]]
  }
  data.frame(bmi, lunit, prevtb, age, previm, smoke, hiv, lp, smult = unname(smult))
}

# --- minimization (protocol Section 6.4, Table 13)
minimize_alloc <- function(d) {
  f <- data.frame(block = d$block,
                  bmi = cut(d$bmi, c(-Inf, 18.5, 25, Inf), labels = FALSE),
                  lun = 1 + (d$lunit >= 40), age = d$age, prevtb = d$prevtb + 1,
                  previm = d$previm + 1, smoke = d$smoke, hiv = d$hiv)
  f$block <- as.integer(factor(f$block))
  grp <- integer(nrow(d))
  for (p in unique(d$prison)) {
    idx <- which(d$prison == p)
    idx <- idx[order(d$enrol[idx])]
    cnt <- lapply(f, function(x) matrix(0L, max(x), 3))
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k <= 15) {
        g <- sample.int(3, 1)
      } else {
        sc <- numeric(3)
        for (j in seq_along(f)) sc <- sc + cnt[[j]][f[[j]][i], ]
        mn <- which(sc == min(sc))
        pr <- if (length(mn) == 1) { x <- rep(0.1, 3); x[mn] <- 0.8; x }
              else if (length(mn) == 2) { x <- rep(0.1, 3); x[mn] <- 0.45; x }
              else rep(1 / 3, 3)
        g <- sample.int(3, 1, prob = pr)
      }
      grp[i] <- g
      for (j in seq_along(f)) cnt[[j]][f[[j]][i], g] <- cnt[[j]][f[[j]][i], g] + 1L
    }
  }
  grp   # 1 placebo, 2 B, 3 B+MMS
}

# hazard-ratio multiplier for onset in each period
PERIODS <- c(0, WINDOW, DAY_YR / 4, DAY_YR / 2, 3 * DAY_YR / 4, DAY_YR, 2 * DAY_YR + 365)
period_loghr <- function(loghr, effect, window_eff = 1) {
  # columns: [0,30], (30,M3], (M3,M6], (M6,M9], (M9,M12], (M12, end]
  w <- switch(effect,
    constant      = c(0, 1, 1, 1, 1, 1),
    wane          = c(0, 1, 1, 0.75, 0.25, 0),
    half_after_m6 = c(0, 1, 1, 0.5, 0.5, 0.5),
    delayed_m3    = c(0, 0, 1, 1, 1, 1))
  if (effect != "delayed_m3") w[1] <- window_eff
  outer(loghr, w)
}

gen_trial <- function(par, h0) {
  ns <- par$site_n; P <- length(ns)
  prison <- rep(seq_len(P), ns)
  n <- length(prison)
  cv <- do.call(rbind, lapply(seq_len(P), function(p) gen_covariates(ns[p], par, p)))
  # enrolment date (days from trial start)
  dur <- (par$site_dur + par$enrol_delay) * DAY_YR / 12
  enrol <- par$site_start[prison] * DAY_YR / 12 + runif(n) * dur[prison]
  # blocks within prison
  nb <- pmax(1, round(ns / par$block_size))
  block <- paste0(prison, "_", unlist(lapply(seq_len(P), function(p) sample.int(nb[p], ns[p], TRUE))))
  ub <- tapply(seq_along(block), block, length)
  beff <- setNames(rnorm(length(ub), -par$block_sd^2 / 2, par$block_sd), names(ub))
  d <- data.frame(id = seq_len(n), prison, block, enrol, cv)
  d$grp <- minimize_alloc(d)
  # individual treatment log-HR (dilution: non-adherence, contamination)
  lhr <- c(0, log(par$hr_b), log(par$hr_bmms))[d$grp]
  if (par$nonadh > 0) {
    low <- d$grp > 1 & runif(n) < par$nonadh
    lhr[low] <- lhr[low] * par$nonadh_eff
  }
  if (par$contam > 0) {
    ct <- d$grp == 1 & runif(n) < par$contam
    lhr[ct] <- log(par$hr_bmms) * par$contam_eff
  }
  rate <- h0 * exp(d$lp) * d$smult * exp(beff[d$block]) / DAY_YR      # per day
  plh  <- period_loghr(lhr, par$effect, par$window_eff)
  len  <- diff(PERIODS)
  cumH <- t(apply(rate * exp(plh) * matrix(len, n, length(len), byrow = TRUE), 1, cumsum))
  E <- rexp(n)
  k <- max.col(cbind(E <= cumH, TRUE), ties.method = "first")       # period index (7 = none)
  onset <- rep(Inf, n)
  ok <- k <= length(len)
  prevH <- cbind(0, cumH)[cbind(which(ok), k[ok])]
  onset[ok] <- PERIODS[k[ok]] + (E[ok] - prevH) /
               (rate[ok] * exp(plh[cbind(which(ok), k[ok])]))
  # detection: symptomatic presentation vs next scheduled visit
  symp <- onset + rexp(n, 1 / par$symp_mean)
  vdet <- rep(Inf, n)
  for (v in rev(VISITS)) {
    hit <- onset < v & runif(n) < par$visit_sens
    vdet[hit] <- v
  }
  diag <- pmin(symp, vdet)
  # window-period (prevalent) cases: diagnosed within 30 days, or leak to M3
  pv <- runif(n) < par$prev_frac
  pdiag <- ifelse(runif(n) < par$prev_leak, VISITS[1], runif(n, 1, WINDOW))
  diag[pv] <- pmin(diag[pv], pdiag[pv])
  # loss to follow-up
  ann <- par$loss + par$loss_diff * (d$grp == 3)
  r <- -log(1 - ann) * exp(par$loss_risk * (d$lp - mean(d$lp)) / sd(d$lp))
  r <- r / mean(exp(par$loss_risk * (d$lp - mean(d$lp)) / sd(d$lp)))
  if (par$loss_shape == "constant") {
    cens <- rexp(n, r) * DAY_YR
  } else {
    # front-loaded: same 2-year retention, half of the 2-year cumulative
    # hazard falls in the first 6 months
    tot <- 2 * r
    r1 <- (0.5 * tot) / 0.5; r2 <- (0.5 * tot) / 1.5
    u <- rexp(n)
    cens <- ifelse(u <= r1 * 0.5, u / r1, 0.5 + (u - r1 * 0.5) / r2) * DAY_YR
  }
  d$diag <- diag
  d$cens <- pmin(cens, T_END)
  d
}
