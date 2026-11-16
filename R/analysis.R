# EVITAR-TB primary analysis (protocol Section 10.4) applied to one dataset.
# Discrete-time cloglog model, quarterly intervals from day 30, log-exposure
# offset, 6-month baseline-hazard periods, prison, natural splines for BMI and
# Lunit (knots at fixed pooled percentiles), previous TB.
# SEs: model-based, CR1 and CR2 (Bell-McCaffrey) clustered by cell block,
# CR2 tested on Bell-McCaffrey/Satterthwaite df.

library(splines)

INT_BREAKS <- c(WINDOW, VISITS)          # (30, M3], (M3, M6], ..., (M21, M24]

make_long <- function(d, cutoff = Inf) {
  # cutoff: calendar day of analysis (interim); follow-up truncated there
  fu_end <- pmin(d$cens, cutoff - d$enrol)
  ev_t   <- ifelse(d$diag <= fu_end, d$diag, Inf)
  exit   <- pmin(fu_end, ev_t)
  keep   <- exit > WINDOW                 # mITT: at-risk time after day 30
  d <- d[keep, ]; exit <- exit[keep]; ev_t <- ev_t[keep]
  K <- length(INT_BREAKS) - 1
  i <- rep(seq_len(nrow(d)), each = K); k <- rep(seq_len(K), nrow(d))
  st <- INT_BREAKS[k]; en <- INT_BREAKS[k + 1]
  ok <- st < exit[i]
  i <- i[ok]; k <- k[ok]; st <- st[ok]; en <- en[ok]
  expo <- pmin(en, exit[i]) - st
  ev   <- as.integer(ev_t[i] > st & ev_t[i] <= en)
  data.frame(id = d$id[i], ev, expo = expo / DAY_YR, period = factor((k + 1) %/% 2),
             grp = factor(d$grp[i], 1:3, c("PBO", "B", "BMMS")),
             prison = factor(d$prison[i]), block = d$block[i],
             bmi = d$bmi[i], lunit = d$lunit[i], prevtb = d$prevtb[i])
}

fit_primary <- function(L, knots_bmi, knots_lun, cr2 = TRUE) {
  f <- ev ~ grp + prison + period +
       ns(bmi, knots = knots_bmi[2], Boundary.knots = knots_bmi[c(1, 3)]) +
       ns(lunit, knots = knots_lun[2], Boundary.knots = knots_lun[c(1, 3)]) + prevtb
  fit <- suppressWarnings(glm(f, family = binomial("cloglog"), data = L,
                              offset = log(L$expo), control = glm.control(maxit = 50)))
  X <- model.matrix(fit); mu <- fitted(fit); eta <- fit$linear.predictors
  mu.eta <- binomial("cloglog")$mu.eta(eta)
  w <- mu.eta^2 / (mu * (1 - mu))
  e <- (L$ev - mu) / mu.eta                     # working residuals
  B <- X * sqrt(w); u <- e * sqrt(w)
  G <- chol2inv(chol(crossprod(B)))
  cl <- split(seq_len(nrow(L)), L$block)
  nG <- length(cl); n <- nrow(L); p <- ncol(X)
  idx <- match(c("grpB", "grpBMMS"), colnames(X))
  est <- coef(fit)[idx]
  # CR0/CR1 meat
  S <- t(vapply(cl, function(ix) colSums(B[ix, , drop = FALSE] * u[ix]), numeric(p)))
  V0 <- G %*% crossprod(S) %*% G
  V1 <- V0 * nG / (nG - 1) * (n - 1) / (n - p)
  out <- list(est = est, se_mod = sqrt(diag(vcov(fit))[idx]),
              se_cr1 = sqrt(diag(V1)[idx]), conv = fit$converged, nG = nG,
              ev = tapply(L$ev, L$grp, sum))
  if (cr2) {
    # CR2: A_g = (I - H_gg)^(-1/2), H_gg = B_g G B_g' (rank <= p); computed via
    # the eigen-decomposition of the p x p matrix G^(1/2) B_g'B_g G^(1/2).
    eg <- eigen(G, symmetric = TRUE); Gh <- eg$vectors %*% (sqrt(eg$values) * t(eg$vectors))
    S2 <- matrix(0, nG, p); Pl <- vector("list", nG)
    for (j in seq_len(nG)) {
      ix <- cl[[j]]; Bg <- B[ix, , drop = FALSE]
      M <- Gh %*% crossprod(Bg) %*% Gh
      em <- eigen(M, symmetric = TRUE)
      lam <- pmin(em$values, 1 - 1e-10); keepv <- lam > 1e-12
      # A_g v = v + Bg Gh U diag(((1-l)^-1/2 - 1)/l) U' Gh Bg' v
      U <- em$vectors[, keepv, drop = FALSE]; l <- lam[keepv]
      Dm <- U %*% (((1 - l)^-0.5 - 1) / l * t(U))
      Ag_u <- u[ix] + Bg %*% (Gh %*% (Dm %*% (Gh %*% crossprod(Bg, u[ix]))))
      S2[j, ] <- crossprod(Bg, Ag_u)
      Pl[[j]] <- list(ix = ix, Dm = Dm)
    }
    V2 <- G %*% crossprod(S2) %*% G
    out$se_cr2 <- sqrt(diag(V2)[idx])
    # Bell-McCaffrey df for each coefficient (working-independence target)
    out$df_cr2 <- sapply(idx, function(jj) {
      cvec <- G[, jj]
      pg <- lapply(seq_len(nG), function(j) {
        ix <- Pl[[j]]$ix; Bg <- B[ix, , drop = FALSE]
        v <- Bg %*% cvec
        v + Bg %*% (Gh %*% (Pl[[j]]$Dm %*% (Gh %*% crossprod(Bg, v))))   # A_g symmetric
      })
      d1 <- vapply(seq_len(nG), function(j) sum(pg[[j]]^2), 0)
      Q  <- vapply(seq_len(nG), function(j) crossprod(B[cl[[j]], , drop = FALSE], pg[[j]])[, 1], numeric(p))
      Mx <- diag(d1) - t(Q) %*% G %*% Q
      sum(diag(Mx))^2 / sum(Mx^2)
    })
  }
  out
}

pooled_knots <- function(d) list(bmi = quantile(d$bmi, c(.10, .50, .90)),
                                 lun = quantile(d$lunit, c(.10, .50, .90)))

# interim (protocol 4.4): 92 primary events in B+MMS + placebo, or 12 months
# after last randomization; conditional power via B-value approximation
interim <- function(d, kn, target = 184, trig = 92) {
  obs <- d$diag <= d$cens & d$diag > WINDOW & d$grp %in% c(1, 3)
  cal <- sort(d$enrol[obs] + d$diag[obs])
  t_ev <- if (length(cal) >= trig) cal[trig] else Inf
  t_ia <- min(t_ev, max(d$enrol) + DAY_YR)
  L <- make_long(d, cutoff = t_ia)
  r <- fit_primary(L, kn$bmi, kn$lun, cr2 = FALSE)
  z <- -r$est[2] / r$se_cr1[2]
  dia <- sum(L$ev[L$grp %in% c("PBO", "BMMS")])
  tf <- min(dia / target, 0.99)
  b <- z * sqrt(tf)
  cp <- function(theta) 1 - pnorm((qnorm(0.975) - b - theta * (1 - tf)) / sqrt(1 - tf))
  cp_obs <- cp(z / sqrt(tf)); cp_des <- cp(qnorm(0.975) + qnorm(0.8))
  c(ia_month = t_ia / (DAY_YR / 12), ia_events = dia, ia_hr = exp(r$est[2]),
    cp_obs = cp_obs, cp_des = cp_des, futile = as.numeric(cp_obs < 0.15 & cp_des < 0.15))
}

analyze_trial <- function(d, do_interim = TRUE) {
  kn <- pooled_knots(d)
  L <- make_long(d)
  r <- fit_primary(L, kn$bmi, kn$lun, cr2 = TRUE)
  zc <- qnorm(0.975)
  p_mod <- 2 * pnorm(-abs(r$est / r$se_mod))
  p_cr1 <- 2 * pnorm(-abs(r$est / r$se_cr1))
  p_cr2 <- 2 * pt(-abs(r$est / r$se_cr2), df = r$df_cr2)
  res <- c(hr_b = exp(r$est[[1]]), hr_bmms = exp(r$est[[2]]),
           lhr_b = r$est[[1]], lhr_bmms = r$est[[2]],
           se_mod_b = r$se_mod[[1]], se_mod_bmms = r$se_mod[[2]],
           se_cr1_b = r$se_cr1[[1]], se_cr1_bmms = r$se_cr1[[2]],
           se_cr2_b = r$se_cr2[[1]], se_cr2_bmms = r$se_cr2[[2]],
           df_cr2_bmms = r$df_cr2[[2]], df_cr2_b = r$df_cr2[[1]],
           p_cr1_b = p_cr1[[1]], p_cr1_bmms = p_cr1[[2]],
           p_cr2_b = p_cr2[[1]], p_cr2_bmms = p_cr2[[2]],
           p_mod_b = p_mod[[1]], p_mod_bmms = p_mod[[2]],
           ev_pbo = r$ev[["PBO"]], ev_b = r$ev[["B"]], ev_bmms = r$ev[["BMMS"]],
           conv = as.numeric(r$conv), n_blocks = r$nG, n_mitt = length(unique(L$id)))
  if (do_interim) res <- c(res, interim(d, kn))
  res
}
