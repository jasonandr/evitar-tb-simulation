# Large-sample check of which DGM features pull the estimated HR toward 1.
source("dgm.R"); source("analysis.R")
cfg <- expand.grid(window_eff = c(0, 1), block_sd = c(0, 0.35))
set.seed(7)
for (j in seq_len(nrow(cfg))) {
  par <- modifyList(default_par(), list(window_eff = cfg$window_eff[j], block_sd = cfg$block_sd[j],
                                        hr_b = 0.656, hr_bmms = 0.656))
  par$site_n <- par$site_n * 10; par$block_size <- 900   # same number of blocks, 10x size
  h0 <- calibrate_h0(par)
  est <- replicate(3, {
    d <- gen_trial(par, h0); kn <- pooled_knots(d)
    r <- fit_primary(make_long(d), kn$bmi, kn$lun, cr2 = FALSE)
    c(hr = exp(r$est[[2]]), pbo = r$ev[["PBO"]] / 10, bmms = r$ev[["BMMS"]] / 10)
  })
  cat(sprintf("window_eff=%g block_sd=%.2f: HR %.3f; events per 8,650: placebo %.1f, B+MMS %.1f\n",
              cfg$window_eff[j], cfg$block_sd[j], mean(est["hr", ]), mean(est["pbo", ]), mean(est["bmms", ])))
}
