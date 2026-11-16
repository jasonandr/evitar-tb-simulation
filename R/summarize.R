# Summarize ../results/*.rds into ../results/summary.csv and print a table.
OUT <- Sys.getenv("RESULTS_DIR", "../results")
fs <- list.files(OUT, "\\.rds$", full.names = TRUE)
mc <- function(p, n) sqrt(p * (1 - p) / n)
rows <- lapply(fs, function(f) {
  o <- readRDS(f); r <- as.data.frame(o$res); n <- nrow(r)
  rej1_cr2 <- r$p_cr2_bmms < 0.05 & r$hr_bmms < 1
  rej2_cr2 <- rej1_cr2 & r$p_cr2_b < 0.05 & r$hr_b < 1
  rej1_cr1 <- r$p_cr1_bmms < 0.05 & r$hr_bmms < 1
  rej2_cr1 <- rej1_cr1 & r$p_cr1_b < 0.05 & r$hr_b < 1
  t1 <- log(o$par$hr_bmms)
  fut <- r$futile == 1
  data.frame(
    scenario = o$name, reps = n,
    hr_bmms = o$par$hr_bmms, hr_b = o$par$hr_b,
    ev_pbo = mean(r$ev_pbo), ev_bmms = mean(r$ev_bmms), ev_b = mean(r$ev_b),
    power_H1 = mean(rej1_cr2), mcse = mc(mean(rej1_cr2), n),
    power_H2 = mean(rej2_cr2),
    power_H1_cr1 = mean(rej1_cr1), power_H2_cr1 = mean(rej2_cr1),
    # two-sided rejection rates (type I error under nulls)
    rej2s_bmms_cr2 = mean(r$p_cr2_bmms < 0.05), rej2s_bmms_cr1 = mean(r$p_cr1_bmms < 0.05),
    rej2s_bmms_mod = mean(r$p_mod_bmms < 0.05),
    rej2s_b_cr2 = mean(r$p_cr2_b < 0.05),
    fwer_cr2 = mean(rej1_cr2 & o$par$hr_bmms == 1 | rej2_cr2 & o$par$hr_b == 1),
    bias_lhr = mean(r$lhr_bmms) - t1,
    emp_se = sd(r$lhr_bmms), mean_se_cr1 = mean(r$se_cr1_bmms),
    mean_se_cr2 = mean(r$se_cr2_bmms), mean_se_mod = mean(r$se_mod_bmms),
    cover_cr2 = mean(abs(r$lhr_bmms - t1) <= qt(.975, r$df_cr2_bmms) * r$se_cr2_bmms),
    df_cr2 = mean(r$df_cr2_bmms), blocks = mean(r$n_blocks),
    nonconv = mean(r$conv == 0),
    ia_month = median(r$ia_month), ia_events = mean(r$ia_events),
    p_futile = mean(fut),
    power_H1_if_rule = mean(rej1_cr2 & !fut))
})
s <- do.call(rbind, rows)
ord <- c(grep("^N", s$scenario), grep("^P", s$scenario), grep("^R_", s$scenario))
s <- s[ord, ]
write.csv(s, file.path(OUT, "summary.csv"), row.names = FALSE)
options(width = 200)
print(format(s[, c("scenario", "reps", "ev_pbo", "ev_bmms", "power_H1", "power_H2",
                   "rej2s_bmms_cr2", "rej2s_bmms_cr1", "rej2s_b_cr2", "fwer_cr2",
                   "bias_lhr", "emp_se", "mean_se_cr2", "cover_cr2", "df_cr2",
                   "ia_month", "p_futile", "power_H1_if_rule")], digits = 3), row.names = FALSE)
