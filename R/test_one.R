source("dgm.R"); source("analysis.R")
par <- default_par()
t0 <- Sys.time(); h0 <- calibrate_h0(par); cat("h0", h0, "calib s", format(Sys.time() - t0), "\n")
set.seed(11)
t0 <- Sys.time(); d <- gen_trial(par, h0); cat("gen s", format(Sys.time() - t0), "\n")
cat("alloc", table(d$grp), "\n")
cat("blocks", length(unique(d$block)), " last enrol month", max(d$enrol) / 30.44, "\n")
t0 <- Sys.time(); r <- analyze_trial(d); cat("analyze s", format(Sys.time() - t0), "\n")
print(round(r, 4))
# validate CR2 against clubSandwich
kn <- pooled_knots(d); L <- make_long(d)
f <- ev ~ grp + prison + period +
  ns(bmi, knots = kn$bmi[2], Boundary.knots = kn$bmi[c(1, 3)]) +
  ns(lunit, knots = kn$lun[2], Boundary.knots = kn$lun[c(1, 3)]) + prevtb
fit <- glm(f, family = binomial("cloglog"), data = L, offset = log(L$expo))
t0 <- Sys.time()
ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = L$block, test = "Satterthwaite",
                              coefs = c("grpB", "grpBMMS"))
cat("clubSandwich s", format(Sys.time() - t0), "\n"); print(ct)
cat("CR1 sandwich:", sqrt(diag(sandwich::vcovCL(fit, cluster = L$block, type = "HC1")))[2:3], "\n")
# diagnosis timing: share at scheduled visits among counted events
ev <- d$diag <= d$cens & d$diag > WINDOW
cat("share of counted events on a visit day:", mean(d$diag[ev] %in% VISITS), "\n")
cat("prevalent window exclusions:", sum(d$diag <= WINDOW), "\n")
