# Scenario definitions (SAP Appendix A). Each entry: changes to default_par() and reps.
sc <- function(reps, ...) list(reps = reps, mod = list(...))
H <- 0.656
SCENARIOS <- list(
  # type I error
  N1          = sc(3000, hr_bmms = 1, hr_b = 1),
  N1_sd0      = sc(2000, hr_bmms = 1, hr_b = 1, block_sd = 0),
  N1_sd05     = sc(2000, hr_bmms = 1, hr_b = 1, block_sd = 0.5),
  N2          = sc(2000, hr_bmms = 1, hr_b = H),
  N3          = sc(3000, hr_bmms = H, hr_b = 1),
  # power at design
  P1          = sc(2000, hr_bmms = H, hr_b = H),
  P2          = sc(1500, hr_bmms = 0.65, hr_b = 0.65),
  P3          = sc(1500, hr_bmms = 0.65, hr_b = 0.70),
  P4          = sc(1000, hr_bmms = 0.70, hr_b = 0.75),
  # cell blocks of ~30 (285 blocks)
  N1_b30      = sc(3000, hr_bmms = 1, hr_b = 1, block_size = 30),
  N1_sd05_b30 = sc(2000, hr_bmms = 1, hr_b = 1, block_sd = 0.5, block_size = 30),
  P1_b30      = sc(2000, hr_bmms = H, hr_b = H, block_size = 30),
  # robustness (from P1)
  R_risk2     = sc(1000, risk = 0.02),
  R_risk25    = sc(1000, risk = 0.025),
  R_risk4     = sc(1000, risk = 0.04),
  R_state     = sc(1000, state_risk = c(MS = 1.15, RS = 0.67, RO = 0.67)),
  R_loss20    = sc(1000, loss = 0.20),
  R_loss45    = sc(1000, loss = 0.45),
  R_lossfront = sc(1000, loss_shape = "front"),
  R_lossrisk  = sc(1000, loss_risk = 0.3),
  R_lossdiff  = sc(1000, loss_diff = 0.03),
  R_sd0       = sc(1000, block_sd = 0),
  R_sd02      = sc(1000, block_sd = 0.2),
  R_sd05      = sc(1000, block_sd = 0.5),
  R_wane      = sc(1000, effect = "wane"),
  R_half      = sc(1000, effect = "half_after_m6"),
  R_delay     = sc(1000, effect = "delayed_m3"),
  R_lag30     = sc(1500, window_eff = 0),       # no effect on onset in days 0-30
  R_nonadh    = sc(1000, nonadh = 0.25),
  R_contam5   = sc(1000, contam = 0.05),
  R_contam10  = sc(1000, contam = 0.10),
  R_leak      = sc(1000, prev_leak = 0.2),
  R_visit60   = sc(1000, symp_mean = 70),
  R_visit95   = sc(1000, symp_mean = 1500),
  R_slow      = sc(1000, enrol_delay = 3)
)
for (nm in names(SCENARIOS)[grepl("^R_", names(SCENARIOS))])
  SCENARIOS[[nm]]$mod <- modifyList(list(hr_bmms = H, hr_b = H), SCENARIOS[[nm]]$mod)

