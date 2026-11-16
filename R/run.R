# Usage: Rscript run.R [scenario ...]   (default: all). Results in ../results/<name>.rds
source("dgm.R"); source("analysis.R"); source("scenarios.R")
library(parallel)
args <- commandArgs(TRUE)
todo <- if (length(args)) args else names(SCENARIOS)
cores <- max(1, detectCores() - 1)
OUT <- Sys.getenv("RESULTS_DIR", "../results")
dir.create(OUT, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")
for (nm in todo) {
  f <- file.path(OUT, paste0(nm, ".rds"))
  if (file.exists(f)) next
  s <- SCENARIOS[[nm]]
  par <- modifyList(default_par(), s$mod)
  h0 <- calibrate_h0(par)
  set.seed(sum(utf8ToInt(nm)) * 1000 + 20261112)
  t0 <- Sys.time()
  out <- mclapply(seq_len(s$reps), function(i) {
    r <- try(analyze_trial(gen_trial(par, h0)), silent = TRUE)
    if (inherits(r, "try-error")) NULL else r
  }, mc.cores = cores, mc.set.seed = TRUE)
  ok <- !vapply(out, is.null, TRUE)
  res <- do.call(rbind, out[ok])
  saveRDS(list(name = nm, par = par, h0 = h0, res = res, failed = sum(!ok)), f)
  cat(sprintf("%s: %d reps (%d failed) in %.1f min\n", nm, sum(ok), sum(!ok),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
