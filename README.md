# EVITAR-TB: simulation of the operating characteristics of the primary analysis

Trial-level simulation supporting the statistical analysis plan (SAP) of EVITAR-TB, a three-arm, individually randomized, double-blind, placebo-controlled trial of B vitamins (B1, B2, B3, B6), alone or with multiple micronutrients (MMS), to prevent pulmonary tuberculosis in adult men in eight Brazilian prisons (N = 8,650).

The simulation follows the ADEMP framework (Morris, White and Crowther, *Stat Med* 2019). Each simulated trial reproduces:

- enrollment in the eight prisons;
- cell blocks within prisons, with a log-normal block frailty;
- Pocock–Simon minimization within prison (protocol Section 6.4);
- TB onset, detected at quarterly visits or by symptomatic presentation;
- the 30-day window period;
- censoring at exit from the prison system.

Each trial is analyzed with the prespecified primary analysis (protocol Section 10.4):

- **Model:** a discrete-time cloglog model with a log-exposure offset and 6-month baseline-hazard periods. It adjusts for prison, BMI and Lunit score (restricted cubic splines), and previous TB.
- **Variance:** CR2 cluster-robust variance by cell block, tested on Satterthwaite (Bell–McCaffrey) degrees of freedom.
- **Testing:** the fixed sequence (H1 B+MMS vs placebo, then H2 B vs placebo).
- **Interim:** the event-driven interim with its non-binding futility rule.

Baseline TB risk follows the TB R01 prison cohort (Lunit score and previous TB drawn from the cohort, with the cohort's spline-by-previous-TB risk model; BMI HR 1.21 per kg/m² decrease). Several inputs remain placeholders until site data are available: the cell-block census and between-block variation, site activation dates, and the release pattern. They are flagged `PLACEHOLDER` in `R/dgm.R`.

## Files

| File | Contents |
|---|---|
| `R/dgm.R` | Data-generating mechanism: sites, blocks, covariates, minimization, onset, detection, censoring |
| `R/analysis.R` | Primary analysis, CR1/CR2 variance and degrees of freedom, interim conditional power |
| `R/scenarios.R` | Scenarios (null, hypothesized effect, robustness; `_b30`: cell blocks of ~30) and number of replicates |
| `R/run.R` | Runs scenarios in parallel; writes `results/<scenario>.rds` |
| `R/summarize.R` | Performance measures; writes `results/summary.csv` |
| `R/test_one.R` | Single-trial check, including validation of the CR2 implementation against `clubSandwich` |
| `sample_size/protocol_tables.py` | Closed-form sample size and power calculations in the protocol and SAP, with the 30-day window period (Python 3 with NumPy and SciPy) |
| `inputs/cohort_covariates.rds` | Aggregate summaries of the TB R01 prison cohort used by all scenarios (see inputs/README.md) |
| `results/summary.csv` | Summary of all scenarios (per-replicate `.rds` files are not tracked) |

## Running

Sample size tables: `python3 sample_size/protocol_tables.py`.

Simulation:

Requires R (developed with 4.5.1). The simulation itself uses only base R plus `splines` and `parallel`. `test_one.R` also uses `clubSandwich` and `sandwich`. Package versions are recorded in `renv.lock`.

```sh
cd R
Rscript run.R            # all scenarios; or: Rscript run.R P1 N1
Rscript summarize.R
```

`run.R` uses the L'Ecuyer-CMRG generator with a seed derived from each scenario name. It uses all cores but one and skips scenarios whose results file already exists. The full set takes about 3.5 hours on 9 cores.

## Author and license

Jason Andrews, Stanford University. Code released under the MIT License (see `LICENSE`).
