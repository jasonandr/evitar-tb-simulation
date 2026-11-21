# Simulation inputs

`cohort_covariates.rds` holds aggregate summaries of the TB R01 prison cohort: 3,414 men in Mato Grosso do Sul with negative baseline sputum microbiology, followed for 5,850 person-years with 241 microbiologically confirmed TB events. These are the data in the cohort TB-incidence manuscript. The file contains no individual-level records, only:

- the proportion with previous TB (0.186);
- the 0.5%-step quantiles of the baseline Lunit INSIGHT CXR (v3.1.5.1) score, separately for participants with and without previous TB;
- the coefficients of a Cox model for incident TB with a natural cubic spline of the Lunit score, its interaction with previous TB, and previous TB itself. The spline has knots at 6.65, 27.68, 49.48 and 92.63, the same knots as the manuscript.

The summaries were produced from the cohort's analysis dataset, which is not public. They are used when `cov_model = "cohort"`, the default in `R/dgm.R`. BMI uses an HR of 1.21 per kg/m² decrease, the cohort estimate cited in the protocol.

The investigators have approved public release of these aggregate summaries.
