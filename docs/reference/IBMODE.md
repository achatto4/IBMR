# Multidimensional Mode-Based Estimation for Mendelian Randomization

Performs mode-based estimation (MBE) for multiple correlated outcomes
simultaneously using multidimensional kernel density estimation. This
approach identifies the modal causal effect across multiple traits while
accounting for pleiotropy.

## Usage

``` r
IBMODE(
  BetaXG,
  BetaYG_matrix,
  seBetaXG,
  seBetaYG_matrix,
  phi = c(1, 0.5, 0.25),
  n_boot = 10000,
  alpha = 0.05,
  cov_ratio = NULL,
  ldsc_intercept = NULL,
  seed = NULL
)
```

## Arguments

- BetaXG:

  Numeric vector of SNP-exposure associations

- BetaYG_matrix:

  Matrix of SNP-outcome associations (rows = SNPs, columns =
  outcomes/traits)

- seBetaXG:

  Numeric vector of standard errors for BetaXG

- seBetaYG_matrix:

  Matrix of standard errors for BetaYG_matrix

- phi:

  Numeric vector of bandwidth multipliers for sensitivity analysis
  (default: c(1, 0.5, 0.25))

- n_boot:

  Integer, number of bootstrap iterations for standard error estimation
  (default: 10000)

- alpha:

  Numeric, significance level for confidence intervals (default: 0.05)

- cov_ratio:

  Optional numeric vector of per-instrument cross-trait covariances of
  the ratio estimates (sigma_12,k); when supplied (or derived from
  `ldsc_intercept`), the bootstrap draws each instrument's outcome pair
  jointly from a bivariate normal, accounting for correlation induced by
  overlapping outcome GWAS samples. `NULL` (default) or all-zero gives
  the independent bootstrap. Two-outcome analyses only.

- ldsc_intercept:

  Optional scalar cross-trait LD-score regression intercept for the
  (primary, auxiliary) pair; if `cov_ratio` is `NULL` and this is
  supplied, `cov_ratio` is formed internally.

- seed:

  Optional integer random seed for a reproducible bootstrap (default:
  `NULL`).

## Value

A data frame containing:

- Method:

  Method name

- phi:

  Bandwidth multiplier used

- Estimate.1, Estimate.2, ...:

  Causal effect estimates for each outcome

- SE.1, SE.2, ...:

  Standard errors for each outcome

- CI_low.1, CI_low.2, ...:

  Lower confidence interval bounds

- CI_upp.1, CI_upp.2, ...:

  Upper confidence interval bounds

- P.1, P.2, ...:

  P-values for each outcome

## Details

The function uses weighted multidimensional kernel density estimation to
find the mode of the joint distribution of causal effects across
multiple outcomes. Bootstrap resampling is used to estimate standard
errors; when the primary and auxiliary outcome GWAS share samples,
supply `cov_ratio` (or `ldsc_intercept`) so the bootstrap is drawn from
the appropriate bivariate distribution. Multiple bandwidth parameters
(phi) can be tested for sensitivity analysis.
