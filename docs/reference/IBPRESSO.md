# MR-PRESSO with Instrument Borrowing

Performs Mendelian Randomization Pleiotropy RESidual Sum and Outlier
(MR-PRESSO) test extended for instrument borrowing scenarios. Detects
and corrects for horizontal pleiotropy using an auxiliary trait to
identify outlier instruments.

## Usage

``` r
IBPRESSO(
  BetaOutcome,
  BetaExposure,
  BetaAux,
  SdOutcome,
  SdExposure,
  SdAux,
  data,
  OUTLIERtest = FALSE,
  DISTORTIONtest = FALSE,
  SignifThreshold = 0.05,
  NbDistribution = 1000,
  seed = NULL,
  n_cores = 1
)
```

## Arguments

- BetaOutcome:

  Character, column name for primary outcome SNP effects

- BetaExposure:

  Character or vector, column name(s) for exposure SNP effects

- BetaAux:

  Character, column name for auxiliary trait SNP effects (borrowed
  instrument)

- SdOutcome:

  Character, column name for standard errors of outcome effects

- SdExposure:

  Character or vector, column name(s) for standard errors of exposure
  effects

- SdAux:

  Character, column name for standard errors of auxiliary trait effects

- data:

  Data frame containing all the above columns

- OUTLIERtest:

  Logical, whether to perform outlier detection (default: FALSE)

- DISTORTIONtest:

  Logical, whether to test if outliers cause distortion (default: FALSE)

- SignifThreshold:

  Numeric, significance threshold for tests (default: 0.05)

- NbDistribution:

  Integer, number of null distributions to generate (default: 1000)

- seed:

  Integer, random seed for reproducibility (default: NULL)

- n_cores:

  Integer, number of cores for parallel processing (default: 1)

## Value

A list containing:

- raw_beta:

  Raw causal effect estimate before outlier correction

- raw_se:

  Standard error of raw estimate

- corrected_beta:

  Outlier-corrected causal effect estimate (if outliers detected)

- corrected_se:

  Standard error of corrected estimate

- p_value:

  Global test p-value for horizontal pleiotropy

- p_distort:

  P-value for distortion test (if DISTORTIONtest = TRUE)

- outlier_idx:

  Indices of detected outlier instruments

- n_instruments:

  Number of instruments used

- n_outliers:

  Number of outliers detected

## Details

This function extends MR-PRESSO to handle instrument borrowing by
incorporating an auxiliary trait. It uses leave-one-out cross-validation
and Mahalanobis distance in a bivariate residual space (outcome +
auxiliary trait) to detect pleiotropic outliers. The global test
assesses overall horizontal pleiotropy, while the distortion test
evaluates whether removing outliers significantly changes the causal
estimate.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
n_snps <- 50
dat <- data.frame(
  beta_exposure = rnorm(n_snps, 0.1, 0.02),
  beta_outcome = rnorm(n_snps, 0.05, 0.03),
  beta_aux = rnorm(n_snps, 0.04, 0.025),
  se_exposure = abs(rnorm(n_snps, 0.01, 0.002)),
  se_outcome = abs(rnorm(n_snps, 0.015, 0.003)),
  se_aux = abs(rnorm(n_snps, 0.012, 0.002))
)

# Run MR-PRESSO with instrument borrowing
results <- mrpresso_ib(
  BetaOutcome = "beta_outcome",
  BetaExposure = "beta_exposure",
  BetaAux = "beta_aux",
  SdOutcome = "se_outcome",
  SdExposure = "se_exposure",
  SdAux = "se_aux",
  data = dat,
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  NbDistribution = 1000
)
} # }
```
