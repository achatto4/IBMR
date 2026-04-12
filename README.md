# IBMR

`IBMR` provides tools for instrument borrowing in Mendelian randomization when you have:

- one exposure of interest
- one primary outcome trait
- one auxiliary outcome trait, or several candidate auxiliary traits

The package is built around a simple workflow:

1. Start with SNP-level summary statistics for the exposure and outcome traits.
2. Measure coheterogeneity between the primary outcome and one or more auxiliary traits.
3. If several auxiliaries are available, pick the auxiliary trait with the strongest coheterogeneity with the primary outcome.
4. Re-use that auxiliary trait in downstream instrument-borrowing analyses with `IBMODE()` or `IBPRESSO()`.

## What Problem This Package Solves

In many Mendelian randomization settings, a primary outcome may share pleiotropic structure with another trait. `IBMR` helps you quantify that shared structure and then exploit it in downstream estimation and outlier detection.

The package currently exposes three main functions:

- `coheterogeneity_Q()`: estimates pairwise coheterogeneity correlation across traits using a guarded theoretical delta-exact method
- `IBMODE()`: performs multidimensional mode-based estimation across outcomes
- `IBPRESSO()`: performs MR-PRESSO with an auxiliary trait for instrument borrowing

## Installation

You can install the package locally from the repository root with:

```r
install.packages(c("MASS", "ks"))

# from the package directory
devtools::install(".")
```

Or load it directly during development:

```r
library(IBMR)
```

## Required Inputs

At minimum, the package expects SNP-level summary statistics for:

- the exposure
- the primary outcome
- one auxiliary trait

If you want to compare several candidate auxiliaries, you provide the same exposure summary statistics together with the primary outcome and all candidate auxiliary outcomes.

### Core objects

For most analyses you will prepare:

- `BetaXG`: numeric vector of SNP-exposure effects
- `seBetaXG`: numeric vector of SNP-exposure standard errors
- `BetaYG_matrix`: matrix of SNP-outcome effects
- `seBetaYG_matrix`: matrix of SNP-outcome standard errors

The rows of all objects must refer to the same SNPs in the same order.

### Recommended matrix layout

For coheterogeneity screening, put the primary outcome in one column and each candidate auxiliary trait in the remaining columns.

For example:

```r
colnames(BetaYG_matrix)
# [1] "primary_trait" "aux_trait_1" "aux_trait_2" "aux_trait_3"
```

The same column order should be used in `seBetaYG_matrix`.

## Typical Workflow

### Scenario 1: One primary outcome and one auxiliary trait

This is the simplest use case. You already know which auxiliary trait you want to borrow from.

1. Build a two-column outcome matrix containing the primary outcome and the auxiliary trait.
2. Run `coheterogeneity_Q()`.
3. Read the off-diagonal value of `rho` as the coheterogeneity estimate between the two traits.
4. Use the same primary/auxiliary pair in `IBMODE()` or `IBPRESSO()`.

Example:

```r
library(IBMR)

# Example objects:
# BetaXG           : length K
# seBetaXG         : length K
# beta_primary     : length K
# se_primary       : length K
# beta_aux         : length K
# se_aux           : length K

BetaYG_matrix <- cbind(
  primary_trait = beta_primary,
  auxiliary_trait = beta_aux
)

seBetaYG_matrix <- cbind(
  primary_trait = se_primary,
  auxiliary_trait = se_aux
)

cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix
)

cohet_res$rho["primary_trait", "auxiliary_trait"]
cohet_res$p_value["primary_trait", "auxiliary_trait"]
cohet_res$flag["primary_trait", "auxiliary_trait"]
```

Interpretation:

- `rho` is the coheterogeneity correlation estimate
- `p_value` gives the pairwise significance level
- `flag` tells you whether the pair was estimable or whether a guardrail was triggered

### Scenario 2: One primary outcome and multiple candidate auxiliary traits

This is the main screening use case.

1. Put the primary outcome and all candidate auxiliaries into the outcome matrices.
2. Run `coheterogeneity_Q()`.
3. Extract the row or column corresponding to the primary outcome.
4. Rank candidate auxiliaries by the magnitude of coheterogeneity, typically `abs(rho)`.
5. Optionally prioritize significant pairs using `p_value`.
6. Choose the top auxiliary trait and carry it into `IBMODE()` or `IBPRESSO()`.

Example:

```r
library(IBMR)

BetaYG_matrix <- cbind(
  primary_trait = beta_primary,
  aux_trait_1 = beta_aux1,
  aux_trait_2 = beta_aux2,
  aux_trait_3 = beta_aux3
)

seBetaYG_matrix <- cbind(
  primary_trait = se_primary,
  aux_trait_1 = se_aux1,
  aux_trait_2 = se_aux2,
  aux_trait_3 = se_aux3
)

cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix
)

rho_row <- cohet_res$rho["primary_trait", ]
p_row <- cohet_res$p_value["primary_trait", ]

rho_row["primary_trait"] <- NA
p_row["primary_trait"] <- NA

ranking <- data.frame(
  aux_trait = names(rho_row),
  rho = unname(rho_row),
  abs_rho = abs(unname(rho_row)),
  p_value = unname(p_row),
  flag = unname(cohet_res$flag["primary_trait", ])
)

ranking <- ranking[!is.na(ranking$rho), ]
ranking <- ranking[order(-ranking$abs_rho), ]
ranking
```

A common decision rule is:

- first prefer auxiliaries with `flag == "OK"`
- then prioritize pairs with small `p_value`
- among those, choose the largest `abs(rho)`

## Understanding `coheterogeneity_Q()`

`coheterogeneity_Q()` now uses a guarded theoretical method with several built-in protections:

- optional SNP restriction through `SNP_keep`
- minimum exposure effect filtering via `bx_min`
- weak-instrument filtering via `F_min`
- optional winsorization of ratio estimates via `winsor_theta`
- optional LDSC intercept adjustment via `ldsc_intercepts`
- pair-level flags when a comparison is not estimable

### Main arguments

```r
coheterogeneity_Q(
  BetaXG,
  BetaYG_matrix = NULL,
  seBetaXG,
  seBetaYG_matrix = NULL,
  BetaIV_matrix = NULL,
  seBetaIV_matrix = NULL,
  ldsc_intercepts = NULL,
  SNP_keep = NULL,
  use_ldsc = TRUE,
  alpha = 0.05,
  step_mult = 1e-3,
  eps = 1e-12,
  bx_min = 1e-6,
  F_min = 0,
  winsor_theta = 20,
  min_K_pair = 50,
  return_diagnostics = FALSE
)
```

### Returned object

The function returns a list with:

- `rho`: coheterogeneity correlation matrix
- `se`: standard errors for `rho`
- `z_statistic`: z statistics for each trait pair
- `p_value`: p values for each trait pair
- `K`: number of SNPs used for each pair
- `flag`: pair-specific status labels such as `"OK"`, `"tau0"`, or `"se NA"`
- `guards`: the guardrail settings used in the run

### Important note

The current `coheterogeneity_Q()` output is no longer the old `Q_matrix` / `Q_corr_matrix` interface. Downstream code should use `rho`, `p_value`, and `flag`.

## Using the Selected Auxiliary in `IBMODE()`

After screening, take the primary outcome together with the chosen auxiliary trait and fit a two-outcome multidimensional mode-based estimator.

Example:

```r
library(IBMR)

BetaYG_mode <- cbind(
  primary_trait = beta_primary,
  chosen_aux = beta_aux2
)

seBetaYG_mode <- cbind(
  primary_trait = se_primary,
  chosen_aux = se_aux2
)

ibmode_res <- IBMODE(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_mode,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_mode,
  phi = c(1, 0.5, 0.25),
  n_boot = 1000
)

ibmode_res
```

Use `IBMODE()` when you want joint mode-based estimation across the primary and auxiliary traits.

## Using the Selected Auxiliary in `IBPRESSO()`

`IBPRESSO()` expects a data frame and column names rather than matrices. Once you have selected the auxiliary trait, build a SNP-level table containing:

- exposure effect and exposure SE
- primary outcome effect and primary outcome SE
- auxiliary trait effect and auxiliary trait SE

Example:

```r
library(IBMR)

dat_ibpresso <- data.frame(
  beta_exposure = BetaXG,
  se_exposure = seBetaXG,
  beta_primary = beta_primary,
  se_primary = se_primary,
  beta_aux = beta_aux2,
  se_aux = se_aux2
)

ibpresso_res <- IBPRESSO(
  BetaOutcome = "beta_primary",
  BetaExposure = "beta_exposure",
  BetaAux = "beta_aux",
  SdOutcome = "se_primary",
  SdExposure = "se_exposure",
  SdAux = "se_aux",
  data = dat_ibpresso,
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  NbDistribution = 1000,
  seed = 123
)

ibpresso_res
```

Use `IBPRESSO()` when you want to detect pleiotropic outliers and assess whether removing them changes the primary MR estimate.

## End-to-End Example

This pattern is usually what users will do in practice:

```r
library(IBMR)

# Step 1: screen auxiliaries
cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  F_min = 20,
  winsor_theta = 20
)

# Step 2: rank candidate auxiliaries for the primary outcome
primary_name <- "primary_trait"
rho_primary <- cohet_res$rho[primary_name, ]
rho_primary[primary_name] <- NA

best_aux <- names(which.max(abs(rho_primary)))
best_aux

# Step 3: use the selected auxiliary in IBMODE
ibmode_res <- IBMODE(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix[, c(primary_name, best_aux), drop = FALSE],
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix[, c(primary_name, best_aux), drop = FALSE]
)
```

If you prefer a more conservative screening rule, rank only auxiliaries with:

- `flag == "OK"`
- finite `p_value`
- `p_value < 0.05`

and then choose the largest `abs(rho)` among them.

## Input Quality Checks

Before running the package, it is a good idea to confirm that:

- SNP order matches across all vectors and matrices
- effect alleles are harmonized across exposure, primary outcome, and auxiliary outcomes
- standard errors are positive and finite
- very weak instruments have been removed or filtered using `F_min`
- missing data have been handled consistently

## Practical Notes

- If you only have one auxiliary trait, you do not need a ranking step.
- If you have many candidate auxiliaries, `coheterogeneity_Q()` can be used as a screening layer before downstream modeling.
- A large absolute coheterogeneity estimate does not automatically imply significance; check `p_value` and `flag`.
- If a pair returns `NA`, inspect `flag` and `K` to understand why that pair was not estimable.
- `IBMODE()` uses bootstrap resampling and can be slow when `n_boot` is large.
- `IBPRESSO()` can also be computationally expensive when `NbDistribution` is large.

## Minimal Reproducible Template

```r
library(IBMR)

# Exposure
BetaXG <- ...
seBetaXG <- ...

# Primary outcome
beta_primary <- ...
se_primary <- ...

# Candidate auxiliaries
beta_aux1 <- ...
se_aux1 <- ...
beta_aux2 <- ...
se_aux2 <- ...

BetaYG_matrix <- cbind(
  primary_trait = beta_primary,
  aux_trait_1 = beta_aux1,
  aux_trait_2 = beta_aux2
)

seBetaYG_matrix <- cbind(
  primary_trait = se_primary,
  aux_trait_1 = se_aux1,
  aux_trait_2 = se_aux2
)

cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix
)

cohet_res$rho
cohet_res$p_value
cohet_res$flag
```

## Citation

If you use this package in applied work, please cite the repository and the relevant method paper(s) associated with your analysis.
