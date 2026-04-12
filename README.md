# IBMR

`IBMR` implements an instrument borrowing framework for Mendelian randomization
(MR) using summary-level genetic association data.

The package is designed for settings in which the analyst has:

- one exposure of interest
- one primary outcome trait
- one auxiliary outcome trait, or several candidate auxiliary traits

The central idea is that closely related outcome traits may share overlapping
sets of valid instruments, or may reflect similar underlying mechanisms of
instrument invalidity, with respect to a common exposure. When such overlap is
present, borrowing information across traits can improve robustness and
efficiency in MR estimation.

The package is organized around the following workflow:

1. Start with SNP-level summary statistics for the exposure and outcome traits.
2. Measure coheterogeneity between the primary outcome and one or more auxiliary traits.
3. If several auxiliaries are available, identify the auxiliary trait with the strongest coheterogeneity with the primary outcome.
4. Re-use that auxiliary trait in downstream instrument-borrowing analyses with `IBMODE()` or `IBPRESSO()`.

## Overview

Mendelian randomization is widely used to estimate causal effects of exposures
on outcomes using genetic variants as instrumental variables. In practice,
however, standard MR procedures can lose power or become biased when a
substantial fraction of candidate instruments are invalid.

`IBMR` addresses this problem by introducing a coheterogeneity-based screening
step for identifying auxiliary traits that share relevant heterogeneity
structure with a primary outcome of interest. The selected auxiliary trait can
then be incorporated into downstream robust MR procedures through:

- `IBMODE()`, an instrument-borrowing extension of mode-based estimation
- `IBPRESSO()`, an instrument-borrowing extension of MR-PRESSO

This package is therefore intended for analyses in which the goal is not only
to estimate a causal effect for a primary outcome, but also to strengthen that
analysis by leveraging information from related outcome traits.

## Graphical Overview

![IBMR graphical summary](man/figures/ibmr-graphical-summary.png)

The figure summarizes the conceptual and methodological workflow implemented in
`IBMR`.

- Panel A shows the motivating setting: some candidate instruments for the
  exposure are valid, while others are invalid because of pleiotropic or
  confounded pathways. The primary outcome and an auxiliary outcome may still
  share useful structure that can be exploited in joint analysis.
- Panel B shows the auxiliary-trait screening step. For a given primary trait,
  `IBMR` calculates coheterogeneity with each candidate auxiliary trait and
  uses this information to identify the most suitable auxiliary trait for
  borrowing.
- Panel C shows the SNP-level summary-statistics inputs used by the methods:
  SNP-exposure effects, SNP-outcome effects, and ratio estimates for the
  primary and auxiliary traits.
- Panel D shows the downstream joint-analysis idea. Once an auxiliary trait has
  been selected, the paired traits can be analyzed jointly to improve robust
  estimation and outlier detection using `IBMODE()` and `IBPRESSO()`.

## What Problem This Package Solves

In many Mendelian randomization applications, a primary outcome may share
pleiotropic structure with one or more related traits. `IBMR` is designed to
quantify that shared structure and then exploit it in downstream estimation and
outlier detection.

The package currently exposes three main functions:

- `coheterogeneity_Q()`: estimates pairwise coheterogeneity across traits using
  a guarded theoretical delta-exact method
- `IBMODE()`: performs multidimensional mode-based estimation across outcomes
- `IBPRESSO()`: performs MR-PRESSO with an auxiliary trait for instrument
  borrowing

## Installation

You can install `IBMR` directly from GitHub.

### Install from GitHub

First install `devtools` from CRAN:

```r
install.packages("devtools")
```

Then load `devtools`:

```r
library(devtools)
```

Install `IBMR` from GitHub:

```r
devtools::install_github("achatto4/IBMR")
```

Finally, load the package:

```r
library(IBMR)
```

### Local installation during development

If you are working from a local clone of the repository, you can also install
the package from the package directory:

```r
install.packages(c("MASS", "ks"))
devtools::install(".")
library(IBMR)
```

## Required Inputs

At minimum, the package requires SNP-level summary statistics for:

- the exposure
- the primary outcome
- one auxiliary trait

If several candidate auxiliary traits are available, the same exposure summary
statistics are combined with the primary outcome and all candidate auxiliary
outcomes in order to screen for the most informative auxiliary trait.

### Core objects

For most analyses you will prepare:

- `BetaXG`: numeric vector of SNP-exposure effects
- `seBetaXG`: numeric vector of SNP-exposure standard errors
- `BetaYG_matrix`: matrix of SNP-outcome effects
- `seBetaYG_matrix`: matrix of SNP-outcome standard errors

The rows of all objects must refer to the same SNPs in the same order.

### Recommended matrix layout

For coheterogeneity screening, place the primary outcome in one column and each
candidate auxiliary trait in the remaining columns.

For example:

```r
colnames(BetaYG_matrix)
# [1] "primary_trait" "aux_trait_1" "aux_trait_2" "aux_trait_3"
```

The same column order should be used in `seBetaYG_matrix`.

## Typical Workflow

### Scenario 1: One primary outcome and one auxiliary trait

This is the simplest use case. The auxiliary trait has already been selected.

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

This is the principal screening use case.

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

A practical decision rule is:

- first prefer auxiliaries with `flag == "OK"`
- then prioritize pairs with small `p_value`
- among those, choose the largest `abs(rho)`

## Understanding `coheterogeneity_Q()`

`coheterogeneity_Q()` uses a guarded theoretical method with several built-in
protections:

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

After screening, the primary outcome and the selected auxiliary trait can be
analyzed jointly using a two-outcome multidimensional mode-based estimator.

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

`IBMODE()` is appropriate when the goal is robust joint mode-based estimation
across the primary and auxiliary traits.

## Using the Selected Auxiliary in `IBPRESSO()`

`IBPRESSO()` expects a data frame and column names rather than matrices. Once
an auxiliary trait has been selected, build a SNP-level table containing:

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

`IBPRESSO()` is appropriate when the goal is to detect pleiotropic outliers and
assess whether removing them materially changes the primary MR estimate.

## End-to-End Example

This is the typical end-to-end workflow in practice:

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

For a more conservative screening rule, rank only auxiliaries with:

- `flag == "OK"`
- finite `p_value`
- `p_value < 0.05`

and then choose the largest `abs(rho)` among them.

## Toy Example With Simulated Summary Statistics

The package includes a ready-to-use toy dataset called `toy_ibmr_example`.

Load it with:

```r
library(IBMR)
data("toy_ibmr_example")
```

This object is a named list containing:

- one exposure
- one primary outcome
- two candidate auxiliary traits

You can inspect its contents with:

```r
names(toy_ibmr_example)
str(toy_ibmr_example$simulation_setup)
```

### Simulation setup

The toy dataset is generated from summary-level effects for 80 SNPs.

- `BetaXG` is simulated as the SNP-exposure association vector.
- `primary_trait` is generated to depend on both the exposure signal and a
  shared latent component.
- `aux_trait_1` is generated to share that same latent component, so it should
  usually have stronger coheterogeneity with the primary trait.
- `aux_trait_2` is generated with much weaker shared structure, so it acts as a
  less suitable auxiliary candidate.

This setup is deliberately simple, but it reflects the central screening idea
of the package: among multiple auxiliary outcomes, we seek the trait that most
closely tracks the heterogeneity pattern observed for the primary outcome.

### Step 1: Extract the toy data

```r
BetaXG <- toy_ibmr_example$BetaXG
seBetaXG <- toy_ibmr_example$seBetaXG
BetaYG_matrix <- toy_ibmr_example$BetaYG_matrix
seBetaYG_matrix <- toy_ibmr_example$seBetaYG_matrix

colnames(BetaYG_matrix)
toy_ibmr_example$primary_trait
toy_ibmr_example$candidate_auxiliaries
```

### Step 2: Compute coheterogeneity

```r
cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  F_min = 5,
  min_K_pair = 20
)

round(cohet_res$rho, 3)
round(cohet_res$p_value, 3)
cohet_res$flag
```

A typical result will show that `aux_trait_1` has a stronger coheterogeneity
with `primary_trait` than `aux_trait_2`, because `aux_trait_1` was simulated to
share more of the same latent structure.

### Step 3: Rank auxiliary traits for the primary outcome

```r
ranking <- data.frame(
  aux_trait = colnames(cohet_res$rho),
  rho = cohet_res$rho["primary_trait", ],
  p_value = cohet_res$p_value["primary_trait", ],
  flag = cohet_res$flag["primary_trait", ],
  stringsAsFactors = FALSE
)

ranking <- subset(ranking, aux_trait != "primary_trait")
ranking$abs_rho <- abs(ranking$rho)
ranking <- ranking[order(-ranking$abs_rho), ]
ranking
```

The auxiliary trait would typically be chosen to satisfy the following:

- has the largest absolute `rho`
- has a usable flag such as `"OK"`
- and, if desired, also has a small `p_value`

In this toy dataset, that should usually be `toy_ibmr_example$recommended_auxiliary`,
which is `aux_trait_1`.

### Step 4: Run `IBMODE()` with the selected auxiliary

```r
chosen_aux <- toy_ibmr_example$recommended_auxiliary

BetaYG_mode <- BetaYG_matrix[, c("primary_trait", chosen_aux), drop = FALSE]
seBetaYG_mode <- seBetaYG_matrix[, c("primary_trait", chosen_aux), drop = FALSE]

ibmode_res <- IBMODE(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_mode,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_mode,
  phi = c(1, 0.5),
  n_boot = 200
)

ibmode_res
```

This returns joint mode-based estimates for the primary trait and the selected
auxiliary trait.

### Step 5: Run `IBPRESSO()` with the same auxiliary

```r
dat_ibpresso <- toy_ibmr_example$dat_ibpresso_aux1

ibpresso_res <- IBPRESSO(
  BetaOutcome = "beta_primary",
  BetaExposure = "beta_exposure",
  BetaAux = "beta_aux",
  SdOutcome = "se_primary",
  SdExposure = "se_exposure",
  SdAux = "se_aux",
  data = dat_ibpresso,
  OUTLIERtest = TRUE,
  DISTORTIONtest = FALSE,
  NbDistribution = 200,
  seed = 123
)

ibpresso_res
```

This gives an instrument-borrowing MR-PRESSO analysis using the same auxiliary
trait that was selected by coheterogeneity screening.

### What this toy example is meant to demonstrate

This example is not intended to represent a realistic full-scale GWAS
simulation. Its purpose is to provide a transparent and reproducible
illustration of the package workflow:

1. load summary statistics
2. compare the primary outcome against candidate auxiliaries
3. choose the auxiliary with the strongest coheterogeneity signal
4. use that same auxiliary in `IBMODE()` or `IBPRESSO()`

This is the principal practical use case that `IBMR` is designed to support.

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
