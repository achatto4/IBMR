# IBMR

`IBMR` is an R package for instrument borrowing in Mendelian randomization
(MR) using summary-level genetic association data.

The motivating setting is one in which a primary outcome may share overlapping
valid instruments, or similar mechanisms of instrument invalidity, with one or
more related auxiliary outcome traits for a given exposure. `IBMR` provides
tools to identify such auxiliary traits and to incorporate them into robust MR
analyses.

The package currently provides three main functions:

- `coheterogeneity_Q()` for coheterogeneity-based auxiliary trait screening
- `IBMODE()` for multidimensional mode-based estimation with instrument borrowing
- `IBPRESSO()` for MR-PRESSO with an auxiliary trait

## Overview

Mendelian randomization is widely used to estimate causal effects of exposures
on outcomes using genetic variants as instrumental variables. Standard MR
methods can lose power or become biased when many candidate instruments are
invalid. `IBMR` addresses this problem by leveraging related outcome traits that
share informative heterogeneity structure with the primary outcome.

The typical workflow is:

1. Assemble SNP-level summary statistics for the exposure, the primary outcome,
   and one or more candidate auxiliary outcomes.
2. Use `coheterogeneity_Q()` to quantify coheterogeneity between the primary
   outcome and each auxiliary trait.
3. Select the most informative auxiliary trait.
4. Carry the selected auxiliary trait into `IBMODE()` or `IBPRESSO()` for
   downstream robust MR analysis.

## Graphical Overview

![IBMR graphical summary](man/figures/ibmr-graphical-summary.png)

The figure summarizes the conceptual workflow implemented in `IBMR`.

- Panel A illustrates the motivating MR setting with valid and invalid
  candidate instruments.
- Panel B illustrates auxiliary-trait selection using coheterogeneity.
- Panel C shows the summary-statistics inputs and ratio estimates used by the
  methods.
- Panel D illustrates the downstream joint-analysis setting used by `IBMODE()`
  and `IBPRESSO()`.

## Installation

Install from GitHub with:

```r
install.packages("devtools")
library(devtools)
devtools::install_github("achatto4/IBMR")
library(IBMR)
```

For local development:

```r
install.packages(c("MASS", "ks"))
devtools::install(".")
library(IBMR)
```

## Example: body mass index (BMI) and coronary artery disease (CAD)

We estimate the effect of BMI on CAD, borrowing instruments from type 2 diabetes
(T2D) as an auxiliary trait. The bundled dataset `ibmr_example` is a simulated
illustration (not real GWAS data): BMI has a true positive effect on CAD, and a
subset of invalid instruments show pleiotropy shared with T2D, so T2D is an
informative auxiliary trait.

#### Step 1. Screen the auxiliary trait with coheterogeneity

```r
library(IBMR)
data("ibmr_example")
dat <- ibmr_example

cohet_res <- coheterogeneity_Q(
  BetaXG          = dat$BetaXG,
  BetaYG_matrix   = dat$BetaYG_matrix,
  seBetaXG        = dat$seBetaXG,
  seBetaYG_matrix = dat$seBetaYG_matrix,
  F_min = 5,
  min_K_pair = 20
)

round(cohet_res$rho, 3)   # coheterogeneity of CAD with the T2D auxiliary
cohet_res$flag
#> coheterogeneity(CAD, T2D) ~ 0.82  (strong shared pleiotropic structure)
```

A large, significant coheterogeneity indicates that CAD and T2D share invalid
(pleiotropic) instruments for BMI, so T2D is an informative auxiliary trait.

#### Step 2. Estimate the causal effect with instrument borrowing

```r
ibp <- IBPRESSO(
  BetaOutcome  = "BetaOutcome",
  BetaExposure = "BetaExposure",
  BetaAux      = "BetaAux",
  SdOutcome    = "SdOutcome",
  SdExposure   = "SdExposure",
  SdAux        = "SdAux",
  data         = dat$dat_ibpresso_aux1,   # CAD primary, T2D auxiliary
  OUTLIERtest  = TRUE,
  seed         = 1
)

c(estimate = ibp$corrected_beta,
  se        = ibp$corrected_se,
  p         = ibp$p_value,
  n_outliers = ibp$n_outliers)
#> corrected_beta ~ 0.40 after removing 60 shared-pleiotropy instruments
#> flagged by borrowing from T2D (recovers the true +0.40 simulated effect).
```

`ibp$corrected_beta` is the outlier-corrected effect of BMI on CAD after
borrowing instruments from T2D; `ibp$outlier_idx` lists the flagged pleiotropic
instruments.

#### Step 3 (alternative). Mode-based estimator

`IBMODE()` provides a mode-based alternative that takes the same summary
statistics directly (primary and auxiliary outcomes as the two columns of
`BetaYG_matrix`):

```r
ibmode <- IBMODE(
  BetaXG          = dat$BetaXG,
  BetaYG_matrix   = dat$BetaYG_matrix,   # columns: CAD (primary), T2D (auxiliary)
  seBetaXG        = dat$seBetaXG,
  seBetaYG_matrix = dat$seBetaYG_matrix,
  phi   = c(1, 0.5),
  n_boot = 200,
  seed  = 1
)

ibmode[, c("phi", "Estimate_CAD", "SE_CAD", "P_CAD")]
#> Estimate_CAD ~ 0.40 at both phi values (agrees with IBPRESSO and the truth).
```

## Included Example Data

The package bundles two datasets:

- `ibmr_example` — simulated illustrative example (BMI &rarr; CAD, T2D
  auxiliary), used in the example above.
- `toy_ibmr_example` — a small simulated example with two candidate auxiliary
  traits, for quick testing.

```r
data("ibmr_example")
names(ibmr_example)
```

This example is intended to illustrate the core workflow rather than to serve
as a realistic full-scale GWAS simulation.

## Vignettes

Longer tutorials are provided in the package vignettes:

- `vignettes/auxiliary-selection.Rmd`: coheterogeneity-based screening of
  candidate auxiliary traits
- `vignettes/instrument-borrowing-workflow.Rmd`: end-to-end workflow using the
  selected auxiliary trait in `IBMODE()` and `IBPRESSO()`

These vignettes expand on the screening logic, toy data structure,
interpretation of outputs, and downstream workflow.

## Input Checks

Before running the package, it is good practice to confirm that:

- SNP order matches across all vectors and matrices
- alleles are harmonized across exposure and outcome traits
- standard errors are finite and positive
- weak instruments are handled appropriately
- missing data are addressed consistently

## Citation

If you use this package in applied work, please cite the repository and the
relevant associated method paper.
