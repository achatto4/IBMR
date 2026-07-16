# Build the bundled illustrative example dataset for IBMR.
#
# SIMULATED example (not real GWAS data): one exposure, one primary outcome, and
# two candidate auxiliary outcomes. Auxiliary_1 shares most of its invalid-
# instrument (pleiotropic) structure with the primary outcome (high
# coheterogeneity, the informative auxiliary); Auxiliary_2 shares little (low
# coheterogeneity). Instruments are strong so both IBMODE and IBPRESSO behave
# stably. Raw values come from data-raw/gen_example.py (seed 2025). Run once from
# the package root:  Rscript data-raw/build_ibmr_example.R

d <- read.csv("data-raw/sim_ibmr_example.csv", stringsAsFactors = FALSE)

BetaYG_matrix   <- as.matrix(d[, c("Beta_Primary", "Beta_Aux1", "Beta_Aux2")])
seBetaYG_matrix <- as.matrix(d[, c("se_Primary", "se_Aux1", "se_Aux2")])
colnames(BetaYG_matrix)   <- c("Primary", "Auxiliary_1", "Auxiliary_2")
colnames(seBetaYG_matrix) <- c("Primary", "Auxiliary_1", "Auxiliary_2")
rownames(BetaYG_matrix)   <- d$SNP

# IBPRESSO input using the recommended (most coheterogeneous) auxiliary, Auxiliary_1
dat_ibpresso_aux1 <- data.frame(
  SNP          = d$SNP,
  BetaOutcome  = d$Beta_Primary,   # primary outcome
  BetaExposure = d$BetaXG,         # exposure
  BetaAux      = d$Beta_Aux1,      # selected auxiliary outcome
  SdOutcome    = d$se_Primary,
  SdExposure   = d$seBetaXG,
  SdAux        = d$se_Aux1,
  stringsAsFactors = FALSE
)

ibmr_example <- list(
  BetaXG                = d$BetaXG,
  seBetaXG              = d$seBetaXG,
  BetaYG_matrix         = BetaYG_matrix,
  seBetaYG_matrix       = seBetaYG_matrix,
  primary_trait         = "Primary",
  candidate_auxiliaries = c("Auxiliary_1", "Auxiliary_2"),
  recommended_auxiliary = "Auxiliary_1",
  dat_ibpresso_aux1     = dat_ibpresso_aux1,
  source_note = paste(
    "Simulated illustrative example (not real GWAS data): one exposure, one",
    "primary outcome, and two candidate auxiliary outcomes;", nrow(d),
    "strong instruments. Auxiliary_1 shares pleiotropic structure with the",
    "primary outcome (informative); Auxiliary_2 does not."
  )
)

usethis::use_data(ibmr_example, overwrite = TRUE)
cat("Saved data/ibmr_example.rda with", nrow(d), "instruments\n")
