## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(IBMR)


## -----------------------------------------------------------------------------
data("toy_ibmr_example")

names(toy_ibmr_example)
toy_ibmr_example$primary_trait
toy_ibmr_example$candidate_auxiliaries


## -----------------------------------------------------------------------------
str(toy_ibmr_example$simulation_setup)


## -----------------------------------------------------------------------------
BetaXG <- toy_ibmr_example$BetaXG
seBetaXG <- toy_ibmr_example$seBetaXG
BetaYG_matrix <- toy_ibmr_example$BetaYG_matrix
seBetaYG_matrix <- toy_ibmr_example$seBetaYG_matrix

colnames(BetaYG_matrix)


## -----------------------------------------------------------------------------
cohet_res <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  F_min = 5,
  min_K_pair = 20
)


## -----------------------------------------------------------------------------
round(cohet_res$rho, 3)
round(cohet_res$p_value, 3)
cohet_res$flag
cohet_res$K


## -----------------------------------------------------------------------------
# Example: equal weights across instruments (a simple fixed-weight scheme).
custom_w <- rep(1, length(BetaXG))

cohet_fixed <- coheterogeneity_Q(
  BetaXG = BetaXG,
  BetaYG_matrix = BetaYG_matrix,
  seBetaXG = seBetaXG,
  seBetaYG_matrix = seBetaYG_matrix,
  weights = custom_w,
  F_min = 5,
  min_K_pair = 20
)

round(cohet_fixed$se, 4)   # fixed-weight closed-form SE
cohet_fixed$method


## -----------------------------------------------------------------------------
primary_name <- toy_ibmr_example$primary_trait
cohet_res$rho[primary_name, ]
cohet_res$p_value[primary_name, ]
cohet_res$flag[primary_name, ]


## -----------------------------------------------------------------------------
ranking <- data.frame(
  aux_trait = colnames(cohet_res$rho),
  rho = cohet_res$rho[primary_name, ],
  p_value = cohet_res$p_value[primary_name, ],
  flag = cohet_res$flag[primary_name, ],
  stringsAsFactors = FALSE
)

ranking <- subset(ranking, aux_trait != primary_name)
ranking$abs_rho <- abs(ranking$rho)
ranking <- ranking[order(-ranking$abs_rho), ]
ranking


## -----------------------------------------------------------------------------
toy_ibmr_example$recommended_auxiliary

