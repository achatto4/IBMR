## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(IBMR)


## -----------------------------------------------------------------------------
data("toy_ibmr_example")

BetaXG <- toy_ibmr_example$BetaXG
seBetaXG <- toy_ibmr_example$seBetaXG
BetaYG_matrix <- toy_ibmr_example$BetaYG_matrix
seBetaYG_matrix <- toy_ibmr_example$seBetaYG_matrix
primary_name <- toy_ibmr_example$primary_trait
candidate_aux <- toy_ibmr_example$candidate_auxiliaries


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
rho_primary <- cohet_res$rho[primary_name, ]
p_primary <- cohet_res$p_value[primary_name, ]
flag_primary <- cohet_res$flag[primary_name, ]

rho_primary
p_primary
flag_primary


## -----------------------------------------------------------------------------
ranking <- data.frame(
  aux_trait = candidate_aux,
  rho = rho_primary[candidate_aux],
  p_value = p_primary[candidate_aux],
  flag = flag_primary[candidate_aux],
  stringsAsFactors = FALSE
)

ranking$abs_rho <- abs(ranking$rho)
ranking <- ranking[order(-ranking$abs_rho), ]
ranking


## -----------------------------------------------------------------------------
chosen_aux <- ranking$aux_trait[1]
chosen_aux


## -----------------------------------------------------------------------------
toy_ibmr_example$recommended_auxiliary


## -----------------------------------------------------------------------------
BetaYG_mode <- BetaYG_matrix[, c(primary_name, chosen_aux), drop = FALSE]
seBetaYG_mode <- seBetaYG_matrix[, c(primary_name, chosen_aux), drop = FALSE]


## ----eval = FALSE-------------------------------------------------------------
# ibmode_res <- IBMODE(
#   BetaXG = BetaXG,
#   BetaYG_matrix = BetaYG_mode,
#   seBetaXG = seBetaXG,
#   seBetaYG_matrix = seBetaYG_mode,
#   phi = c(1, 0.5),
#   n_boot = 200,
#   seed = 123
# )
# 
# ibmode_res


## -----------------------------------------------------------------------------
dat_ibpresso <- toy_ibmr_example$dat_ibpresso_aux1
head(dat_ibpresso)


## ----eval = FALSE-------------------------------------------------------------
# ibpresso_res <- IBPRESSO(
#   BetaOutcome = "beta_primary",
#   BetaExposure = "beta_exposure",
#   BetaAux = "beta_aux",
#   SdOutcome = "se_primary",
#   SdExposure = "se_exposure",
#   SdAux = "se_aux",
#   data = dat_ibpresso,
#   OUTLIERtest = TRUE,
#   DISTORTIONtest = FALSE,
#   NbDistribution = 200,
#   seed = 123
# )
# 
# ibpresso_res

