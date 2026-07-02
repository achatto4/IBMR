# `Weights` and `Weights1` are data-frame columns referenced via non-standard
# evaluation inside model formulas (the `weights =` argument of `lm()`) in
# IBPRESSO(). Declaring them here silences the spurious "no visible binding for
# global variable" NOTE from R CMD check without changing any behaviour.
globalVariables(c("Weights", "Weights1"))
