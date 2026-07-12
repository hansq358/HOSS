# MSPE and variance-component study (Table 4, Supplementary Tables S11 and S14).
# Run from this folder:  Rscript MSPE_Variance.R
# Results are written to ./results/. Optional flags, e.g. --B 300 --jobs 8,
# are listed in run_prediction() in common.R.
source("common.R")
run_prediction_case("case4")  # replace "case4" with "case1", "case2", or "case3" for the other cases
