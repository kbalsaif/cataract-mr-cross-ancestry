# 09_mvmr_strength.R
# BST 693 - conditional instrument strength and heterogeneity for the
#            multivariable MR (T2D + BMI) fitted in script 08
#
# In multivariable MR the question is not whether the SNPs predict each
# exposure on its own, but whether they predict each exposure CONDITIONAL on
# the other (conditional F-statistic; Sanderson et al. 2019, 2021). Values
# below 10 signal weak-instrument bias. The MVMR package computes this, plus
# a Q-statistic for heterogeneity. It reads the mvdat objects saved by script 08.
# Install once:  remotes::install_github("WSpiller/MVMR")

library(MVMR)
base_dir <- "~/Desktop/BST693"
dir_mv   <- file.path(base_dir, "results", "mvmr")
outcomes <- c("bbj-a-94", "ukb-d-20002_1278", "finn-b-H7_CATARACTSENILE", "ebi-a-GCST90018814")

out <- list()
for (oc in outcomes) {
  mvdat <- readRDS(file.path(dir_mv, paste0("mvdat__", oc, ".rds")))
  # columns of exposure_beta are in the order of mvdat$expname
  r_in <- format_mvmr(BXGs = mvdat$exposure_beta, BYG = mvdat$outcome_beta,
                      seBXGs = mvdat$exposure_se, seBYG = mvdat$outcome_se,
                      RSID = rownames(mvdat$exposure_beta))
  fs  <- strength_mvmr(r_in, gencov = 0)      # gencov = 0 assumes no sample overlap between exposure GWAS
  q   <- pleiotropy_mvmr(r_in, gencov = 0)
  est <- ivw_mvmr(r_in)                        # cross-check of the script-08 estimates
  out[[oc]] <- list(cond_F = fs, Q = q, ivw = est)
  cat("\n==", oc, "==\n")
  cat("Exposure order:", paste(mvdat$expname$exposure, collapse = " | "), "\n")
  cat("Conditional F:\n"); print(fs)
  cat("Heterogeneity Q:\n"); print(q)
  cat("IVW-MVMR estimates (should match script 08):\n"); print(est)
}
saveRDS(out, file.path(dir_mv, "mvmr_strength_heterogeneity.rds"))
