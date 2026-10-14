# 03_steiger_direction.R
# BST 693 - Steiger direction test for the 8 main pairs (T2D and BMI x 4 outcomes)
#
# What this does:
#   The Steiger test asks: do the SNPs explain more of the exposure than of the
#   outcome? If yes, the causal direction exposure -> cataract is supported.
#   For binary traits (T2D, and the cataract outcomes on the log-odds scale) the
#   SNP-trait correlation is computed from the log odds ratio, allele frequency,
#   case/control counts and an assumed population prevalence (get_r_from_lor).
#   For continuous or linear-scale traits (BMI; UKB self-report cataract) it is
#   computed from the p-value and sample size (get_r_from_pn).
#   The test is then repeated over a grid of prevalence assumptions.

library(TwoSampleMR)

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_sens <- file.path(dir_res, "sensitivity")

dat_all <- readRDS(file.path(dir_res, "harmonised_all.rds"))
main_exposures <- c("ebi-a-GCST006867", "ieu-b-40")
dat <- dat_all[dat_all$id.exposure %in% main_exposures & dat_all$mr_keep, ]

# ---- Case/control counts (from ao; two fixes: T2D ncontrol corrected to
#      655,666 - 61,714 = 593,952; FinnGen total filled as cases + controls) ----
n_info <- data.frame(
  id       = c("ebi-a-GCST006867", "bbj-a-94", "ukb-d-20002_1278",
               "finn-b-H7_CATARACTSENILE", "ebi-a-GCST90018814"),
  ncase    = c(61714, 24622, 5045, 26758, 39519),
  ncontrol = c(593952, 187831, 356096, 189604, 452358),
  stringsAsFactors = FALSE)
n_info$n <- n_info$ncase + n_info$ncontrol

# ---- Prevalence assumptions -------------------------------------------------
# T2D: 0.092 = IDF Diabetes Atlas 2021, Europe region, adults 20-79 (Sun et al. 2022, Table 3)
# Cataract: 0.172 = global age-standardised pooled prevalence of any cataract
#           (Hashemi et al. 2020, Eye). Replace with regional values if obtained.
prev_t2d_main <- 0.092
prev_cat_main <- c("bbj-a-94" = 0.172,
                   "finn-b-H7_CATARACTSENILE" = 0.172,
                   "ebi-a-GCST90018814" = 0.172)
lor_outcomes <- names(prev_cat_main)          # outcomes on the log-odds scale
lin_outcomes <- "ukb-d-20002_1278"            # linear-scale outcome

# ---- Data checks -----------------------------------------------------------
cat("Rows per pair:\n"); print(table(dat$id.exposure, dat$id.outcome))
cat("\nMissing allele frequency in outcome, by outcome:\n")
print(tapply(is.na(dat$eaf.outcome), dat$id.outcome, sum))
cat("\nMissing sample size in exposure, by exposure (before fill):\n")
print(tapply(is.na(dat$samplesize.exposure), dat$id.exposure, sum))

# ---- Function: add SNP-trait correlations ----------------------------------
add_r <- function(d, prev_t2d, prev_cat) {
  for (i in seq_len(nrow(n_info))) {
    d$samplesize.exposure[d$id.exposure == n_info$id[i]] <- n_info$n[i]
    d$samplesize.outcome[d$id.outcome == n_info$id[i]]   <- n_info$n[i]
  }
  d$r.exposure <- NA_real_
  d$r.outcome  <- NA_real_

  # Exposure: T2D (binary, log odds) vs BMI (continuous)
  i <- d$id.exposure == "ebi-a-GCST006867"
  ni <- n_info[n_info$id == "ebi-a-GCST006867", ]
  d$r.exposure[i] <- get_r_from_lor(d$beta.exposure[i], d$eaf.exposure[i],
                                    ni$ncase, ni$ncontrol, prev_t2d)
  i <- d$id.exposure == "ieu-b-40"
  d$r.exposure[i] <- get_r_from_pn(d$pval.exposure[i], d$samplesize.exposure[i])

  # Outcomes on the log-odds scale
  for (oc in lor_outcomes) {
    i  <- d$id.outcome == oc
    ni <- n_info[n_info$id == oc, ]
    af <- ifelse(is.na(d$eaf.outcome[i]), d$eaf.exposure[i], d$eaf.outcome[i])
    d$r.outcome[i] <- get_r_from_lor(d$beta.outcome[i], af, ni$ncase, ni$ncontrol, prev_cat[[oc]])
  }
  # Linear-scale outcome (UKB self-report)
  i <- d$id.outcome %in% lin_outcomes
  d$r.outcome[i] <- get_r_from_pn(d$pval.outcome[i], d$samplesize.outcome[i])

  d
}

# ---- Main run ----------------------------------------------------------------
d1 <- add_r(dat, prev_t2d_main, prev_cat_main)
drop <- is.na(d1$r.exposure) | is.na(d1$r.outcome)
cat("\nSNPs dropped because a correlation could not be computed:", sum(drop), "\n")
d1 <- d1[!drop, ]

st <- directionality_test(d1)
st <- as.data.frame(st)
write.csv(st, file.path(dir_sens, "steiger_main.csv"), row.names = FALSE)

# ---- Sensitivity over prevalence assumptions -------------------------------
grid <- expand.grid(prev_t2d = c(0.082, 0.092, 0.105),
                    prev_cat = c(0.10, 0.172, 0.25))
sens <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
  pc <- setNames(rep(grid$prev_cat[k], length(lor_outcomes)), lor_outcomes)
  dk <- add_r(dat, grid$prev_t2d[k], pc)
  dk <- dk[!(is.na(dk$r.exposure) | is.na(dk$r.outcome)), ]
  s  <- as.data.frame(suppressMessages(directionality_test(dk)))
  cbind(prev_t2d = grid$prev_t2d[k], prev_cat = grid$prev_cat[k],
        s[, c("id.exposure", "id.outcome", "snp_r2.exposure", "snp_r2.outcome",
              "correct_causal_direction", "steiger_pval")])
}))
write.csv(sens, file.path(dir_sens, "steiger_sensitivity.csv"), row.names = FALSE)

# ---- Printout to paste back -------------------------------------------------
cat("\nSteiger test, main assumptions (T2D prevalence 9.2%, cataract 17.2%):\n")
print(st[, c("id.exposure", "id.outcome", "snp_r2.exposure", "snp_r2.outcome",
             "correct_causal_direction", "steiger_pval")], row.names = FALSE, digits = 3)

cat("\nAcross all 9 prevalence settings, was the direction always exposure -> cataract?\n")
print(aggregate(correct_causal_direction ~ id.exposure + id.outcome, data = sens, FUN = all))

cat("\nRange of variance explained in the outcome across settings:\n")
print(aggregate(snp_r2.outcome ~ id.exposure + id.outcome, data = sens,
                FUN = function(x) signif(range(x), 3)))
