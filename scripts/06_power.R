# 06_power.R
# BST 693 - statistical power for each cataract outcome, and the UKB estimates
#            converted from the linear scale to the log-odds scale
#
# Power: approximation for a binary outcome from Burgess S. Int J Epidemiol
#   2014;43:922-929 (the same inputs the mRnd website asks for):
#   z = beta_per_SD * sqrt(N * R2 * K * (1 - K));  power = Phi(z - 1.96)
#   N = outcome sample size, K = proportion of cases, R2 = variance in the
#   exposure explained by the instruments (taken from steiger_main.csv).
#   For BMI, beta_per_SD = log(OR per SD of BMI).
#   For T2D the exposure is genetic liability on the log-odds scale, whose SD
#   is pi/sqrt(3) = 1.81 under a logistic model; so OR per SD of liability =
#   (OR per unit log-odds)^1.81. This is an approximation and is stated as such.
#
# UKB conversion: the Neale-lab GWAS (ukb-d-20002_1278) used a linear model on
#   a 0/1 outcome. Its betas can be converted to the log-odds scale
#   approximately by dividing by K(1-K) (Lloyd-Jones et al. 2018, Genetics).
#   This lets the UKB estimate be placed next to FinnGen on the same scale.

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")

outcomes <- data.frame(
  id    = c("bbj-a-94", "ukb-d-20002_1278", "finn-b-H7_CATARACTSENILE", "ebi-a-GCST90018814"),
  label = c("Biobank Japan", "UKB self-report", "FinnGen R5", "Sakaue 2021 EUR"),
  N     = c(212453, 361141, 216362, 491877),
  cases = c(24622, 5045, 26758, 39519),
  stringsAsFactors = FALSE)
outcomes$K <- outcomes$cases / outcomes$N

st  <- read.csv(file.path(dir_res, "sensitivity", "steiger_main.csv"), stringsAsFactors = FALSE)
res <- read.csv(file.path(dir_res, "mr_all_methods.csv"), stringsAsFactors = FALSE)

sd_liab <- pi / sqrt(3)
# effect sizes to evaluate, on the scale each exposure is reported on
effects <- list(
  "ebi-a-GCST006867" = c("OR 1.05 per log-odds" = 1.05, "OR 1.09 (FinnGen-size)" = 1.09,
                         "OR 1.14 (BBJ-size)" = 1.14, "OR 1.20 per log-odds" = 1.20),
  "ieu-b-40"         = c("OR 1.05 per SD" = 1.05, "OR 1.13 (FinnGen-size)" = 1.13,
                         "OR 1.21 (BBJ-size)" = 1.21, "OR 1.30 per SD" = 1.30))
exp_label <- c("ebi-a-GCST006867" = "T2D", "ieu-b-40" = "BMI")

rows <- list()
for (ex in names(effects)) {
  for (i in seq_len(nrow(outcomes))) {
    oc <- outcomes$id[i]
    r2 <- st$snp_r2.exposure[st$id.exposure == ex & st$id.outcome == oc]
    if (length(r2) != 1) stop("R2 not found for ", ex, " / ", oc)
    N <- outcomes$N[i]; K <- outcomes$K[i]
    scale_fac <- if (ex == "ebi-a-GCST006867") sd_liab else 1
    for (nm in names(effects[[ex]])) {
      b_sd <- log(effects[[ex]][[nm]]) * scale_fac
      z <- b_sd * sqrt(N * r2 * K * (1 - K))
      rows[[length(rows) + 1]] <- data.frame(
        exposure = exp_label[[ex]], outcome = outcomes$label[i], outcome_id = oc,
        N = N, cases = outcomes$cases[i], K = round(K, 4), R2_instruments = signif(r2, 3),
        effect = nm, beta_per_SD = round(b_sd, 3), power = round(pnorm(z - qnorm(0.975)), 2))
    }
    # minimum detectable OR at 80% power
    b_min <- (qnorm(0.975) + qnorm(0.80)) / sqrt(N * r2 * K * (1 - K))
    or_min <- if (ex == "ebi-a-GCST006867") exp(b_min / sd_liab) else exp(b_min)
    rows[[length(rows) + 1]] <- data.frame(
      exposure = exp_label[[ex]], outcome = outcomes$label[i], outcome_id = oc,
      N = N, cases = outcomes$cases[i], K = round(K, 4), R2_instruments = signif(r2, 3),
      effect = "minimum detectable OR at 80% power", beta_per_SD = round(b_min, 3), power = round(or_min, 3))
  }
}
power_tab <- do.call(rbind, rows)
write.csv(power_tab, file.path(dir_res, "power_calculations.csv"), row.names = FALSE)

# ---- UKB linear-scale IVW estimates converted to the log-odds scale ----------
K_ukb <- outcomes$K[outcomes$id == "ukb-d-20002_1278"]
ukb <- res[res$id.outcome == "ukb-d-20002_1278" & res$method == "Inverse variance weighted", ]
ukb_conv <- data.frame(
  exposure = sub(" \\|\\|.*", "", ukb$exposure),
  b_linear = ukb$b, se_linear = ukb$se, p = ukb$pval,
  b_logodds = ukb$b / (K_ukb * (1 - K_ukb)),
  se_logodds = ukb$se / (K_ukb * (1 - K_ukb)))
ukb_conv$OR <- exp(ukb_conv$b_logodds)
ukb_conv$OR_low  <- exp(ukb_conv$b_logodds - 1.96 * ukb_conv$se_logodds)
ukb_conv$OR_high <- exp(ukb_conv$b_logodds + 1.96 * ukb_conv$se_logodds)
write.csv(ukb_conv, file.path(dir_res, "ukb_logodds_conversion.csv"), row.names = FALSE)

# ---- Printout ----------------------------------------------------------------
cat("Power by outcome (Burgess 2014 approximation):\n")
print(power_tab[, c("exposure", "outcome", "cases", "R2_instruments", "effect", "power")], row.names = FALSE)
cat("\nUKB self-report estimates on the log-odds scale (approximate conversion):\n")
print(ukb_conv[, c("exposure", "OR", "OR_low", "OR_high", "p")], row.names = FALSE, digits = 3)
