# 07_eas_arm_sensitivity.R
# BST 693 - checks specific to using European instruments in an East Asian outcome
#
# For T2D and BMI on Biobank Japan (and FinnGen for comparison), four checks:
#   A. Drop palindromic SNPs (strand was inferred from allele frequency, which
#      differs between ancestries) and drop proxy SNPs (chosen with a European
#      LD reference, so they may not tag the target SNP in East Asians).
#   B. Prune to one SNP per 1 Mb window (European-independent SNPs can be in
#      LD in East Asians).
#   C. Restrict BOTH arms to instruments with demonstrated East Asian exposure
#      effects (p < 0.05, and separately p < 5e-8, in the East Asian exposure
#      GWAS from script 05), then re-test the ancestry difference.
#   D. East Asian-scale estimate: same instruments, but with the East Asian
#      SNP-exposure effects as the denominator instead of the European ones.
#      Uses SNPs with East Asian exposure p < 5e-8 (to limit weak-instrument
#      bias) and reports IVW and RAPS (RAPS allows for error in the exposure
#      effects). Caveat: the East Asian exposure GWAS overlap Biobank Japan.
# Everything is computed on the same harmonised data saved by scripts 01 and 05.

library(TwoSampleMR)
base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_tr   <- file.path(dir_res, "transferability")
dir_out  <- file.path(dir_res, "eas_arm_sensitivity")
dir.create(dir_out, showWarnings = FALSE)

dat_all <- readRDS(file.path(dir_res, "harmonised_all.rds"))
dat_all <- dat_all[dat_all$mr_keep, ]
eas_exp <- list("ebi-a-GCST006867" = readRDS(file.path(dir_tr, "harmonised__ebi-a-GCST010118.rds")),
                "ieu-b-40"         = readRDS(file.path(dir_tr, "harmonised__bbj-a-1.rds")))
labs <- c("ebi-a-GCST006867" = "T2D", "ieu-b-40" = "BMI")
methods <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")

fit <- function(d, tag) {
  if (nrow(d) < 3) return(NULL)
  r <- mr(d, method_list = methods)
  data.frame(analysis = tag, exposure = labs[[d$id.exposure[1]]], outcome = d$id.outcome[1],
             method = r$method, nsnp = r$nsnp, b = r$b, se = r$se, pval = r$pval)
}
prune_1mb <- function(d) {
  d <- d[order(d$chr.exposure, d$pos.exposure), ]
  keep <- rep(TRUE, nrow(d)); last <- 1
  for (i in seq_len(nrow(d))[-1]) {
    same <- d$chr.exposure[i] == d$chr.exposure[last] && (d$pos.exposure[i] - d$pos.exposure[last]) < 1e6
    if (same) {
      if (abs(d$beta.exposure[i] / d$se.exposure[i]) > abs(d$beta.exposure[last] / d$se.exposure[last])) { keep[last] <- FALSE; last <- i } else keep[i] <- FALSE
    } else last <- i
  }
  d[keep, ]
}

res <- list()
for (ex in names(labs)) {
  te <- eas_exp[[ex]]; te <- te[te$mr_keep, ]
  eas_p   <- setNames(te$pval.outcome, te$SNP)     # East Asian exposure p-value per SNP
  eas_b   <- setNames(te$beta.outcome, te$SNP)     # East Asian exposure effect (aligned to the European effect allele)
  eas_se  <- setNames(te$se.outcome, te$SNP)
  for (oc in c("bbj-a-94", "finn-b-H7_CATARACTSENILE")) {
    d <- dat_all[dat_all$id.exposure == ex & dat_all$id.outcome == oc, ]
    res[[length(res) + 1]] <- fit(d, "0. all instruments (main analysis)")
    res[[length(res) + 1]] <- fit(d[!d$palindromic, ], "A1. no palindromic SNPs")
    res[[length(res) + 1]] <- fit(d[!(d$proxy.outcome %in% TRUE), ], "A2. no proxy SNPs")
    res[[length(res) + 1]] <- fit(prune_1mb(d), "B. one SNP per 1 Mb")
    for (thr in c(0.05, 5e-8)) {
      ok <- d$SNP %in% names(eas_p)[eas_p < thr]
      res[[length(res) + 1]] <- fit(d[ok, ], sprintf("C. East Asian exposure p < %s", format(thr)))
    }
    if (oc == "bbj-a-94") {
      ok <- d$SNP %in% names(eas_p)[eas_p < 5e-8]
      dd <- d[ok, ]
      dd$beta.exposure <- eas_b[dd$SNP]; dd$se.exposure <- eas_se[dd$SNP]
      dd$id.exposure <- paste0(ex, "_EASscale"); dd$exposure <- paste(labs[[ex]], "(East Asian effects)")
      r <- mr(dd, method_list = c("mr_ivw", "mr_raps", "mr_weighted_median"))
      res[[length(res) + 1]] <- data.frame(analysis = "D. East Asian-scale (EAS exposure effects, p < 5e-8)",
        exposure = labs[[ex]], outcome = oc, method = r$method, nsnp = r$nsnp, b = r$b, se = r$se, pval = r$pval)
    }
  }
}
out <- do.call(rbind, res)
out$OR <- exp(out$b)
write.csv(out, file.path(dir_out, "eas_arm_sensitivity.csv"), row.names = FALSE)

# Ancestry difference (IVW) under each restriction that applies to both arms
cmp <- list()
for (ex in labs) for (an in unique(out$analysis)) {
  a <- out[out$exposure == ex & out$analysis == an & out$method == "Inverse variance weighted" & out$outcome == "bbj-a-94", ]
  b <- out[out$exposure == ex & out$analysis == an & out$method == "Inverse variance weighted" & out$outcome == "finn-b-H7_CATARACTSENILE", ]
  if (nrow(a) == 1 && nrow(b) == 1) {
    z <- (a$b - b$b) / sqrt(a$se^2 + b$se^2)
    cmp[[length(cmp) + 1]] <- data.frame(exposure = ex, analysis = an, nsnp_bbj = a$nsnp, nsnp_finn = b$nsnp,
                                         b_bbj = a$b, b_finn = b$b, z = z, p_diff = 2 * pnorm(-abs(z)))
  }
}
cmp <- do.call(rbind, cmp)
write.csv(cmp, file.path(dir_out, "eas_arm_ancestry_difference.csv"), row.names = FALSE)

cat("IVW under each check:\n")
print(out[out$method == "Inverse variance weighted", c("analysis", "exposure", "outcome", "nsnp", "b", "se", "pval")], row.names = FALSE, digits = 3)
cat("\nEast Asian-scale estimates (all methods run):\n")
print(out[grepl("^D", out$analysis), c("exposure", "method", "nsnp", "b", "se", "pval")], row.names = FALSE, digits = 3)
cat("\nAncestry difference (Biobank Japan vs FinnGen, IVW) under each restriction:\n")
print(cmp, row.names = FALSE, digits = 3)
