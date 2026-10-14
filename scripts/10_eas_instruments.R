# 10_eas_instruments.R
# BST 693 - East Asian-discovered T2D instruments against Biobank Japan cataract
#
# Why: our East Asian-scale T2D estimate (script 07, analysis D) is OR 1.18,
# above the ~1.10 reported by Zhang 2022, who used Biobank Japan T2D instruments
# (36,614 cases; OpenGWAS bbj-a-77 matches that count). This script runs that
# design and the same with the larger East Asian T2D meta-analysis
# (Spracklen 2020, ebi-a-GCST010118), all methods, so the two can be compared.
# Caveat: both exposure GWAS overlap the outcome cohort (bbj-a-77 fully; Spracklen
# partly), so estimates are pulled toward the observational association; with
# strong instruments the bias is small but must be stated.
# Server is contacted twice for instruments and twice for outcome data; each
# download is saved and skipped on re-run.

library(TwoSampleMR)
base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_out  <- file.path(dir_res, "eas_instruments")
dir.create(dir_out, showWarnings = FALSE)

ao <- readRDS(file.path(base_dir, "data", "ao.rds"))
exposures <- c("bbj-a-77", "ebi-a-GCST010118")
stopifnot(all(exposures %in% ao$id))
cols <- intersect(c("id", "trait", "population", "ncase", "ncontrol", "author", "year"), names(ao))
print(ao[ao$id %in% exposures, cols], row.names = FALSE)
outcome <- "bbj-a-94"

all_res <- list()
for (ex in exposures) {
  f_exp <- file.path(dir_out, paste0("instruments__", ex, ".rds"))
  if (!file.exists(f_exp)) {
    # East Asian LD reference for clumping East Asian instruments
    ed <- tryCatch(extract_instruments(outcomes = ex, clump = FALSE), error = function(e) e)
    if (inherits(ed, "error")) { message("Instrument download failed for ", ex, ": ", conditionMessage(ed)); next }
    ed <- tryCatch(clump_data(ed, clump_r2 = 0.001, clump_kb = 10000, pop = "EAS"), error = function(e) e)
    if (inherits(ed, "error")) { message("Clumping failed for ", ex, ": ", conditionMessage(ed)); next }
    saveRDS(ed, f_exp); Sys.sleep(10)
  }
  ed <- readRDS(f_exp)
  cat("\n", ex, ": ", nrow(ed), " East Asian instruments after EAS clumping\n", sep = "")

  f_out <- file.path(dir_out, paste0("outcome__", ex, "__", outcome, ".rds"))
  if (!file.exists(f_out)) {
    od <- tryCatch(extract_outcome_data(snps = ed$SNP, outcomes = outcome), error = function(e) e)
    if (inherits(od, "error")) { message("Outcome download failed: ", conditionMessage(od)); next }
    saveRDS(od, f_out); Sys.sleep(10)
  }
  od  <- readRDS(f_out)
  dat <- harmonise_data(ed, od)
  saveRDS(dat, file.path(dir_out, paste0("harmonised__", ex, "__", outcome, ".rds")))

  set.seed(693)
  r <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode", "mr_raps"))
  r$Q_pval   <- mr_heterogeneity(dat, method_list = "mr_ivw")$Q_pval[1]
  r$egger_intercept_p <- mr_pleiotropy_test(dat)$pval
  r$mean_F <- mean((dat$beta.exposure[dat$mr_keep] / dat$se.exposure[dat$mr_keep])^2)
  all_res[[length(all_res) + 1]] <- r
}
res <- do.call(rbind, all_res)
res$OR <- exp(res$b); res$OR_low <- exp(res$b - 1.96 * res$se); res$OR_high <- exp(res$b + 1.96 * res$se)
write.csv(res, file.path(dir_out, "eas_instruments_results.csv"), row.names = FALSE)

cat("\nEast Asian-instrument T2D -> Biobank Japan cataract:\n")
print(res[, c("id.exposure", "method", "nsnp", "b", "se", "pval", "OR", "OR_low", "OR_high", "Q_pval", "egger_intercept_p", "mean_F")],
      row.names = FALSE, digits = 3)
cat("\nCompare with: European-instrument main analysis OR 1.14 (1.10-1.19); East Asian-scale (script 07 D) OR 1.18 (1.14-1.23); Zhang 2022 about OR 1.10 (1.06-1.17).\n")
