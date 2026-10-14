# 08_mvmr_t2d_bmi.R
# BST 693 - multivariable MR: T2D and BMI together, on each cataract outcome
#
# Why: Arruda 2025 (T2D effect survives adjustment for BMI, Europeans) and
# Liu 2026 (BMI effect disappears after adjustment for T2D, Europeans) suggest
# T2D is the proximal pathway. Running the same model in both ancestries shows
# whether that holds in East Asians. Uses the TwoSampleMR multivariable functions
# (mv_extract_exposures -> extract_outcome_data -> mv_harmonise_data -> mv_multiple).
# Server is contacted once per outcome; each outcome is saved and skipped on re-run.

library(TwoSampleMR)
base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_mv   <- file.path(dir_res, "mvmr")
dir.create(dir_mv, showWarnings = FALSE)

exposures <- c("ebi-a-GCST006867", "ieu-b-40")   # T2D, BMI (IDs from the handoff; checked in ao by script 01)
outcomes  <- c("bbj-a-94", "ukb-d-20002_1278", "finn-b-H7_CATARACTSENILE", "ebi-a-GCST90018814")

exp_file <- file.path(dir_mv, "mv_exposures.rds")
if (file.exists(exp_file)) {
  mv_exp <- readRDS(exp_file)
} else {
  mv_exp <- mv_extract_exposures(exposures)   # union of instruments, re-clumped, both exposures' effects
  saveRDS(mv_exp, exp_file)
}
cat("SNPs in the multivariable instrument set:", length(unique(mv_exp$SNP)), "\n")

all_res <- list()
for (oc in outcomes) {
  f <- file.path(dir_mv, paste0("outcome__", oc, ".rds"))
  if (!file.exists(f)) {
    od <- tryCatch(extract_outcome_data(snps = unique(mv_exp$SNP), outcomes = oc), error = function(e) e)
    if (inherits(od, "error")) { message("Download failed for ", oc, ": ", conditionMessage(od)); next }
    saveRDS(od, f); Sys.sleep(10)
  }
  od <- readRDS(f)
  mvdat <- mv_harmonise_data(mv_exp, od)
  saveRDS(mvdat, file.path(dir_mv, paste0("mvdat__", oc, ".rds")))
  r <- mv_multiple(mvdat)$result
  r$outcome_id <- oc
  all_res[[length(all_res) + 1]] <- r
}
res <- do.call(rbind, all_res)
res$OR <- exp(res$b); res$OR_low <- exp(res$b - 1.96 * res$se); res$OR_high <- exp(res$b + 1.96 * res$se)
write.csv(res, file.path(dir_mv, "mvmr_results.csv"), row.names = FALSE)

cat("\nMultivariable MR (each exposure adjusted for the other):\n")
print(res[, c("exposure", "outcome_id", "nsnp", "b", "se", "pval", "OR", "OR_low", "OR_high")], row.names = FALSE, digits = 3)
cat("\nNote: UKB outcome is on a linear scale (ORs not meaningful there); East Asian arm uses European exposure effects.\n")
