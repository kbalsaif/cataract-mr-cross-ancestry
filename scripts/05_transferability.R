# 05_transferability.R
# BST 693 - do the European instruments work the same way in East Asians?
#
# What this does:
#   Takes the European T2D and BMI instruments and asks how strongly the same
#   SNPs affect T2D and BMI when measured in East Asian samples. The IVW slope
#   is the scaling factor: 1 means the per-allele effects match; below 1 means
#   the SNPs have weaker effects in East Asians; above 1, stronger.
#   Each download is saved, so re-running skips anything already fetched.

library(TwoSampleMR)

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_tr   <- file.path(dir_res, "transferability")
dir.create(dir_tr, showWarnings = FALSE)

exp_dat <- readRDS(file.path(base_dir, "data", "exposure_instruments.rds"))

pairs <- data.frame(
  exposure = c("ebi-a-GCST006867", "ebi-a-GCST006867", "ieu-b-40"),
  outcome  = c("bbj-a-153", "ebi-a-GCST010118", "bbj-a-1"),
  label    = c("T2D EUR instruments -> T2D in BBJ (Ishigaki 2019)",
               "T2D EUR instruments -> T2D East Asian meta (Spracklen 2020)",
               "BMI EUR instruments -> BMI in BBJ (Ishigaki 2019)"),
  stringsAsFactors = FALSE)

# check the IDs exist in ao before touching the server
ao <- readRDS(file.path(base_dir, "data", "ao.rds"))
stopifnot(all(pairs$outcome %in% ao$id))

all_res <- list()
for (i in seq_len(nrow(pairs))) {
  f  <- file.path(dir_tr, paste0(pairs$exposure[i], "__", pairs$outcome[i], ".rds"))
  ed <- exp_dat[exp_dat$id.exposure == pairs$exposure[i], ]
  if (!file.exists(f)) {
    message("\nDownloading: ", pairs$label[i])
    od <- tryCatch(extract_outcome_data(snps = unique(ed$SNP), outcomes = pairs$outcome[i]),
                   error = function(e) e)
    if (inherits(od, "error")) {
      message("  Download failed: ", conditionMessage(od), "\n  Re-run later to retry this pair.")
      next
    }
    saveRDS(od, f)
    Sys.sleep(10)
  }
  od  <- readRDS(f)
  dat <- harmonise_data(ed, od)
  saveRDS(dat, file.path(dir_tr, paste0("harmonised__", pairs$outcome[i], ".rds")))

  r   <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  het <- mr_heterogeneity(dat, method_list = "mr_ivw")
  r$label      <- pairs$label[i]
  r$Q_pval_ivw <- het$Q_pval[1]
  r$snps_found <- sum(dat$mr_keep)
  r$snps_total <- nrow(ed)
  all_res[[length(all_res) + 1]] <- r
}

res <- do.call(rbind, all_res)
write.csv(res, file.path(dir_tr, "transferability_results.csv"), row.names = FALSE)

cat("\nTransferability slopes (1 = same per-allele effect in both ancestries):\n")
print(res[, c("label", "method", "snps_found", "snps_total", "b", "se", "pval", "Q_pval_ivw")],
      row.names = FALSE, digits = 3)
