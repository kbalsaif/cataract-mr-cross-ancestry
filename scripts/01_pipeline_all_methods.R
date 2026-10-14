# 01_pipeline_all_methods_v2.R
# BST 693 - T2D / BMI / BP / alcohol / smoking -> cataract, all MR methods
#
# What this does:
#   1. Checks every dataset ID against the `ao` table and prints the trait names.
#   2. Downloads exposure instruments once and saves them.
#   3. Downloads outcome data ONE exposure-outcome pair at a time and saves each
#      pair to its own file. If the server times out, the script stops cleanly.
#      Re-running it skips pairs already saved and continues where it left off.
#   4. Harmonises, runs every method in mr_method_list() for every pair,
#      plus Cochran's Q and the Egger intercept on the same data.
#   5. Saves everything (harmonised SNP-level data + all results) to CSV and RDS.

library(TwoSampleMR)
library(ieugwasr)

# ---- 0. Settings -------------------------------------------------------------
base_dir <- "~/Desktop/BST693"
dir_data <- file.path(base_dir, "data")
dir_raw  <- file.path(dir_data, "outcome_pieces")
dir_res  <- file.path(base_dir, "results")
for (d in c(dir_data, dir_raw, dir_res)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# IDs as listed in the handoff document; step 0b checks them against `ao`
exposure_ids <- c("ebi-a-GCST006867",  # T2D
                  "ieu-b-40",          # BMI
                  "ieu-b-38",          # systolic BP
                  "ieu-b-39",          # diastolic BP
                  "ieu-b-73",          # alcohol drinks/week
                  "ieu-b-4877")        # smoking initiation
outcome_ids  <- c("bbj-a-94",                  # cataract, Biobank Japan
                  "ukb-d-20002_1278",          # UKB self-report, 5,045 cases
                  "finn-b-H7_CATARACTSENILE",  # FinnGen R5 (primary European)
                  "ebi-a-GCST90018814")        # Sakaue 2021 European

# mr_raps is one of the methods; it needs the mr.raps package
if (!requireNamespace("mr.raps", quietly = TRUE)) {
  stop("Install mr.raps first:\n install.packages('mr.raps', repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))")
}

# ---- 0b. Check IDs against ao ------------------------------------------------
ao_file <- file.path(dir_data, "ao.rds")
if (!exists("ao")) {
  if (file.exists(ao_file)) ao <- readRDS(ao_file) else ao <- available_outcomes()
}
if (!file.exists(ao_file)) saveRDS(ao, ao_file)
all_ids <- c(exposure_ids, outcome_ids)
cols <- intersect(c("id", "trait", "population", "ncase", "sample_size"), names(ao))
print(ao[ao$id %in% all_ids, cols])
missing_ids <- setdiff(all_ids, ao$id)
if (length(missing_ids) > 0) stop("Not found in ao: ", paste(missing_ids, collapse = ", "))

# ---- Helper: try a server query, wait once, try again, then give up ----------
try_query <- function(fun, tries = 2, wait_min = 15) {
  for (k in seq_len(tries)) {
    res <- tryCatch(fun(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    message("   Attempt ", k, " failed: ", conditionMessage(res))
    if (k < tries) {
      message("   Waiting ", wait_min, " minutes before retrying...")
      Sys.sleep(wait_min * 60)
    }
  }
  res
}

# ---- 1. Exposure instruments (downloaded once, then read from disk) ----------
exp_file <- file.path(dir_data, "exposure_instruments.rds")
if (file.exists(exp_file)) {
  exp_dat <- readRDS(exp_file)
} else {
  exp_dat <- try_query(function() extract_instruments(outcomes = exposure_ids))
  if (inherits(exp_dat, "error")) stop("Instrument download failed. Try again later.")
  saveRDS(exp_dat, exp_file)
}
cat("\nInstruments per exposure:\n")
print(table(exp_dat$id.exposure))

# ---- 2. Outcome data, one pair at a time -------------------------------------
pairs <- expand.grid(exposure = exposure_ids, outcome = outcome_ids,
                     stringsAsFactors = FALSE)
piece_file <- function(ex, oc) file.path(dir_raw, paste0(ex, "__", oc, ".rds"))

for (i in seq_len(nrow(pairs))) {
  ex <- pairs$exposure[i]; oc <- pairs$outcome[i]
  f  <- piece_file(ex, oc)
  if (file.exists(f)) next
  message("\n== Outcome data ", i, " of ", nrow(pairs), ": ", ex, " -> ", oc, " ==")
  snps <- unique(exp_dat$SNP[exp_dat$id.exposure == ex])
  od <- try_query(function() extract_outcome_data(snps = snps, outcomes = oc))
  if (inherits(od, "error")) {
    message("\nServer not responding. Stopping downloads here.",
            "\nRe-run this script later; saved pairs will be skipped.")
    break
  }
  saveRDS(od, f)
  Sys.sleep(10)   # short pause between pairs to go easy on the server
}

pairs$downloaded <- file.exists(mapply(piece_file, pairs$exposure, pairs$outcome))
cat("\nDownload status:\n")
print(pairs)

# ---- 3. Harmonise ------------------------------------------------------------
dat_list <- list()
for (i in which(pairs$downloaded)) {
  od <- readRDS(piece_file(pairs$exposure[i], pairs$outcome[i]))
  if (is.null(od) || nrow(od) == 0) {
    message("No outcome rows for ", pairs$exposure[i], " -> ", pairs$outcome[i])
    next
  }
  ed <- exp_dat[exp_dat$id.exposure == pairs$exposure[i], ]
  dat_list[[length(dat_list) + 1]] <- harmonise_data(ed, od)
}
dat_all <- dplyr::bind_rows(dat_list)
saveRDS(dat_all, file.path(dir_res, "harmonised_all.rds"))
write.csv(dat_all, file.path(dir_res, "harmonised_all.csv"), row.names = FALSE)

# ---- 4. Every method, every pair ---------------------------------------------
methods <- mr_method_list()$obj
cat("\nMethods in mr_method_list():", length(methods), "\n")
print(mr_method_list()[, c("obj", "name")])

set.seed(693)  # bootstrap- and mode-based methods use random draws
combos   <- unique(dat_all[, c("id.exposure", "id.outcome")])
res_list <- list()
fail_list <- list()
for (j in seq_len(nrow(combos))) {
  d <- dat_all[dat_all$id.exposure == combos$id.exposure[j] &
               dat_all$id.outcome  == combos$id.outcome[j], ]
  n_keep <- sum(d$mr_keep)
  for (m in methods) {
    if (m == "mr_wald_ratio" && n_keep > 1) next  # Wald ratio is for 1-SNP only
    r <- tryCatch(mr(d, method_list = m), error = function(e) e)
    if (inherits(r, "error")) {
      fail_list[[length(fail_list) + 1]] <- data.frame(
        id.exposure = combos$id.exposure[j], id.outcome = combos$id.outcome[j],
        method = m, error = conditionMessage(r))
    } else {
      res_list[[length(res_list) + 1]] <- r
    }
  }
}
res_all  <- dplyr::bind_rows(res_list)
fail_all <- dplyr::bind_rows(fail_list)

# Sensitivity statistics computed on exactly the same harmonised data
het_all   <- mr_heterogeneity(dat_all)
pleio_all <- mr_pleiotropy_test(dat_all)

# ---- 5. Save -----------------------------------------------------------------
saveRDS(res_all, file.path(dir_res, "mr_all_methods.rds"))
write.csv(res_all,   file.path(dir_res, "mr_all_methods.csv"),   row.names = FALSE)
write.csv(fail_all,  file.path(dir_res, "mr_method_failures.csv"), row.names = FALSE)
write.csv(het_all,   file.path(dir_res, "heterogeneity_Q.csv"),  row.names = FALSE)
write.csv(pleio_all, file.path(dir_res, "egger_intercept.csv"),  row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(dir_res, "sessionInfo.txt"))

# ---- 6. Short printout to paste back -----------------------------------------
cat("\nPairs analysed:", nrow(combos), "of", nrow(pairs), "\n")
cat("Result rows:", nrow(res_all), " | Method failures:", nrow(fail_all), "\n")
if (nrow(fail_all) > 0) print(fail_all)
cat("\nIVW results:\n")
print(res_all[res_all$method == "Inverse variance weighted",
              c("exposure", "id.outcome", "nsnp", "b", "se", "pval")])
