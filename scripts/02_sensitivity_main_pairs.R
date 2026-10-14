# 02_sensitivity_main_pairs.R
# BST 693 - sensitivity analyses for the 8 main pairs (T2D and BMI x 4 outcomes)
#
# What this does, for each main pair, on the SAME harmonised data as script 01:
#   1. MR-PRESSO: looks for individual SNPs that behave as outliers, removes
#      them, and re-estimates. Also tests whether removing them changes the answer.
#   2. Leave-one-out: re-runs IVW dropping one SNP at a time, to see whether
#      any single SNP is driving the result.
#   3. Single-SNP estimates, used for the funnel plot.
#   4. Scatter, funnel and leave-one-out plots saved as PDF.
# Everything is saved to results/sensitivity/.

library(TwoSampleMR)

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_sens <- file.path(dir_res, "sensitivity")
dir_plot <- file.path(dir_sens, "plots")
for (d in c(dir_sens, dir_plot)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dat_all <- readRDS(file.path(dir_res, "harmonised_all.rds"))
res_all <- readRDS(file.path(dir_res, "mr_all_methods.rds"))

main_exposures <- c("ebi-a-GCST006867", "ieu-b-40")   # T2D, BMI (IDs from handoff)
dat <- dat_all[dat_all$id.exposure %in% main_exposures, ]
combos <- unique(dat[, c("id.exposure", "id.outcome")])
cat("Main pairs:", nrow(combos), "\n")

set.seed(693)

# ---- 1. MR-PRESSO ------------------------------------------------------------
# The outlier test needs NbDistribution >= (number of SNPs / 0.05), otherwise
# MR-PRESSO warns the outlier test is unstable. So it is set per pair.
# T2D pairs (~115 SNPs) run first and are fast; BMI pairs (~480 SNPs) are slow.
# Each pair is saved as soon as it finishes; re-running skips finished pairs.
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  stop("Install MR-PRESSO first:\n install.packages('remotes'); remotes::install_github('rondolab/MR-PRESSO')")
}
dir_presso <- file.path(dir_sens, "presso_pieces")
dir.create(dir_presso, showWarnings = FALSE)
combos <- combos[order(combos$id.exposure != "ebi-a-GCST006867"), ]  # T2D first

presso_list <- list()
for (i in seq_len(nrow(combos))) {
  nm <- paste(combos$id.exposure[i], combos$id.outcome[i], sep = "__")
  f  <- file.path(dir_presso, paste0(nm, ".rds"))
  if (!file.exists(f)) {
    d <- dat[dat$id.exposure == combos$id.exposure[i] & dat$id.outcome == combos$id.outcome[i], ]
    n_snp <- sum(d$mr_keep)
    nb <- max(3000, ceiling(n_snp / 0.05 / 1000) * 1000 + 1000)
    message("\nMR-PRESSO ", i, " of ", nrow(combos), ": ", nm,
            " (", n_snp, " SNPs, ", nb, " simulations) started ", format(Sys.time(), "%H:%M"))
    set.seed(693)
    p <- tryCatch(run_mr_presso(d, NbDistribution = nb, SignifThreshold = 0.05),
                  error = function(e) e)
    if (inherits(p, "error")) {
      message("  MR-PRESSO failed: ", conditionMessage(p))
      out <- list(error = conditionMessage(p))
    } else {
      out <- p[[1]]
    }
    saveRDS(out, f)
  }
  presso_list[[nm]] <- readRDS(f)
}
saveRDS(presso_list, file.path(dir_sens, "mr_presso_raw.rds"))

# Pull the key numbers into one table
presso_tab <- do.call(rbind, lapply(names(presso_list), function(nm) {
  p <- presso_list[[nm]]
  if (!is.null(p$error)) {
    return(data.frame(pair = nm, raw_b = NA, raw_se = NA, raw_p = NA,
                      corrected_b = NA, corrected_se = NA, corrected_p = NA,
                      global_p = NA, n_outliers = NA, distortion_p = NA, note = p$error))
  }
  mr <- p$`Main MR results`
  gl <- p$`MR-PRESSO results`$`Global Test`$Pvalue
  out <- p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  dp  <- p$`MR-PRESSO results`$`Distortion Test`$Pvalue
  data.frame(pair = nm,
             raw_b = mr$`Causal Estimate`[1], raw_se = mr$Sd[1], raw_p = mr$`P-value`[1],
             corrected_b = mr$`Causal Estimate`[2], corrected_se = mr$Sd[2], corrected_p = mr$`P-value`[2],
             global_p = if (is.null(gl)) NA else gl,
             n_outliers = if (is.null(out) || is.character(out)) 0 else length(out),
             distortion_p = if (is.null(dp)) NA else dp,
             note = "")
}))
write.csv(presso_tab, file.path(dir_sens, "mr_presso_summary.csv"), row.names = FALSE)

# ---- 2. Leave-one-out --------------------------------------------------------
loo <- mr_leaveoneout(dat)
write.csv(loo, file.path(dir_sens, "leave_one_out.csv"), row.names = FALSE)

# Summary: for each pair, the weakest result when any single SNP is dropped
loo_sum <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
  x <- loo[loo$id.exposure == combos$id.exposure[i] & loo$id.outcome == combos$id.outcome[i] & loo$SNP != "All", ]
  data.frame(id.exposure = combos$id.exposure[i], id.outcome = combos$id.outcome[i],
             min_b = min(x$b), max_b = max(x$b), max_p = max(x$p),
             snp_at_max_p = x$SNP[which.max(x$p)])
}))
write.csv(loo_sum, file.path(dir_sens, "leave_one_out_summary.csv"), row.names = FALSE)

# ---- 3. Single-SNP and plots -------------------------------------------------
single <- mr_singlesnp(dat)
write.csv(single, file.path(dir_sens, "single_snp.csv"), row.names = FALSE)

main_res <- res_all[res_all$id.exposure %in% main_exposures &
                    res_all$method %in% c("Inverse variance weighted", "MR Egger",
                                          "Weighted median", "Weighted mode"), ]
p_scatter <- mr_scatter_plot(main_res, dat)
p_funnel  <- mr_funnel_plot(single)
p_loo     <- mr_leaveoneout_plot(loo)

pdf(file.path(dir_plot, "scatter_plots.pdf"), width = 7, height = 6)
for (p in p_scatter) print(p)
dev.off()
pdf(file.path(dir_plot, "funnel_plots.pdf"), width = 7, height = 6)
for (p in p_funnel) print(p)
dev.off()
pdf(file.path(dir_plot, "leave_one_out_plots.pdf"), width = 7, height = 30)
for (p in p_loo) print(p)
dev.off()

# ---- 4. Printout to paste back -----------------------------------------------
cat("\nMR-PRESSO summary:\n")
print(presso_tab[, c("pair", "raw_b", "raw_p", "corrected_b", "corrected_p", "global_p", "n_outliers", "distortion_p", "note")],
      row.names = FALSE, digits = 3)
cat("\nLeave-one-out summary:\n")
print(loo_sum, row.names = FALSE, digits = 3)
