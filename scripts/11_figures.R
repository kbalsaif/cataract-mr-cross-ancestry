# 11_figures.R
# BST 693 - figures for the paper draft, drawn from the saved results only
#
# What this does (no server contact; nothing is re-estimated):
#   Figure 1  Forest plot of the main estimates: T2D and BMI on each cataract
#             outcome, by method (IVW, weighted median, MR-Egger, weighted mode,
#             RAPS, MR-PRESSO outlier-corrected), plus the UK Biobank estimate
#             converted to the log-odds scale.
#   Figure 2  T2D on cataract: the East Asian estimate on a like-for-like scale
#             (five routes) next to the FinnGen estimate.
#   Figure 3  Instrument transferability: each SNP's effect on T2D (or BMI) in
#             Europeans against its effect in East Asians, with the IVW slope.
#   Figure 4  Multivariable MR: T2D and BMI, each alone and adjusted for the other.
#   Figure S1 Scatter plots for the four primary pairs (TwoSampleMR).
#   Figure S2 Funnel plots for the four primary pairs (TwoSampleMR).
# Each figure is written as PDF (for the knitted PDF) and PNG (for Word) to
# results/figures/. The numbers plotted are printed at the end so they can be
# checked against the results workbook.
#
# Run with:  source("~/Desktop/BST693/scripts/11_figures.R")

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
dir_fig  <- file.path(dir_res, "figures")
dir.create(dir_fig, showWarnings = FALSE)

need <- c(file.path(dir_res, "mr_all_methods.csv"),
          file.path(dir_res, "sensitivity", "mr_presso_summary.csv"),
          file.path(dir_res, "ukb_logodds_conversion.csv"),
          file.path(dir_res, "eas_arm_sensitivity", "eas_arm_sensitivity.csv"),
          file.path(dir_res, "eas_instruments", "eas_instruments_results.csv"),
          file.path(dir_res, "mvmr", "mvmr_results.csv"),
          file.path(dir_res, "transferability", "harmonised__bbj-a-153.rds"),
          file.path(dir_res, "transferability", "harmonised__bbj-a-1.rds"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) stop("Missing result files:\n  ", paste(miss, collapse = "\n  "))

res    <- read.csv(file.path(dir_res, "mr_all_methods.csv"), stringsAsFactors = FALSE)
presso <- read.csv(file.path(dir_res, "sensitivity", "mr_presso_summary.csv"), stringsAsFactors = FALSE)
ukb    <- read.csv(file.path(dir_res, "ukb_logodds_conversion.csv"), stringsAsFactors = FALSE)
eas    <- read.csv(file.path(dir_res, "eas_arm_sensitivity", "eas_arm_sensitivity.csv"), stringsAsFactors = FALSE)
easi   <- read.csv(file.path(dir_res, "eas_instruments", "eas_instruments_results.csv"), stringsAsFactors = FALSE)
mvmr   <- read.csv(file.path(dir_res, "mvmr", "mvmr_results.csv"), stringsAsFactors = FALSE)

T2D <- "ebi-a-GCST006867"; BMI <- "ieu-b-40"
BBJ <- "bbj-a-94"; UKB <- "ukb-d-20002_1278"; FINN <- "finn-b-H7_CATARACTSENILE"; SAK <- "ebi-a-GCST90018814"
oc_lab <- c("Biobank Japan (East Asian)", "FinnGen R5 (European)", "Sakaue 2021 (European)")
names(oc_lab) <- c(BBJ, FINN, SAK)
ivw <- "Inverse variance weighted"

# ---- helper: draw one figure to PDF and PNG --------------------------------------
save_fig <- function(name, draw, width, height) {
  pdf(file.path(dir_fig, paste0(name, ".pdf")), width = width, height = height); draw(); dev.off()
  png(file.path(dir_fig, paste0(name, ".png")), width = width, height = height, units = "in", res = 600); draw(); dev.off()
  cat("Wrote", name, "\n")
}

# ---- helper: forest plot (base graphics) -----------------------------------------
# tab: data frame with columns label, b, se, header (TRUE = group heading row),
#      pch, col. Rows are drawn top to bottom in the order given. Labels go in the
#      left margin, "OR (95% CI)" text in the right margin.
forest <- function(tab, xlab, xlim, ref = 1, title = "", mar_left = 16, mar_right = 8, cex_lab = 0.8) {
  tab$or <- exp(tab$b); tab$lo <- exp(tab$b - 1.96 * tab$se); tab$hi <- exp(tab$b + 1.96 * tab$se)
  n <- nrow(tab); y <- rev(seq_len(n))
  op <- par(mar = c(4, mar_left, if (nzchar(title)) 2.5 else 1, mar_right)); on.exit(par(op))
  plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), log = "x", xaxt = "n", yaxt = "n", ylab = "", xlab = xlab,
       bty = "n", main = title, font.main = 1, cex.main = 0.95)
  ticks <- pretty(xlim, n = 6); ticks <- ticks[ticks >= xlim[1] & ticks <= xlim[2]]
  axis(1, at = ticks, labels = format(ticks))
  abline(v = ref, lty = 3, col = "grey50")
  ok <- !tab$header & !is.na(tab$b)
  lo_c <- pmax(tab$lo, xlim[1]); hi_c <- pmin(tab$hi, xlim[2])
  segments(lo_c[ok], y[ok], hi_c[ok], y[ok], col = tab$col[ok], lwd = 1.3)
  points(tab$or[ok], y[ok], pch = tab$pch[ok], col = tab$col[ok], bg = "white", cex = 1.1)
  hd <- tab$header
  mtext(tab$label[hd],  side = 2, line = mar_left - 1, at = y[hd],  las = 1, adj = 0, cex = cex_lab, font = 2)
  mtext(tab$label[!hd], side = 2, line = mar_left - 1, at = y[!hd], las = 1, adj = 0, cex = cex_lab, font = 1)
  txt <- ifelse(ok, sprintf("%.2f (%.2f-%.2f)", tab$or, tab$lo, tab$hi), "")
  mtext(txt, side = 4, line = 0.5, at = y, las = 1, adj = 0, cex = cex_lab * 0.95)
}
row_hd <- function(label) data.frame(label = label, b = NA, se = NA, header = TRUE, pch = NA, col = NA)
row_pt <- function(label, b, se, pch = 15, col = "black") data.frame(label = paste0("   ", label), b = b, se = se, header = FALSE, pch = pch, col = col)

pick <- function(ex, oc, m) {
  r <- res[res$id.exposure == ex & res$id.outcome == oc & res$method == m, ]
  if (nrow(r) != 1) return(c(NA, NA))
  c(r$b, r$se)
}
methods_show <- c("Inverse variance weighted" = "IVW (random effects)",
                  "Weighted median" = "Weighted median",
                  "MR Egger" = "MR-Egger",
                  "Weighted mode" = "Weighted mode",
                  "Robust adjusted profile score (RAPS)" = "RAPS")

# ---- Figure 1: main estimates by method ---------------------------------------------
build_main <- function(ex) {
  rows <- list()
  for (oc in c(BBJ, FINN, SAK)) {
    rows[[length(rows) + 1]] <- row_hd(oc_lab[[oc]])
    for (m in names(methods_show)) { v <- pick(ex, oc, m); rows[[length(rows) + 1]] <- row_pt(methods_show[[m]], v[1], v[2]) }
    p <- presso[presso$pair == paste0(ex, "__", oc), ]
    if (nrow(p) == 1 && !is.na(p$corrected_b)) {
      rows[[length(rows) + 1]] <- row_pt(sprintf("MR-PRESSO (%d outlier%s removed)", p$n_outliers, ifelse(p$n_outliers == 1, "", "s")),
                                         p$corrected_b, p$corrected_se)
    }
  }
  rows[[length(rows) + 1]] <- row_hd("UK Biobank self-report (European)")
  u <- ukb[grepl(if (ex == T2D) "diabetes" else "body mass", ukb$exposure, ignore.case = TRUE), ]
  rows[[length(rows) + 1]] <- row_pt("IVW, linear beta converted to log-odds", u$b_logodds[1], u$se_logodds[1], pch = 22, col = "grey40")
  do.call(rbind, rows)
}
tab_t2d <- build_main(T2D); tab_bmi <- build_main(BMI)
fig1 <- function() {
  layout(matrix(1:2, nrow = 2))
  forest(tab_t2d, xlab = "OR for cataract per unit log-odds of T2D liability", title = "A. Type 2 diabetes", xlim = c(0.85, 1.45))
  forest(tab_bmi, xlab = "OR for cataract per SD of BMI", title = "B. Body mass index", xlim = c(0.80, 1.90))
}
save_fig("fig1_forest_main", fig1, width = 8.5, height = 10)

# ---- Figure 2: like-for-like T2D comparison ------------------------------------------
r_main  <- eas[eas$exposure == "T2D" & eas$outcome == BBJ  & grepl("^0", eas$analysis) & eas$method == ivw, ]
r_C     <- eas[eas$exposure == "T2D" & eas$outcome == BBJ  & grepl("^C", eas$analysis) & grepl("5e-0?8", eas$analysis) & eas$method == ivw, ]
r_D     <- eas[eas$exposure == "T2D" & eas$outcome == BBJ  & grepl("^D", eas$analysis) & eas$method == ivw, ]
r_Cfin  <- eas[eas$exposure == "T2D" & eas$outcome == FINN & grepl("^C", eas$analysis) & grepl("5e-0?8", eas$analysis) & eas$method == ivw, ]
r_bbj77 <- easi[easi$id.exposure == "bbj-a-77" & easi$method == ivw, ]
r_spr   <- easi[easi$id.exposure == "ebi-a-GCST010118" & easi$method == ivw, ]
fin     <- pick(T2D, FINN, ivw)
n_fin   <- res$nsnp[res$id.exposure == T2D & res$id.outcome == FINN & res$method == ivw]
stopifnot(nrow(r_main) == 1, nrow(r_C) == 1, nrow(r_D) == 1, nrow(r_Cfin) == 1, nrow(r_bbj77) == 1, nrow(r_spr) == 1)

tab2 <- rbind(
  row_hd("Biobank Japan (East Asian)"),
  row_pt(sprintf("European instruments, European scale (pre-specified; %d SNPs)", r_main$nsnp), r_main$b, r_main$se),
  row_pt(sprintf("European instruments with East Asian p < 5e-8 (%d SNPs)", r_C$nsnp), r_C$b, r_C$se),
  row_pt(sprintf("Same %d SNPs, East Asian SNP-T2D effects as denominators", r_D$nsnp), r_D$b, r_D$se),
  row_pt(sprintf("East Asian instruments, Biobank Japan T2D GWAS (%d SNPs)", r_bbj77$nsnp), r_bbj77$b, r_bbj77$se),
  row_pt(sprintf("East Asian instruments, Spracklen 2020 GWAS (%d SNPs)", r_spr$nsnp), r_spr$b, r_spr$se),
  row_hd("FinnGen R5 (European)"),
  row_pt(sprintf("European instruments (pre-specified; %d SNPs)", n_fin), fin[1], fin[2], pch = 22, col = "grey40"),
  row_pt(sprintf("European instruments with East Asian p < 5e-8 (%d SNPs)", r_Cfin$nsnp), r_Cfin$b, r_Cfin$se, pch = 22, col = "grey40"))
fig2 <- function() {
  forest(tab2, xlab = "OR for cataract per unit log-odds of T2D liability (IVW)", xlim = c(1.0, 1.30), mar_left = 25, cex_lab = 0.78)
}
save_fig("fig2_like_for_like_t2d", fig2, width = 10.5, height = 4)

# ---- Figure 3: transferability scatter -------------------------------------------------
tr_t2d <- readRDS(file.path(dir_res, "transferability", "harmonised__bbj-a-153.rds"))
tr_bmi <- readRDS(file.path(dir_res, "transferability", "harmonised__bbj-a-1.rds"))
tr_t2d <- tr_t2d[tr_t2d$mr_keep, ]; tr_bmi <- tr_bmi[tr_bmi$mr_keep, ]
ivw_slope <- function(d) {                       # IVW slope through the origin (same point estimate as mr_ivw)
  w <- 1 / d$se.outcome^2
  sum(w * d$beta.exposure * d$beta.outcome) / sum(w * d$beta.exposure^2)
}
s_t2d <- ivw_slope(tr_t2d); s_bmi <- ivw_slope(tr_bmi)
scatter_tr <- function(d, slope, main, xlab, ylab) {
  lim <- range(c(d$beta.exposure, d$beta.outcome, 0)) * 1.05
  plot(d$beta.exposure, d$beta.outcome, xlim = lim, ylim = lim, type = "n",
       xlab = xlab, ylab = ylab, main = main, font.main = 1, cex.main = 0.95, bty = "l")
  abline(h = 0, v = 0, col = "grey85", lty = 3)
  segments(d$beta.exposure, d$beta.outcome - 1.96 * d$se.outcome, d$beta.exposure, d$beta.outcome + 1.96 * d$se.outcome, col = "grey75", lwd = 0.5)
  points(d$beta.exposure, d$beta.outcome, pch = 16, cex = 0.6, col = "grey30")
  abline(0, 1, lty = 2, col = "grey50"); abline(0, slope, lwd = 1.6)
  legend("topleft", bty = "n", cex = 0.8, lty = c(2, 1), lwd = c(1, 1.6), col = c("grey50", "black"),
         legend = c("Slope 1 (same effect in both ancestries)", sprintf("IVW slope = %.2f", slope)))
}
fig3 <- function() {
  layout(matrix(1:2, nrow = 1)); op <- par(mar = c(4.2, 4.2, 2.5, 0.8))
  scatter_tr(tr_t2d, s_t2d, sprintf("A. T2D: %d European instruments in Biobank Japan", nrow(tr_t2d)),
             "SNP effect on T2D, European GWAS (log-odds)", "SNP effect on T2D, Biobank Japan (log-odds)")
  scatter_tr(tr_bmi, s_bmi, sprintf("B. BMI: %d European instruments in Biobank Japan", nrow(tr_bmi)),
             "SNP effect on BMI, European GWAS (SD units)", "SNP effect on BMI, Biobank Japan (SD units)")
  par(op)
}
save_fig("fig3_transferability", fig3, width = 10, height = 5)

# ---- Figure 4: multivariable MR --------------------------------------------------------
mv_pick <- function(ex, oc) {
  if ("id.exposure" %in% names(mvmr)) {
    r <- mvmr[mvmr$id.exposure == ex & mvmr$outcome_id == oc, ]
  } else {
    r <- mvmr[grepl(ex, mvmr$exposure, fixed = TRUE) & mvmr$outcome_id == oc, ]
  }
  if (nrow(r) != 1) return(c(NA, NA))
  c(r$b, r$se)
}
tab4 <- list()
for (ex in c(T2D, BMI)) {
  exl <- if (ex == T2D) "Type 2 diabetes" else "Body mass index"
  adj <- if (ex == T2D) "adjusted for BMI" else "adjusted for T2D"
  tab4[[length(tab4) + 1]] <- row_hd(exl)
  for (oc in c(BBJ, FINN, SAK)) {
    u <- pick(ex, oc, ivw); m <- mv_pick(ex, oc)
    tab4[[length(tab4) + 1]] <- row_pt(paste0(oc_lab[[oc]], ", univariable"), u[1], u[2], pch = 22, col = "grey40")
    tab4[[length(tab4) + 1]] <- row_pt(paste0(oc_lab[[oc]], ", ", adj), m[1], m[2])
  }
}
tab4 <- do.call(rbind, tab4)
fig4 <- function() {
  forest(tab4, xlab = "OR for cataract (IVW; per unit log-odds of T2D liability, per SD of BMI)", xlim = c(0.90, 1.35), mar_left = 19, cex_lab = 0.78)
}
save_fig("fig4_mvmr", fig4, width = 9.5, height = 5.2)

# ---- Supplementary: scatter and funnel plots for the four primary pairs (TwoSampleMR) --
ok_supp <- tryCatch({
  suppressPackageStartupMessages({ library(TwoSampleMR); library(ggplot2) })
  dat_all <- readRDS(file.path(dir_res, "harmonised_all.rds"))
  res_rds <- readRDS(file.path(dir_res, "mr_all_methods.rds"))
  single  <- read.csv(file.path(dir_res, "sensitivity", "single_snp.csv"), stringsAsFactors = FALSE)
  prim <- expand.grid(ex = c(T2D, BMI), oc = c(BBJ, FINN), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(prim))) {
    ex <- prim$ex[i]; oc <- prim$oc[i]
    d  <- dat_all[dat_all$id.exposure == ex & dat_all$id.outcome == oc, ]
    r  <- res_rds[res_rds$id.exposure == ex & res_rds$id.outcome == oc &
                  res_rds$method %in% c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode"), ]
    tag <- paste0(if (ex == T2D) "T2D" else "BMI", "_", if (oc == BBJ) "BBJ" else "FinnGen")
    p1 <- mr_scatter_plot(r, d)[[1]] + theme_bw(base_size = 10) + theme(legend.position = "bottom") +
      labs(x = sprintf("SNP effect on %s", if (ex == T2D) "T2D (log-odds)" else "BMI (SD)"),
           y = sprintf("SNP effect on cataract, %s (log-odds)", if (oc == BBJ) "Biobank Japan" else "FinnGen R5"))
    ggsave(file.path(dir_fig, paste0("figS1_scatter_", tag, ".pdf")), p1, width = 6, height = 5.5)
    ggsave(file.path(dir_fig, paste0("figS1_scatter_", tag, ".png")), p1, width = 6, height = 5.5, dpi = 300)
    s  <- single[single$id.exposure == ex & single$id.outcome == oc, ]
    p2 <- mr_funnel_plot(s)[[1]] + theme_bw(base_size = 10) + theme(legend.position = "bottom")
    ggsave(file.path(dir_fig, paste0("figS2_funnel_", tag, ".pdf")), p2, width = 6, height = 5.5)
    ggsave(file.path(dir_fig, paste0("figS2_funnel_", tag, ".png")), p2, width = 6, height = 5.5, dpi = 300)
    cat("Wrote figS1/figS2 for", tag, "\n")
  }
  TRUE
}, error = function(e) { message("Supplementary scatter/funnel plots skipped: ", conditionMessage(e)); FALSE })

# ---- Printout: the numbers that were plotted -------------------------------------------
fmt <- function(tab) {
  t <- tab[!tab$header, c("label", "b", "se")]
  t$OR <- round(exp(t$b), 3); t$lo <- round(exp(t$b - 1.96 * t$se), 3); t$hi <- round(exp(t$b + 1.96 * t$se), 3)
  t$b <- round(t$b, 4); t$se <- round(t$se, 4); t
}
cat("\nFigure 1A (T2D):\n"); print(fmt(tab_t2d), row.names = FALSE)
cat("\nFigure 1B (BMI):\n"); print(fmt(tab_bmi), row.names = FALSE)
cat("\nFigure 2:\n"); print(fmt(tab2), row.names = FALSE)
cat("\nFigure 3 IVW slopes (should match transferability_results.csv): T2D", round(s_t2d, 3), " BMI", round(s_bmi, 3), "\n")
cat("\nFigure 4:\n"); print(fmt(tab4), row.names = FALSE)
cat("\nFiles in", dir_fig, ":\n"); print(list.files(dir_fig))
