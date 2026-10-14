# 04_cross_ancestry_test.R
# BST 693 - formal comparison of East Asian vs European estimates
#
# What this does:
#   For T2D and BMI, compares the Biobank Japan estimate with each European
#   estimate on the same (log-odds) scale. Two equivalent tests are reported:
#   a z-test on the difference, and Cochran's Q for the two estimates (Q = z^2
#   with 1 degree of freedom). I^2 describes how much of the spread between the
#   two estimates is beyond chance. Done for every method that has a usable SE.
#   The UKB self-report outcome is excluded: its effects are on a linear scale
#   and cannot be compared with log-odds estimates.

base_dir <- "~/Desktop/BST693"
dir_res  <- file.path(base_dir, "results")
res <- read.csv(file.path(dir_res, "mr_all_methods.csv"), stringsAsFactors = FALSE)

main_exposures <- c("ebi-a-GCST006867", "ieu-b-40")
eas <- "bbj-a-94"
eur <- c("finn-b-H7_CATARACTSENILE", "ebi-a-GCST90018814")

res <- res[res$id.exposure %in% main_exposures & !is.na(res$se) &
           res$method != "Unweighted regression", ]   # no usable SE for that one

out <- list()
for (ex in main_exposures) for (oc in eur) for (m in unique(res$method)) {
  a <- res[res$id.exposure == ex & res$id.outcome == eas & res$method == m, ]
  b <- res[res$id.exposure == ex & res$id.outcome == oc  & res$method == m, ]
  if (nrow(a) != 1 || nrow(b) != 1) next
  diff <- a$b - b$b
  se_d <- sqrt(a$se^2 + b$se^2)
  z    <- diff / se_d
  w    <- c(1 / a$se^2, 1 / b$se^2)
  est  <- c(a$b, b$b)
  Q    <- sum(w * (est - sum(w * est) / sum(w))^2)
  out[[length(out) + 1]] <- data.frame(
    exposure = sub(" \\|\\|.*", "", a$exposure), method = m, european_outcome = oc,
    b_eas = a$b, se_eas = a$se, b_eur = b$b, se_eur = b$se,
    diff = diff, se_diff = se_d, z = z, p_diff = 2 * pnorm(-abs(z)),
    Q = Q, Q_p = pchisq(Q, df = 1, lower.tail = FALSE),
    I2_pct = 100 * max(0, (Q - 1) / Q),
    OR_ratio_eas_vs_eur = exp(diff),
    stringsAsFactors = FALSE)
}
ca <- do.call(rbind, out)
write.csv(ca, file.path(dir_res, "cross_ancestry_test.csv"), row.names = FALSE)

cat("IVW: East Asian vs European (main comparison)\n")
print(ca[ca$method == "Inverse variance weighted",
         c("exposure", "european_outcome", "b_eas", "b_eur", "diff", "se_diff", "z", "p_diff", "I2_pct")],
      row.names = FALSE, digits = 3)

cat("\nAll methods: p-value for an ancestry difference\n")
tab <- ca[, c("exposure", "european_outcome", "method", "p_diff")]
tab$p_diff <- signif(tab$p_diff, 2)
print(tab[order(tab$exposure, tab$european_outcome, tab$p_diff), ], row.names = FALSE)

cat("\nSmallest p across all methods:", signif(min(ca$p_diff), 2),
    " | Number of comparisons:", nrow(ca), "\n")
