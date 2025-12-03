############################################################
# Appendix C: Correlation matrix for FINAL empirical spec
# Variables: ln_gdp_pc, agri_gdp, polstab_l1, exposure_std, post
# Output: table_corr_final.tex (booktabs LaTeX; lower triangle with stars)
############################################################

# ---- Guard: expect 'dd' and 'outdir' from the summary script ----
if (!exists("dd") || !exists("outdir")) {
  stop("Please run the summary table script first so 'dd' and 'outdir' exist.")
}

vars_corr <- c("ln_gdp_pc","agri_gdp","polstab_l1","exposure_std","post")
labels    <- c("ln(GDP per capita)",
               "Agriculture share of GDP (%)",
               "Political Stability (lag, 0--100)",
               "Exposure (baseline, SD units)",
               "Post (t \\ge 2005)")

X <- dd[, vars_corr, drop = FALSE]

# ---- Correlations (pairwise complete) and p-values ----
cor_mat <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))
p_mat   <- matrix(NA_real_, nrow = length(vars_corr), ncol = length(vars_corr),
                  dimnames = list(vars_corr, vars_corr))

for (i in seq_along(vars_corr)) {
  for (j in seq_along(vars_corr)) {
    xi <- X[[vars_corr[i]]]
    xj <- X[[vars_corr[j]]]
    keep <- is.finite(xi) & is.finite(xj)
    if (sum(keep) >= 3) {
      p_mat[i, j] <- tryCatch(cor.test(xi[keep], xj[keep])$p.value,
                              error = function(e) NA_real_)
    } else {
      p_mat[i, j] <- NA_real_
    }
  }
}

# Star formatter
stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.01) "$^{***}$" else if (p < 0.05) "$^{**}$" else if (p < 0.10) "$^{*}$" else ""
}

# ---- Build a display matrix (lower triangle w/ stars; diag = 1.00) ----
disp <- matrix("", nrow = length(vars_corr), ncol = length(vars_corr))
for (i in seq_along(vars_corr)) {
  for (j in seq_along(vars_corr)) {
    if (i == j) {
      disp[i, j] <- sprintf("%.2f", 1.00)
    } else if (i > j) {
      val <- sprintf("%.2f", cor_mat[i, j])
      st  <- stars(p_mat[i, j])
      disp[i, j] <- paste0(val, st)
    } else {
      disp[i, j] <- ""  # upper triangle blank for readability
    }
  }
}
rownames(disp) <- labels
colnames(disp) <- labels

# ---- Write LaTeX table ----
tex_path <- file.path(outdir, "table_corr_final.tex")
con <- file(tex_path, open = "wt", encoding = "UTF-8")

writeLines(c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Pairwise Correlations: Variables Used in Final Empirical Framework}",
  "\\label{tab:corr_final}",
  paste0("\\begin{tabular}{l", paste(rep("c", length(vars_corr)), collapse=""), "}"),
  "\\toprule",
  paste(c(" ", labels), collapse = " & "), " \\\\",
  "\\midrule"
), con)

for (i in seq_len(nrow(disp))) {
  line <- paste(c(rownames(disp)[i], disp[i, ]), collapse = " & ")
  writeLines(paste0(line, " \\\\"), con)
}

writeLines(c(
  "\\bottomrule",
  "\\multicolumn{", as.character(length(vars_corr) + 1),
  "}{l}{\\footnotesize{Notes: Lower triangle shows Pearson correlations; upper triangle blank. ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.}} \\\\",
  "\\end{tabular}",
  "\\end{table}"
), con)

close(con)
cat("Saved LaTeX correlation table to:", tex_path, "\n")
