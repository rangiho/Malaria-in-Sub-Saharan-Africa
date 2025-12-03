############################################################
# Nets vs antimalarial treatment: which is more effective?
# - Builds panel from Working Data_v3 (1).xlsx
# - Nearest-year match (|Delta year| <= 3) for sparse series
# - Standardized outcomes and treatments
# - OLS with HC3 robust SEs and Wald test beta_nets = beta_treat
# Outputs (in outdir):
#   - treat_effects_nearest.csv
#   - treat_effects_nearest.tex
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(readr)
  library(broom)
  library(sandwich)
  library(lmtest)
})

# ---------- Paths ----------
outdir    <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path <- file.path(outdir, "Working Data_v3 (1).xlsx")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helper functions ----------
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")

norm <- function(x) {
  x %>% str_replace_all("\\s+", " ") %>% str_trim() %>% tolower()
}

zfun <- function(v) {
  if (all(is.na(v))) NA_real_ else as.numeric(scale(v))
}

# ---------- Series labels ----------
label_nets_exact  <- "Use of insecticide-treated bed nets (% of under-5 population)"
label_treat_exact <- "Children with fever receiving antimalarial drugs (% of children under age 5 with fever)"

# ---------- Scan sheets and build panel ----------
map_and_make_panel <- function(sheet_name) {
  df <- tryCatch(
    read_xlsx(xlsx_path, sheet = sheet_name, na = na_tokens, guess_max = 1e5),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)
  
  df <- df %>%
    remove_empty("cols") %>%
    select(where(~ !all(is.na(.)))) %>%
    select(-matches("^Unnamed"))
  
  need_cols <- c("Series Name","Country Name")
  if (!all(need_cols %in% names(df))) return(NULL)
  
  # Identify year columns like "2000" or "2000 [YR2000]"
  id_like   <- c("Series Name","Series Code","Country Name","Country Code")
  year_cols <- setdiff(names(df)[str_detect(names(df), "^(\\d{4})(\\s*\\[YR\\1\\])?$")], id_like)
  if (!length(year_cols)) return(NULL)
  
  sN <- norm(df$`Series Name`)
  
  is_nets  <- sN == norm(label_nets_exact) |
    str_detect(sN, "insecticide.*bed net|treated.*net.*under.?5")
  is_treat <- sN == norm(label_treat_exact) |
    str_detect(sN, "fever.*antimalarial|receiv.*antimalarial.*fever")
  
  var_compact <- rep(NA_character_, length(sN))
  
  var_compact[sN == norm("GDP per capita (current US$)")] <- "gdp_pc_cur"
  var_compact[sN == norm("GDP (current US$)")]            <- "gdp_current_usd"
  var_compact[sN == norm("Population, total")]            <- "pop"
  var_compact[sN == norm("School enrollment, primary (% gross)")] <- "enroll_primary_gross"
  var_compact[sN == norm("Primary completion rate, total (% of relevant age group)")] <- "primary_completion"
  var_compact[sN == norm("Educational attainment, at least completed primary, population 25+ years, total (%) (cumulative)")] <- "attain_primary_25p"
  
  var_compact[is_nets]  <- "nets_u5"
  var_compact[is_treat] <- "treat_u5"
  
  df$var_compact <- var_compact
  
  panel <- df %>%
    filter(!is.na(var_compact)) %>%
    pivot_longer(all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
    mutate(
      year  = readr::parse_number(year_raw),
      var   = var_compact,
      value = suppressWarnings(readr::parse_number(as.character(val_raw)))
    ) %>%
    filter(!is.na(var)) %>%
    transmute(country = `Country Name`,
              code    = `Country Code`,
              year, var, value) %>%
    group_by(country, code, year, var) %>%
    summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
    pivot_wider(id_cols = c(country, code, year),
                names_from = var, values_from = value) %>%
    arrange(country, year) %>%
    mutate(
      gdp_pc_cur = coalesce(gdp_pc_cur, gdp_current_usd / pop),
      ln_gdp_pc  = ifelse(is.finite(gdp_pc_cur) & gdp_pc_cur > 0, log(gdp_pc_cur), NA_real_)
    )
  
  cov <- tibble(
    nets  = sum(!is.na(panel$nets_u5)),
    treat = sum(!is.na(panel$treat_u5)),
    gdp   = sum(!is.na(panel$ln_gdp_pc)),
    edu   = sum(!is.na(panel$enroll_primary_gross) |
                  !is.na(panel$primary_completion) |
                  !is.na(panel$attain_primary_25p))
  )
  score <- as.numeric(cov$nets + cov$treat + cov$gdp + cov$edu)
  
  list(sheet = sheet_name, panel = panel, coverage = cov, score = score)
}

sheets  <- excel_sheets(xlsx_path)
scanned <- purrr::map(sheets, map_and_make_panel) %>% purrr::compact()
stopifnot("No usable sheets found" = length(scanned) > 0)

cov_df <- tibble(
  sheet = sapply(scanned, `[[`, "sheet"),
  nets  = sapply(scanned, function(x) x$coverage$nets),
  treat = sapply(scanned, function(x) x$coverage$treat),
  gdp   = sapply(scanned, function(x) x$coverage$gdp),
  edu   = sapply(scanned, function(x) x$coverage$edu),
  score = sapply(scanned, `[[`, "score")
) %>% arrange(desc(score))

cat("\n--- SHEET COVERAGE ---\n"); print(cov_df, n = Inf)

best  <- scanned[[which.max(sapply(scanned, `[[`, "score"))]]
panel <- best$panel
cat("\nChosen sheet:", best$sheet, "\n")

# ---------- Build standardized outcomes ----------
panel <- panel %>%
  mutate(
    z_lngdppc = zfun(ln_gdp_pc),
    z_enroll  = zfun(enroll_primary_gross),
    z_complete= zfun(primary_completion),
    z_attain  = zfun(attain_primary_25p)
  ) %>%
  mutate(
    edu_raw   = rowMeans(cbind(z_enroll, z_complete, z_attain), na.rm = TRUE),
    edu_index = ifelse(is.finite(edu_raw), zfun(edu_raw), NA_real_)
  )

# ---------- Nearest-year matching for GDP and education ----------
max_gap <- 3      # tolerance in years
by_keys <- "country"

nets_long <- panel %>% select(country, year, nets_u5)  %>% drop_na(nets_u5)  %>% rename(year_r = year)
trt_long  <- panel %>% select(country, year, treat_u5) %>% drop_na(treat_u5) %>% rename(year_r = year)

out_gdp <- panel %>% select(country, year, y = z_lngdppc)  %>% drop_na(y) %>% rename(year_l = year)
out_edu <- panel %>% select(country, year, y = edu_index)  %>% drop_na(y) %>% rename(year_l = year)

nearest_for <- function(out_df, right_df, right_val_name) {
  out_df %>%
    inner_join(right_df, by = by_keys, relationship = "many-to-many") %>%
    mutate(gap = abs(year_l - year_r)) %>%
    group_by(country, year_l) %>%
    slice_min(gap, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(gap <= max_gap) %>%
    transmute(country, year = year_l, y, !!right_val_name := .data[[right_val_name]])
}

gdp_nets <- nearest_for(out_gdp, nets_long, "nets_u5")
gdp_trt  <- nearest_for(out_gdp, trt_long,  "treat_u5")

edu_nets <- nearest_for(out_edu, nets_long, "nets_u5")
edu_trt  <- nearest_for(out_edu, trt_long,  "treat_u5")

gdp_s <- full_join(gdp_nets, gdp_trt, by = c("country","year","y")) %>%
  mutate(
    z_y     = zfun(y),
    z_nets  = zfun(nets_u5),
    z_treat = zfun(treat_u5)
  ) %>%
  filter(!is.na(z_y)) %>%
  filter(!(is.na(z_nets) & is.na(z_treat)))

edu_s <- full_join(edu_nets, edu_trt, by = c("country","year","y")) %>%
  mutate(
    z_y     = zfun(y),
    z_nets  = zfun(nets_u5),
    z_treat = zfun(treat_u5)
  ) %>%
  filter(!is.na(z_y)) %>%
  filter(!(is.na(z_nets) & is.na(z_treat)))

cat("Matched samples (rows = country x year): GDP =", nrow(gdp_s),
    " | Education =", nrow(edu_s), "\n")

# ---------- OLS with HC3 robust SEs ----------
safe_fit <- function(df) {
  if (nrow(df) < 5) return(NULL)
  keep <- c("z_nets","z_treat")[colSums(!is.na(df[, c("z_nets","z_treat")])) > 0]
  if (!length(keep)) return(NULL)
  fml <- as.formula(paste0("z_y ~ ", paste(keep, collapse = " + ")))
  m   <- lm(fml, data = df)
  vc  <- sandwich::vcovHC(m, type = "HC3")
  list(model = m, vcov = vc)
}

res_gdp <- safe_fit(gdp_s)
res_edu <- safe_fit(edu_s)

collect_rows <- function(res, outcome_label) {
  if (is.null(res)) return(tibble())
  broom::tidy(res$model, vcov = res$vcov, conf.int = FALSE) %>%
    filter(term %in% c("z_nets","z_treat")) %>%
    transmute(
      Outcome      = outcome_label,
      Intervention = recode(term,
                            z_nets = "Mosquito nets",
                            z_treat = "Antimalarial treatment"),
      Effect    = round(estimate, 3),
      SE        = round(std.error, 3),
      t         = round(statistic, 2),
      p         = signif(p.value, 3),
      Countries = nobs(res$model),
      Method    = sprintf("OLS (HC3), nearest-year |Delta year| <= %d", max_gap)
    )
}

results <- bind_rows(
  collect_rows(res_gdp, "GDP (z)"),
  collect_rows(res_edu, "Education Index (z)")
)

# ---------- Wald test beta_nets = beta_treat ----------
wald_equal <- function(res) {
  if (is.null(res)) return(NA_real_)
  b <- coef(res$model)
  if (!all(c("z_nets","z_treat") %in% names(b))) return(NA_real_)
  V <- res$vcov
  R <- c(z_nets = 1, z_treat = -1)
  diff <- sum(R[names(b)] * b)
  var  <- as.numeric(t(R[names(b)]) %*% V[names(b), names(b)] %*% R[names(b)])
  if (!is.finite(var) || var <= 0) return(NA_real_)
  2 * pnorm(-abs(diff / sqrt(var)))
}

p_wald_gdp <- wald_equal(res_gdp)
p_wald_edu <- wald_equal(res_edu)

cat("\n--- Nets vs treatment (main question) ---\n")
print(results)
cat("\nWald p-value beta_nets = beta_treat: GDP =", p_wald_gdp,
    " | Education =", p_wald_edu, "\n")

# ---------- Save CSV and LaTeX ----------
csv_path <- file.path(outdir, "treat_effects_nearest.csv")
write_csv(results, csv_path)

tex_path <- file.path(outdir, "treat_effects_nearest.tex")
con <- file(tex_path, open = "wt", encoding = "UTF-8")
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  sprintf("\\caption{Nets vs Antimalarial Treatment: nearest-year match (|$\\Delta$year|\\le %d), standardized}", max_gap),
  "\\label{tab:treat_nearest}",
  "\\small",
  "\\begin{tabular}{@{}lccccccc@{}}\\toprule",
  "Outcome & Intervention & Effect & SE & $t$ & $p$ & Countries & Method \\\\ \\midrule"
), con)
if (nrow(results) == 0) {
  writeLines(" \\multicolumn{8}{c}{No estimable effects.} \\\\", con)
} else {
  apply(results, 1, function(r) {
    writeLines(sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                       r[["Outcome"]], r[["Intervention"]], r[["Effect"]],
                       r[["SE"]], r[["t"]], r[["p"]], r[["Countries"]], r[["Method"]]), con)
  })
}
wald_line <- sprintf("Wald test (nets = treatment): GDP $p$ = %s; Education $p$ = %s.",
                     ifelse(is.na(p_wald_gdp),"NA",signif(p_wald_gdp,3)),
                     ifelse(is.na(p_wald_edu),"NA",signif(p_wald_edu,3)))
writeLines(c(
  "\\bottomrule",
  paste0("\\multicolumn{8}{p{0.95\\linewidth}}{\\footnotesize Outcomes and interventions are standardized (z-scores). ",
         "HC3 robust standard errors. ", wald_line, "}"),
  "\\end{tabular}\\end{table}"
), con)
close(con)

cat("\nSaved:\n- ", csv_path, "\n- ", tex_path, "\n", sep = "")
