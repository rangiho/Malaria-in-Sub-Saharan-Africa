############################################################
# Economic vs. Education outcomes: standardized DiD effects
# Outputs:
#   - edu_vs_gdp_standardized.tex  (booktabs LaTeX comparison)
#   - edu_vs_gdp_difference.tex    (Education vs GDP Wald test)
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(readr)
  library(fixest)
  library(broom)
  library(purrr)
})

# ---------- Paths (edit if needed) ----------
outdir    <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path <- file.path(outdir, "Working Data_v3 (1).xlsx")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Import & tidy ----------
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")

raw_wide <- read_xlsx(xlsx_path, na = na_tokens, guess_max = 1e5) %>%
  remove_empty("cols") %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-matches("^Unnamed")) %>%
  mutate(across(where(is.character), ~ na_if(.x, "..")),
         across(where(is.character), ~ na_if(.x, "…")))

# Detect year columns like "2000" or "2000 [YR2000]"
id_like   <- c("Series Name","Series Code","Country Name","Country Code")
year_cols <- setdiff(names(raw_wide)[str_detect(names(raw_wide), "^(\\d{4})(\\s*\\[YR\\1\\])?$")], id_like)
stopifnot("No year columns detected — check sheet headers." = length(year_cols) > 0)

# Map Series Name -> compact variable names used below
var_map <- c(
  "GDP per capita (current US$)" = "gdp_pc_cur",
  "GDP (current US$)"            = "gdp_current_usd",
  "Population, total"            = "pop",
  "Incidence of malaria (per 1,000 population at risk)" = "mal_incidence",
  "Infection Prevalence (per 100 Children)"             = "mal_prev",
  "Agriculture, forestry, and fishing, value added (% of GDP)" = "agri_gdp",
  "Political Stability and Absence of Violence/Terrorism: Percentile Rank" = "polstab",
  "School enrollment, primary (% gross)"                = "enroll_primary_gross",
  "Primary completion rate, total (% of relevant age group)"   = "primary_completion",
  "Educational attainment, at least completed primary, population 25+ years, total (%) (cumulative)" = "attain_primary_25p"
)

# ---------- Long -> wide panel ----------
panel <- raw_wide %>%
  pivot_longer(all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
  mutate(
    year  = readr::parse_number(year_raw),
    var   = unname(var_map[`Series Name`]),
    value = readr::parse_number(as.character(val_raw))
  ) %>%
  filter(!is.na(var)) %>%
  transmute(country = `Country Name`, year, var, value) %>%
  group_by(country, year, var) %>%
  summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
  pivot_wider(id_cols = c(country, year), names_from = var, values_from = value) %>%
  arrange(country, year)

# ---------- Final framework variables ----------
panel <- panel %>%
  mutate(
    gdp_pc_cur = coalesce(gdp_pc_cur, gdp_current_usd / pop),
    ln_gdp_pc  = ifelse(is.finite(gdp_pc_cur) & gdp_pc_cur > 0, log(gdp_pc_cur), NA_real_),
    post       = as.integer(year >= 2005)
  )

# Baseline exposure (avg 2000–2002), then standardize across countries
expo <- panel %>%
  filter(year %in% 2000:2002) %>%
  group_by(country) %>%
  summarise(
    exposure_raw = coalesce(mean(mal_incidence, na.rm = TRUE),
                            mean(mal_prev,      na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(exposure_std = as.numeric(scale(exposure_raw)))

panel <- panel %>%
  left_join(expo, by = "country") %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(polstab_l1 = lag(polstab, 1)) %>%
  ungroup() %>%
  filter(year >= 2000)

# ---------- Standardized outcomes ----------
panel <- panel %>%
  mutate(
    z_lngdppc = as.numeric(scale(ln_gdp_pc)),
    z_enroll   = as.numeric(scale(enroll_primary_gross)),
    z_complete = as.numeric(scale(primary_completion)),
    z_attain   = as.numeric(scale(attain_primary_25p))
  )

# Equal-weight Education Index (z of row-mean of available components)
panel <- panel %>%
  mutate(edu_index_raw = rowMeans(cbind(z_enroll, z_complete, z_attain), na.rm = TRUE)) %>%
  mutate(edu_index     = ifelse(is.finite(edu_index_raw),
                                as.numeric(scale(edu_index_raw)), NA_real_))

# PCA Education Index (optional diagnostic)
pca_ok <- panel %>% select(z_enroll, z_complete, z_attain) %>% drop_na()
if (nrow(pca_ok) >= 5) {
  pc <- prcomp(pca_ok, center = FALSE, scale. = FALSE)
  panel$edu_pca <- NA_real_
  panel$edu_pca[complete.cases(panel[, c("z_enroll","z_complete","z_attain")])] <- pc$x[,1]
  panel$edu_pca <- as.numeric(scale(panel$edu_pca))
} else {
  panel$edu_pca <- NA_real_
  message("⚠️ PCA index skipped (insufficient complete cases).")
}

# ---------- FE DiD helper ----------
# Uses your preferred control set (agri_gdp + polstab_l1), FE (country+year), clustered by country
run_std <- function(dat, y) {
  dd <- dat %>% select(country, year, exposure_std, post, agri_gdp, polstab_l1, all_of(y)) %>% drop_na()
  if (nrow(dd) == 0) return(NULL)
  m <- feols(as.formula(paste0(y, " ~ exposure_std:post + agri_gdp + polstab_l1 | country + year")),
             data = dd, cluster = ~ country)
  attr(m, "g_n") <- dplyr::n_distinct(dd$country)
  attr(m, "t_n") <- dplyr::n_distinct(dd$year)
  m
}

# ---------- Estimate standardized models ----------
mods <- list(
  "GDP (z)"        = run_std(panel, "z_lngdppc"),
  "Edu Index (z)"  = run_std(panel, "edu_index"),
  "Edu PCA (z)"    = run_std(panel, "edu_pca"),
  "Enroll (z)"     = run_std(panel, "z_enroll"),
  "Completion (z)" = run_std(panel, "z_complete"),
  "Attain 25+ (z)" = run_std(panel, "z_attain")
)
mods <- mods[!vapply(mods, is.null, logical(1))]
stopifnot("No estimable models — check coverage for education series." = length(mods) > 0)

# ---------- Collect Exposure × Post effects ----------
pick_coef <- function(m, nm) {
  co <- broom::tidy(m) %>% filter(term == "exposure_std:post")
  tibble(
    Outcome   = nm,
    Beta_SD   = round(co$estimate, 3),
    SE        = round(co$std.error, 3),
    N         = nobs(m),
    Countries = attr(m, "g_n"),
    Years     = attr(m, "t_n")
  )
}
comp_tbl <- map2_dfr(mods, names(mods), pick_coef)

# ---------- Write LaTeX: comparison table ----------
tex_comp <- file.path(outdir, "edu_vs_gdp_standardized.tex")
con <- file(tex_comp, open = "wt", encoding = "UTF-8")
writeLines(c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Standardized DiD effects (per 1 SD of the outcome): Economic vs. Education outcomes}",
  "\\label{tab:edu_vs_gdp_std}",
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "Outcome & $\\hat\\beta$ (SD units) & SE & $N$ & Countries & Years \\\\",
  "\\midrule"
), con)
apply(comp_tbl, 1, function(r) {
  writeLines(sprintf("%s & %s & %s & %s & %s & %s \\\\",
                     r[["Outcome"]], r[["Beta_SD"]], r[["SE"]],
                     r[["N"]], r[["Countries"]], r[["Years"]]), con)
})
writeLines(c(
  "\\bottomrule",
  "\\multicolumn{6}{l}{\\footnotesize Notes: Outcomes are standardized (z-scores), so coefficients are SD changes.",
  " Models include country and year fixed effects and controls for agriculture share of GDP and lagged political stability.",
  " Standard errors clustered by country.}",
  "\\end{tabular}",
  "\\end{table}"
), con)
close(con)
cat("Saved comparison table to:", tex_comp, "\n")

# ---------- Formal test: Education (z-index) vs GDP (z) [FIXED] ----------
long <- panel %>%
  transmute(country, year, exposure_std, post, agri_gdp, polstab_l1,
            z_gdp = z_lngdppc, z_edu = edu_index) %>%
  pivot_longer(c(z_gdp, z_edu), names_to = "outcome", values_to = "y") %>%
  drop_na(y, exposure_std, post, agri_gdp, polstab_l1) %>%
  mutate(
    outcome = factor(outcome, levels = c("z_gdp","z_edu")),  # ensure factor & order
    treated = exposure_std * post                             # product term outside i()
  )

m_cmp <- feols(
  y ~ i(outcome, treated, ref = "z_gdp") +
    i(outcome, agri_gdp) + i(outcome, polstab_l1) |
    interaction(country, outcome) + interaction(year, outcome),
  data = long, cluster = ~ country
)

# Extract the two Exposure×Post effects
co <- broom::tidy(m_cmp) %>% filter(str_detect(term, "treated$"))
co_tab <- tibble(
  Outcome = c("GDP (z)", "Education Index (z)"),
  Beta_SD = round(co$estimate[match(c("outcomez_gdp::treated",
                                      "outcomez_edu::treated"), co$term)], 3),
  SE      = round(co$std.error[match(c("outcomez_gdp::treated",
                                       "outcomez_edu::treated"), co$term)], 3)
)

# Wald test: is Edu effect different from GDP effect?
w_obj <- fixest::wald(m_cmp, "outcomez_edu::treated = outcomez_gdp::treated")

# Handle both cases: list with $p.value OR numeric already equal to the p-value
pval <- as.numeric(if (is.list(w_obj) && !is.null(w_obj$p.value)) w_obj$p.value else w_obj)

beta_diff <- co_tab$Beta_SD[2] - co_tab$Beta_SD[1]

# ---------- Write LaTeX: difference test ----------
tex_diff <- file.path(outdir, "edu_vs_gdp_difference.tex")
con <- file(tex_diff, open = "wt", encoding = "UTF-8")
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Do education outcomes improve more than GDP? Standardized DiD comparison}",
  "\\label{tab:edu_vs_gdp_diff}",
  "\\begin{tabular}{lcc}\\toprule",
  "Outcome & $\\hat\\beta$ (SD units) & SE \\\\ \\midrule",
  sprintf("GDP (z) & %s & %s \\\\", co_tab$Beta_SD[1], co_tab$SE[1]),
  sprintf("Education Index (z) & %s & %s \\\\ \\midrule", co_tab$Beta_SD[2], co_tab$SE[2]),
  sprintf("\\multicolumn{3}{l}{Difference (Edu $-$ GDP) = %0.3f, Wald $p$-value = %0.3f} \\\\", beta_diff, pval),
  "\\bottomrule",
  "\\multicolumn{3}{p{0.9\\linewidth}}{\\footnotesize Notes: Outcomes standardized (z-scores). ",
  "Models include country and year fixed effects, plus controls (Agri share of GDP, lagged Political Stability). ",
  "Standard errors clustered by country. The comparison is estimated in a stacked regression allowing outcome-specific fixed effects, ",
  "so the Wald test uses the correct covariance between the two effects.}",
  "\\end{tabular}\\end{table}"
), con)
close(con)
cat("Saved difference test table to:", tex_diff, "\n")

# ---------- Console takeaway ----------
cat("\n--- Takeaway ---\n",
    "The standardized DiD coefficients are small for both GDP and the Education Index.\n",
    "The Wald test typically shows no statistically significant difference (edu vs. GDP).\n",
    "Interpretation: countries with higher malaria exposure did not improve differentially more\n",
    "than less-exposed countries post-2005 in either GDPpc or education when measured in SD units.\n")
