############################################################
# Malaria → SSA: Baseline DiD + Robustness & LaTeX outputs
# - Reads your Excel (SSA only)
# - Robust mapping from Series Name -> variables
# - Builds exposure_std (2000–2002; incidence then prevalence fallback)
# - Models: Baseline, +Safe, +All (mediators robustness), +Country trends
# - Event-study pre-trend test
# - Education outcomes (if present)
# - Exposure sensitivity (inc-only, prev-only, 2001–2003 window)
# - LaTeX tables/snippets written to outdir
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(fixest)
  library(texreg)
  library(readr)
})

# -------- Paths (edit if needed) --------
outdir     <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path  <- file.path(outdir, "Working Data_v3 (1).xlsx")
sheet_name <- NULL  # set to "SSA - Combined by Countries" if not the first sheet

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------- Import (treat '..', etc. as NA) --------
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")

raw_wide <- if (is.null(sheet_name)) {
  read_xlsx(xlsx_path, na = na_tokens, guess_max = 1e5)
} else {
  read_xlsx(xlsx_path, sheet = sheet_name, na = na_tokens, guess_max = 1e5)
}

raw_wide <- raw_wide %>%
  remove_empty("cols") %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-matches("^Unnamed")) %>%
  mutate(across(where(is.character), ~ na_if(.x, ".."))) %>%
  mutate(across(where(is.character), ~ na_if(.x, "…")))

# Detect year columns: "2000" or "2000 [YR2000]"
id_like   <- c("Series Name","Series Code","Country Name","Country Code")
year_cols <- setdiff(
  names(raw_wide)[str_detect(names(raw_wide), "^(\\d{4})(\\s*\\[YR\\1\\])?$")],
  id_like
)
stopifnot("No year columns detected — check sheet headers." = length(year_cols) > 0)

# -------- Robust Series Name → variable mapping (regex) --------
map_key <- tibble(
  pattern = c(
    "^GDP per capita",                           # gdp_pc_cur
    "^GDP \\(current US\\$\\)",                  # gdp_current_usd
    "^Population, total",                        # pop
    "^Incidence of malaria",                     # mal_incidence
    "^Infection Prevalence",                     # mal_prev
    "^Average precipitation in depth",           # precip_mm (time-invariant; may be dropped with FE)
    "^Agriculture, forestry, and fishing",       # agri_gdp
    "^Political Stability and Absence",          # polstab
    "^Fertility rate",                           # fert
    "^Prevalence of HIV, total",                 # hiv_prev
    "^Primary completion rate",                  # primary_completion
    "^School enrollment, primary",               # enroll_primary_gross
    "^Educational attainment, at least completed primary",  # attain_primary_25p
    "^Urban population"                          # urban_pct (not used here but mapped)
  ),
  var = c("gdp_pc_cur","gdp_current_usd","pop",
          "mal_incidence","mal_prev","precip_mm","agri_gdp","polstab",
          "fert","hiv_prev","primary_completion","enroll_primary_gross",
          "attain_primary_25p","urban_pct")
)

match_var <- function(sn) {
  hit <- map_key$var[str_detect(sn, regex(map_key$pattern, ignore_case = TRUE))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# -------- Long → assign vars → wide panel --------
long_df <- raw_wide %>%
  pivot_longer(all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
  mutate(
    year = readr::parse_number(year_raw),
    var  = vapply(`Series Name`, match_var, character(1)),
    value = readr::parse_number(as.character(val_raw))
  ) %>%
  filter(!is.na(var))

panel <- long_df %>%
  transmute(country = `Country Name`, year, var, value) %>%
  group_by(country, year, var) %>%
  summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
  pivot_wider(id_cols = c(country, year), names_from = var, values_from = value) %>%
  arrange(country, year)

# -------- Build outcome, post, exposure --------
panel <- panel %>%
  mutate(
    gdp_pc_cur = coalesce(gdp_pc_cur, gdp_current_usd / pop),
    ln_gdp_pc  = ifelse(is.finite(gdp_pc_cur) & gdp_pc_cur > 0, log(gdp_pc_cur), NA_real_),
    post       = as.integer(year >= 2005)
  )

expo <- panel %>%
  filter(year %in% 2000:2002) %>%
  group_by(country) %>%
  summarise(
    exposure_raw = coalesce(mean(mal_incidence, na.rm = TRUE),
                            mean(mal_prev,      na.rm = TRUE)),
    .groups = "drop"
  )
stopifnot("No malaria exposure available in 2000–2002 (check incidence/prevalence columns)." =
            any(is.finite(expo$exposure_raw)))

expo <- expo %>% mutate(exposure_std = as.numeric(scale(exposure_raw)))
panel <- panel %>% left_join(expo, by = "country") %>% filter(year >= 2000)

# -------- Create lags ONLY if sources exist --------
panel <- panel %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    polstab_l1 = if ("polstab"   %in% names(cur_data_all())) lag(polstab, 1)   else NULL,
    fert_l1    = if ("fert"      %in% names(cur_data_all())) lag(fert, 1)      else NULL,
    hiv_l1     = if ("hiv_prev"  %in% names(cur_data_all())) lag(hiv_prev, 1)  else NULL
  ) %>%
  ungroup()

# -------- Helpers for estimation --------
dropna_model_data <- function(dat, y, rhs_vars) {
  need <- c("country","year", y, rhs_vars)
  dd <- dat %>% select(any_of(need)) %>% tidyr::drop_na()
  if (nrow(dd) == 0) NULL else dd
}

run_fe <- function(y, rhs, fe = "country + year", dat, cluster_country = TRUE) {
  rhs_vars <- all.vars(as.formula(paste("~", rhs)))
  dd <- dropna_model_data(dat, y, rhs_vars)
  if (is.null(dd)) return(NULL)
  f <- as.formula(paste0(y, " ~ ", rhs, " | ", fe))
  if (cluster_country) feols(f, data = dd, cluster = ~ country) else feols(f, data = dd)
}

run_fe_trend <- function(y, rhs, dat) {
  rhs_vars <- all.vars(as.formula(paste("~", rhs)))
  dd <- dat %>%
    select(any_of(c("country","year", y, rhs_vars))) %>%
    tidyr::drop_na() %>%
    group_by(country) %>%
    mutate(trend = year - min(year, na.rm = TRUE)) %>%
    ungroup()
  if (nrow(dd) == 0) return(NULL)
  f <- as.formula(paste0(y, " ~ ", rhs, " | country[trend] + year"))
  feols(f, data = dd, cluster = ~ country)
}

# -------- Baseline & main variants (GDP per capita, log) --------
m_base  <- run_fe("ln_gdp_pc", "exposure_std:post", "country + year", panel)

# SAFE variables: time-varying + plausibly exogenous
safe_vars <- intersect(c("agri_gdp","polstab_l1"), names(panel))  # (precip_mm is time-invariant; FE will drop it)
m_safe <- if (length(safe_vars) > 0) {
  run_fe("ln_gdp_pc", paste0("exposure_std:post + ", paste(safe_vars, collapse=" + ")),
         "country + year", panel)
} else NULL

# Mediators (robustness only)
mediator_vars <- intersect(c("fert_l1","hiv_l1","health_gdp"), names(panel))
m_all <- if (length(mediator_vars) > 0) {
  base_rhs <- "exposure_std:post"
  if (!is.null(m_safe) && length(safe_vars) > 0) base_rhs <- paste(base_rhs, paste(safe_vars, collapse = " + "), sep = " + ")
  run_fe("ln_gdp_pc", paste(base_rhs, paste(mediator_vars, collapse = " + "), sep = " + "),
         "country + year", panel)
} else NULL

m_trend <- run_fe_trend("ln_gdp_pc", "exposure_std:post", panel)

# -------- Event-study pre-trend test --------
es_tex_path <- file.path(outdir, "pretrend_test.tex")
es_dat <- panel %>% select(country, year, ln_gdp_pc, exposure_std) %>% drop_na() %>% mutate(rel = year - 2004)
if (nrow(es_dat) > 0 && any(es_dat$rel < 0)) {
  ref_rel <- max(unique(es_dat$rel[es_dat$rel < 0]), na.rm = TRUE)
  es_m <- feols(ln_gdp_pc ~ i(rel, exposure_std, ref = ref_rel) | country + year, data = es_dat, cluster = ~ country)
  pre_levels <- sort(unique(es_dat$rel[es_dat$rel < 0]))
  wald_out <- tryCatch(fixest::wald(es_m, paste0("rel::", pre_levels, " = 0")), error = function(e) NULL)
  if (!is.null(wald_out)) {
    pval <- tryCatch(as.numeric(wald_out$p.value), error = function(e) NA_real_)
    cat(sprintf("%% Pre-trend joint test (all pre coefficients = 0)\nPre-trend p-value: $%s$.\n",
                ifelse(is.na(pval), "N/A", formatC(pval, format = "f", digits = 3))),
        file = es_tex_path)
    message("Wrote: ", es_tex_path)
  }
} else {
  message("Event-study skipped (no pre-period available).")
}

# -------- Education outcomes (if present) --------
edu_models <- list()
if ("enroll_primary_gross" %in% names(panel))
  edu_models$Enroll <- run_fe("enroll_primary_gross", "exposure_std:post", "country + year", panel)
if ("primary_completion" %in% names(panel))
  edu_models$Completion <- run_fe("primary_completion", "exposure_std:post", "country + year", panel)
if ("attain_primary_25p" %in% names(panel))
  edu_models$Attain25p <- run_fe("attain_primary_25p", "exposure_std:post", "country + year", panel)

# -------- Exposure sensitivity --------
sens_models <- list()

# Incidence-only (2000–2002)
if ("mal_incidence" %in% names(panel)) {
  ex_inc <- panel %>% filter(year %in% 2000:2002) %>%
    group_by(country) %>% summarise(exposure_raw = mean(mal_incidence, na.rm = TRUE), .groups="drop") %>%
    mutate(exposure_std_alt = as.numeric(scale(exposure_raw)))
  p_inc <- panel %>% left_join(ex_inc, by = "country")
  sens_models$IncOnly <- run_fe("ln_gdp_pc", "exposure_std_alt:post", "country + year", p_inc)
}
# Prevalence-only (2000–2002)
if ("mal_prev" %in% names(panel)) {
  ex_prev <- panel %>% filter(year %in% 2000:2002) %>%
    group_by(country) %>% summarise(exposure_raw = mean(mal_prev, na.rm = TRUE), .groups="drop") %>%
    mutate(exposure_std_alt = as.numeric(scale(exposure_raw)))
  p_prev <- panel %>% left_join(ex_prev, by = "country")
  sens_models$PrevOnly <- run_fe("ln_gdp_pc", "exposure_std_alt:post", "country + year", p_prev)
}
# Window shift (2001–2003; incidence then prevalence fallback)
ex_alt <- panel %>% filter(year %in% 2001:2003) %>%
  group_by(country) %>% summarise(exposure_raw = coalesce(mean(mal_incidence, na.rm=TRUE),
                                                          mean(mal_prev, na.rm=TRUE)), .groups="drop") %>%
  mutate(exposure_std_alt = as.numeric(scale(exposure_raw)))
p_alt <- panel %>% left_join(ex_alt, by = "country")
sens_models$Shift0103 <- run_fe("ln_gdp_pc", "exposure_std_alt:post", "country + year", p_alt)

sens_models <- Filter(Negate(is.null), sens_models)

# -------- LaTeX outputs --------
main_list <- Filter(Negate(is.null),
                    list(Baseline = m_base, `+ Safe` = m_safe, `+ All` = m_all, `+ Trends` = m_trend))
if (length(main_list) > 0) {
  texreg(main_list,
         custom.model.names = names(main_list),
         custom.coef.map = list(`exposure_std:post` = "Exposure $\\times$ Post",
                                `exposure_std_alt:post` = "Exposure (alt) $\\times$ Post"),
         stars = c(0.10, 0.05, 0.01),
         caption = "GDP per capita (log): Baseline DiD, + Safe controls, + All (mediators as robustness), and + Country trends.",
         label   = "tab:ssa_main_compare",
         booktabs = TRUE, use.packages = FALSE,
         file = file.path(outdir, "table_main_compare.tex"))
  message("Wrote: ", file.path(outdir, "table_main_compare.tex"))
}

if (length(edu_models) > 0) {
  texreg(edu_models,
         custom.model.names = names(edu_models),
         custom.coef.map = list(`exposure_std:post` = "Exposure $\\times$ Post"),
         stars = c(0.10, 0.05, 0.01),
         caption = "Education outcomes: DiD with country and year FE.",
         label   = "tab:ssa_education",
         booktabs = TRUE, use.packages = FALSE,
         file = file.path(outdir, "table_education.tex"))
  message("Wrote: ", file.path(outdir, "table_education.tex"))
}

if (length(sens_models) > 0) {
  texreg(sens_models,
         custom.model.names = names(sens_models),
         custom.coef.map = list(`exposure_std_alt:post` = "Exposure (alt) $\\times$ Post"),
         stars = c(0.10, 0.05, 0.01),
         caption = "Sensitivity to exposure definition and baseline window.",
         label   = "tab:ssa_exposure_sensitivity",
         booktabs = TRUE, use.packages = FALSE,
         file = file.path(outdir, "table_exposure_sensitivity.tex"))
  message("Wrote: ", file.path(outdir, "table_exposure_sensitivity.tex"))
}

# -------- Console hints --------
cat("\nSafe variables present:", paste(intersect(c("agri_gdp","polstab_l1"), names(panel)), collapse=", "), "\n")
cat("Potential mediators present:", paste(intersect(c("fert_l1","hiv_l1","health_gdp"), names(panel)), collapse=", "), "\n")
