############################################################
# SSA DiD — Preferred spec across multiple outcomes
# RHS: exposure_std:post + agri_gdp + polstab_l1
# + Appendix: country-trend robustness
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(fixest)
  library(readr)
})

# -------- Paths --------
outdir     <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path  <- file.path(outdir, "Working Data_v3 (1).xlsx")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------- Import & tidy --------
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")

raw_wide <- read_xlsx(xlsx_path, na = na_tokens, guess_max = 1e5) %>%
  remove_empty("cols") %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-matches("^Unnamed")) %>%
  mutate(across(where(is.character), ~ na_if(.x, "..")),
         across(where(is.character), ~ na_if(.x, "…")))

# Detect year columns
id_like   <- c("Series Name","Series Code","Country Name","Country Code")
year_cols <- setdiff(names(raw_wide)[str_detect(names(raw_wide), "^(\\d{4})(\\s*\\[YR\\1\\])?$")], id_like)
stopifnot("No year columns detected — check sheet headers." = length(year_cols) > 0)

# -------- Series-name → variable mapping (regex) --------
map_key <- tibble(
  pattern = c(
    "^GDP per capita", "^GDP \\(current US\\$\\)", "^Population, total",
    "^Incidence of malaria", "^Infection Prevalence",
    "^Agriculture, forestry, and fishing", "^Political Stability and Absence",
    "^Primary completion rate", "^School enrollment, primary",
    "^Educational attainment, at least completed primary", "^Urban population"
  ),
  var = c(
    "gdp_pc_cur","gdp_current_usd","pop",
    "mal_incidence","mal_prev",
    "agri_gdp","polstab",
    "primary_completion","enroll_primary_gross",
    "attain_primary_25p","urban_pct"
  )
)
match_var <- function(sn) {
  hit <- map_key$var[str_detect(sn, regex(map_key$pattern, ignore_case = TRUE))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# -------- Long → wide panel --------
panel <- raw_wide %>%
  pivot_longer(all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
  mutate(year = readr::parse_number(year_raw),
         var  = vapply(`Series Name`, match_var, character(1)),
         value = readr::parse_number(as.character(val_raw))) %>%
  filter(!is.na(var)) %>%
  transmute(country = `Country Name`, year, var, value) %>%
  group_by(country, year, var) %>%
  summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
  pivot_wider(id_cols = c(country, year), names_from = var, values_from = value) %>%
  arrange(country, year)

# -------- Outcome, Post, Exposure (std, 2000–2002) --------
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
stopifnot("Exposure could not be computed for 2000–2002." = any(is.finite(expo$exposure_raw)))
panel <- panel %>%
  left_join(expo %>% mutate(exposure_std = as.numeric(scale(exposure_raw))), by = "country") %>%
  filter(year >= 2000)

# Lag political stability if available
panel <- panel %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(polstab_l1 = if ("polstab" %in% names(cur_data_all())) lag(polstab, 1) else polstab_l1) %>%
  ungroup()

# -------- Preferred RHS --------
rhs <- "exposure_std:post + agri_gdp + polstab_l1"

# Outcomes as columns (only those present will be estimated)
outs <- c("ln_gdp_pc", "enroll_primary_gross", "primary_completion", "attain_primary_25p")
out_labels <- c("ln GDPpc", "Enroll (Primary)", "Primary Completion", "Attain (Primary, 25+)")

run_one <- function(y, fe_spec = "country + year") {
  vars <- all.vars(as.formula(paste("~", rhs)))
  dd <- panel %>% select(any_of(c("country","year", y, vars))) %>% tidyr::drop_na()
  if (nrow(dd) == 0) return(NULL)
  feols(as.formula(paste0(y, " ~ ", rhs, " | ", fe_spec)), data = dd, cluster = ~ country)
}

# ---------- Preferred (country + year FE) ----------
mods <- lapply(outs, run_one, fe_spec = "country + year")
names(mods) <- out_labels
mods <- Filter(Negate(is.null), mods)
stopifnot(length(mods) > 0)

dict <- c(
  "exposure_std:post" = "Exposure × Post",
  "agri_gdp"          = "Agri share of GDP",
  "polstab_l1"        = "Political Stability (lag)"
)

etable(
  mods,
  se        = "cluster",
  dict      = dict,
  fitstat   = c("n","g","r2","pr2","ar2","apr2"),
  drop      = "(Intercept)",
  style.tex = style.tex("aer"),
  tex       = TRUE,
  title     = "SSA: Preferred specification across outcomes (country and year FE; clustered by country).",
  label     = "tab:ssa_multi_outcomes_preferred",
  file      = file.path(outdir, "table_ssa_multi_outcomes_preferred.tex")
)

# Also print ln GDPpc column (if present) for a quick console check
if ("ln GDPpc" %in% names(mods)) print(mods[["ln GDPpc"]]) else print(mods[[1]])

# ---------- Appendix: Country-trend robustness ----------
# Ensure year is numeric (for trend construction)
panel <- panel %>% mutate(year = as.numeric(year))

mods_trend <- lapply(outs, function(y) {
  vars <- all.vars(as.formula(paste("~", rhs)))
  dd <- panel %>% select(any_of(c("country","year", y, vars))) %>% drop_na()
  if (nrow(dd) == 0) return(NULL)
  feols(as.formula(paste0(y, " ~ ", rhs, " | country[year] + year")),
        data = dd, cluster = ~ country)
}) %>% Filter(Negate(is.null), .)

if (length(mods_trend) > 0) {
  etable(
    mods_trend,
    se        = "cluster",
    dict      = dict,
    fitstat   = c("n","g","r2","pr2","ar2","apr2"),
    drop      = "(Intercept)",
    style.tex = style.tex("aer"),
    tex       = TRUE,
    title     = "SSA: Preferred specification with country-specific linear trends.",
    label     = "tab:ssa_multi_outcomes_trend",
    file      = file.path(outdir, "table_ssa_multi_outcomes_trend.tex")
  )
  message("Saved appendix trend table: table_ssa_multi_outcomes_trend.tex")
} else {
  message("Trend models skipped (insufficient data after drops).")
}

