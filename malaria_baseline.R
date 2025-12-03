############################################################
# Baseline DiD: SSA only — Exposure × Post + FE (no controls)
############################################################

# ---- Packages ----
library(tidyverse)
library(readxl)
library(janitor)
library(fixest)
library(readr)
library(textreg)

# ---- Paths ----
outdir     <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path  <- file.path(outdir, "Working Data_v3 (1).xlsx")

# ---- Import (treat .. etc. as NA) ----
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")

raw_wide <- read_xlsx(xlsx_path, na = na_tokens, guess_max = 1e5) %>%
  remove_empty("cols") %>%
  select(-matches("^Unnamed")) %>%
  mutate(across(where(is.character), ~ na_if(.x, ".."))) %>%
  mutate(across(where(is.character), ~ na_if(.x, "…")))

# Detect year columns (handles '2000' and '2000 [YR2000]')
id_like   <- c("Series Name","Series Code","Country Name","Country Code")
all_cols  <- colnames(raw_wide)
year_cols <- setdiff(all_cols[str_detect(all_cols, "^(\\d{4})(\\s*\\[YR\\1\\])?$")], id_like)
stopifnot(length(year_cols) > 0)

# ---- Long → map variables → panel (country-year) ----
var_map <- c(
  "GDP per capita (current US$)" = "gdp_pc_cur",
  "GDP (current US$)" = "gdp_current_usd",
  "Population, total" = "pop",
  "Incidence of malaria (per 1,000 population at risk)" = "mal_incidence",
  "Infection Prevalence (per 100 Children)" = "mal_prev"
)

panel <- raw_wide %>%
  pivot_longer(cols = all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
  mutate(year = readr::parse_number(year_raw),
         series_name = `Series Name`,
         var = unname(var_map[series_name])) %>%
  filter(!is.na(var)) %>%
  mutate(value = readr::parse_number(as.character(val_raw))) %>%
  transmute(country = `Country Name`, code = `Country Code`, year, var, value) %>%
  group_by(country, year, var) %>% summarise(value = dplyr::first(value[!is.na(value)]), .groups="drop") %>%
  pivot_wider(id_cols = c(country, year), names_from = var, values_from = value) %>%
  arrange(country, year)

# Build GDPpc if needed; log outcome; Post; Exposure (avg 2000–2002, std)
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
  ) %>%
  mutate(exposure_std = as.numeric(scale(exposure_raw)))

panel <- panel %>%
  left_join(expo, by = "country") %>%
  filter(year >= 2000)

# ---- Helper to run baseline on any outcome ----
run_baseline <- function(dat, yvar) {
  dd <- dat %>% select(country, year, exposure_std, post, all_of(yvar)) %>% drop_na()
  if (nrow(dd) == 0) return(NULL)
  f <- as.formula(paste0(yvar, " ~ exposure_std:post | country + year"))
  feols(f, data = dd, cluster = ~ country)
}

# ---- Outcomes you care about (add more if needed) ----
mods <- list(
  gdp  = run_baseline(panel, "ln_gdp_pc")
  # If you want education outcomes, map them first in var_map and then add here, e.g.:
  # enroll = run_baseline(panel, "enroll_primary_gross"),
  # comp   = run_baseline(panel, "primary_completion"),
  # attain = run_baseline(panel, "attain_primary_25p")
)

# ---- Print results ----
lapply(mods, function(m) if (!is.null(m)) { print(m); cat("\n") })

# ---- LaTeX export (texreg) ----
# install.packages("texreg")  # run once if needed
library(texreg)

# keep only non-NULL models
valid_mods <- Filter(Negate(is.null), mods)

if (length(valid_mods) > 0) {
  texreg(
    valid_mods,
    custom.model.names = names(valid_mods),                     # e.g., "gdp"
    stars   = c(0.10, 0.05, 0.01),                              # *, **, ***
    file    = file.path(outdir, "table_baseline_ssa.tex"),      # <-- output path
    caption = "Baseline SSA DiD: Exposure\\,\\(\\times\\)\\,Post with country and year FE.",
    label   = "tab:ssa_baseline",
    booktabs = TRUE, use.packages = FALSE
  )
  message("LaTeX table written to: ", file.path(outdir, "table_baseline_ssa.tex"),
          "\nUse '\\input{table_baseline_ssa.tex}' in your .tex file.")
} else {
  message("No valid models to export.")
}
