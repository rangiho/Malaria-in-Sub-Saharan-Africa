############################################################
# Malaria → Education & Economic Outcomes in SSA (R)
# Panel econometrics template (country × year) — SSA ONLY
# Author: <your name>
# Date: <today>
############################################################

# 0) Packages ---------------------------------------------------------------
req <- c(
  "tidyverse","data.table","readxl","janitor","countrycode",
  "fixest","modelsummary","broom","broom.helpers","rlang"
)
invisible(lapply(req, function(p) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)))
invisible(lapply(req, library, character.only = TRUE))

# 0a) Helpers ---------------------------------------------------------------
# Safely pick the first existing column from a list of candidates
col_first <- function(df, ...) {
  nms <- c(...)
  for (nm in nms) if (nm %in% names(df)) return(df[[nm]])
  return(rep(NA, nrow(df)))
}

# 1) Data import ------------------------------------------------------------
# IMPORTANT: set your file path to the XLSX (v3)
path_raw   <- "/Users/rangiho/Desktop/2025/Y4S1/Working Data_v3 (1).xlsx"  # <-- adjust if needed
sheet_name <- Countries (SSA+India)                                                            # or actual sheet name

df0 <- readxl::read_xlsx(path_raw, sheet = sheet_name) |>
  janitor::clean_names()

# 1a) Column mapping (robust) ----------------------------------------------
# Expected concepts (names in your file may differ; we choose the first that exists)
df <- tibble::tibble(
  # IDs
  iso3    = toupper(if ("iso3" %in% names(df0)) df0$iso3 else countrycode(df0$country, "country.name", "iso3c")),
  country = col_first(df0, "country", "country_name"),
  year    = as.integer(col_first(df0, "year")),
  
  # Outcomes
  gdp_pc_cur = suppressWarnings(as.numeric(col_first(df0, "gdp_per_capita_current_usd"))),
  enroll_primary_gross = suppressWarnings(as.numeric(col_first(
    df0,
    "school_enrollment_primary_percent_gross",
    "primary_enroll_gross_pct"
  ))),
  primary_completion = suppressWarnings(as.numeric(col_first(
    df0, "primary_completion_pct"
  ))),
  attain_primary_25p = suppressWarnings(as.numeric(col_first(
    df0,
    "educational_attainment_at_least_completed_primary_population_25_years_total_percent_cumulative",
    "attain_primary_25p"
  ))),
  
  # Malaria burden / intervention
  mal_incidence = suppressWarnings(as.numeric(col_first(
    df0,
    "malaria_incidence_per_1000_population_at_risk",
    "malaria_incidence_per_1000"
  ))),
  mal_prev = suppressWarnings(as.numeric(col_first(
    df0, "infection_prevalence_per_100_children", "malaria_prevalence_pfpr"
  ))),
  mal_mort100k = suppressWarnings(as.numeric(col_first(
    df0, "mortality_rate_per_100_000", "malaria_mort_per_100k"
  ))),
  itn_u5 = suppressWarnings(as.numeric(col_first(
    df0,
    "use_of_insecticide_treated_bed_nets_percent_of_under_5_population",
    "itn_use_u5"
  ))),
  antimal_u5 = suppressWarnings(as.numeric(col_first(
    df0,
    "children_with_fever_receiving_antimalarial_drugs_percent_of_children_under_age_5_with_fever",
    "fever_antimalarial_u5"
  ))),
  
  # Controls
  u5_mort = suppressWarnings(as.numeric(col_first(
    df0, "mortality_rate_under_5_per_1_000_live_births", "u5_mort_per_1000"
  ))),
  fert = suppressWarnings(as.numeric(col_first(
    df0, "fertility_rate_total_births_per_woman", "fert_rate_total"
  ))),
  hiv_prev = suppressWarnings(as.numeric(col_first(
    df0, "prevalence_of_hiv_total_percent_of_population_ages_15_49", "hiv_prev_1549"
  ))),
  health_gdp = suppressWarnings(as.numeric(col_first(
    df0, "domestic_general_government_health_expenditure_percent_of_gdp", "gov_health_exp_gdp"
  ))),
  agri_gdp = suppressWarnings(as.numeric(col_first(
    df0, "agriculture_forestry_and_fishing_value_added_percent_of_gdp", "agri_share_gdp"
  ))),
  precip_mm = suppressWarnings(as.numeric(col_first(
    df0, "average_precipitation_in_depth_mm_per_year", "precip_mm_year"
  ))),
  polstab = suppressWarnings(as.numeric(col_first(
    df0, "political_stability_and_absence_of_violence_terrorism_percentile_rank", "polstab_percentile"
  ))),
  pop = suppressWarnings(as.numeric(col_first(df0, "population_total"))),
  urban_pct = suppressWarnings(as.numeric(col_first(
    df0, "urban_population_percent_of_total_population", "urban_pct_total"
  ))),
  
  # Optional: GDP level current USD (if needed)
  gdp_cur = suppressWarnings(as.numeric(col_first(df0, "gdp_current_usd")))
) |>
  arrange(iso3, year)

# 1b) SSA-only filter -------------------------------------------------------
ssa <- countrycode::countrycode_data |>
  dplyr::filter(region == "Sub-Saharan Africa") |>
  dplyr::distinct(iso3c) |>
  dplyr::pull(iso3c)

df <- df |>
  dplyr::filter(iso3 %in% ssa)

# 2) Basic cleaning & transformations ---------------------------------------
df <- df |>
  mutate(
    ln_gdp_pc = log(pmax(gdp_pc_cur, 1)),           # guard log
    polstab_z = as.numeric(scale(polstab)),
    across(
      c(enroll_primary_gross, primary_completion, attain_primary_25p,
        itn_u5, antimal_u5, urban_pct, hiv_prev),
      ~ ifelse(is.na(.x), NA_real_, pmin(pmax(as.numeric(.x), 0), 100))
    )
  )

# 3) Baseline exposure (2000–2002) & POST (≥2005) ---------------------------
expo <- df |>
  filter(year %in% 2000:2002) |>
  group_by(iso3) |>
  summarise(
    exposure_inc  = mean(mal_incidence, na.rm = TRUE),
    exposure_prev = mean(mal_prev,      na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    exposure_inc_std  = as.numeric(scale(exposure_inc)),
    exposure_prev_std = as.numeric(scale(exposure_prev))
  )

df <- df |>
  left_join(expo, by = "iso3") |>
  mutate(post = as.integer(year >= 2005))

# 4) Lags for core controls --------------------------------------------------
setDT(df)
df <- df[order(iso3, year)]
df[, ln_gdp_pc_l1 := shift(ln_gdp_pc, 1), by = iso3]
df[, fert_l1      := shift(fert, 1),_]()
