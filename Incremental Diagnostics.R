############################################################
# Malaria → SSA: Baseline DiD + Incremental Control Diagnostics
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(fixest)
  library(texreg)
  library(readr)
})

# ========= Paths (EDIT if needed) =========
outdir     <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path  <- file.path(outdir, "Working Data_v3 (1).xlsx")
sheet_name <- NULL  # or "SSA - Combined by Countries"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ========= Import & tidy =========
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

id_like   <- c("Series Name","Series Code","Country Name","Country Code")
year_cols <- setdiff(
  names(raw_wide)[str_detect(names(raw_wide), "^(\\d{4})(\\s*\\[YR\\1\\])?$")],
  id_like
)
stopifnot("No year columns detected — check sheet headers." = length(year_cols) > 0)

# ========= Robust Series Name → variable mapping (regex) =========
map_key <- tibble(
  pattern = c(
    "^GDP per capita",
    "^GDP \\(current US\\$\\)",
    "^Population, total",
    "^Incidence of malaria",
    "^Infection Prevalence",
    "^Average precipitation in depth",
    "^Agriculture, forestry, and fishing",
    "^Political Stability and Absence",
    "^Fertility rate",
    "^Prevalence of HIV, total",
    "^Primary completion rate",
    "^School enrollment, primary",
    "^Educational attainment, at least completed primary",
    "^Urban population"
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

# ========= Long → assign vars → wide panel =========
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

# ========= Construct outcome, post, exposure =========
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
stopifnot("No malaria exposure in 2000–2002 (check incidence/prevalence columns)." =
            any(is.finite(expo$exposure_raw)))

expo <- expo %>% mutate(exposure_std = as.numeric(scale(exposure_raw)))
panel <- panel %>% left_join(expo, by = "country") %>% filter(year >= 2000)

# ========= Create lags (only if sources exist) =========
panel <- panel %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    polstab_l1 = if ("polstab"   %in% names(cur_data_all())) lag(polstab, 1)   else NULL,
    fert_l1    = if ("fert"      %in% names(cur_data_all())) lag(fert, 1)      else NULL,
    hiv_l1     = if ("hiv_prev"  %in% names(cur_data_all())) lag(hiv_prev, 1)  else NULL,
    u5_mort_l5 = if ("u5_mort"   %in% names(cur_data_all())) lag(u5_mort, 5)   else NULL
  ) %>% ungroup()

# ========= Helpers =========
dropna_model_data <- function(dat, y, rhs_vars) {
  need <- c("country","year", y, rhs_vars)
  dd <- dat %>% select(any_of(need)) %>% tidyr::drop_na()
  if (nrow(dd) == 0) NULL else dd
}

run_fe <- function(y, rhs, dat, fe = "country + year") {
  rhs_vars <- all.vars(as.formula(paste("~", rhs)))
  dd <- dropna_model_data(dat, y, rhs_vars)
  if (is.null(dd)) return(NULL)
  feols(as.formula(paste0(y," ~ ", rhs, " | ", fe)), data = dd, cluster = ~ country)
}

# ========= Baseline DiD =========
m_base <- run_fe("ln_gdp_pc", "exposure_std:post", panel)
stopifnot(!is.null(m_base))
b0  <- coef(m_base)["exposure_std:post"]; se0 <- se(m_base)["exposure_std:post"]

# ========= Candidate control sets =========
# SAFE (time-varying, plausibly exogenous)
safe_candidates       <- intersect(c("agri_gdp","polstab_l1","urban_pct"), names(panel))
# Time-invariant (dropped by country FE; kept only for audit)
time_invariant_struct <- intersect(c("precip_mm"), names(panel))
# Potential MEDIATORS (use for robustness only)
mediator_candidates   <- intersect(c("fert_l1","hiv_l1","health_gdp","u5_mort_l5"), names(panel))

message("Safe candidates: ", paste(safe_candidates, collapse = ", "))
message("Mediator candidates: ", paste(mediator_candidates, collapse = ", "))

# ========= Diagnostic: add each control individually =========
test_one <- function(var) {
  rhs <- paste0("exposure_std:post + ", var)
  m   <- run_fe("ln_gdp_pc", rhs, panel)
  if (is.null(m)) return(NULL)
  b   <- coef(m)["exposure_std:post"]; s <- se(m)["exposure_std:post"]
  tibble(
    var = var,
    beta = as.numeric(b), se = as.numeric(s),
    nobs = nobs(m),
    shift_abs = as.numeric(b - b0),
    shift_pct = ifelse(isTRUE(abs(b0) > 0),
                       as.numeric((b - b0)/abs(b0))*100,
                       NA_real_),
    sign_flip = sign(b0) != sign(b),
    se_ratio  = as.numeric(s / se0),
    type = case_when(
      var %in% safe_candidates ~ "safe",
      var %in% mediator_candidates ~ "mediator",
      var %in% time_invariant_struct ~ "time_invariant",
      TRUE ~ "other"
    )
  )
}

diag_list <- lapply(c(safe_candidates, mediator_candidates, time_invariant_struct), test_one)
diag_tbl  <- bind_rows(Filter(Negate(is.null), diag_list)) %>%
  arrange(desc(abs(shift_pct))) %>%
  mutate(flag_big = (abs(shift_pct) >= 25) | sign_flip)

# ========= Diagnostic: add pairs of safe controls =========
pair_tbl <- tibble()
if (length(safe_candidates) >= 2) {
  pairs <- combn(safe_candidates, 2, simplify = FALSE)
  for (p in pairs) {
    rhs <- paste0("exposure_std:post + ", paste(p, collapse = " + "))
    m   <- run_fe("ln_gdp_pc", rhs, panel)
    if (!is.null(m)) {
      b <- coef(m)["exposure_std:post"]; s <- se(m)["exposure_std:post"]
      pair_tbl <- bind_rows(pair_tbl, tibble(
        var = paste(p, collapse = " + "),
        beta = as.numeric(b), se = as.numeric(s), nobs = nobs(m),
        shift_abs = as.numeric(b - b0),
        shift_pct = ifelse(isTRUE(abs(b0) > 0),
                           as.numeric((b - b0)/abs(b0))*100,
                           NA_real_),
        sign_flip = sign(b0) != sign(b),
        se_ratio  = as.numeric(s / se0),
        type = "safe_pair",
        flag_big = (abs(shift_pct) >= 25) | sign_flip
      ))
    }
  }
}

diag_tbl_all <- bind_rows(diag_tbl, pair_tbl) %>%
  arrange(type, desc(abs(shift_pct)))

# ========= Save diagnostics (CSV + LaTeX) =========
csv_path <- file.path(outdir, "control_shift_diagnostics.csv")
readr::write_csv(diag_tbl_all, csv_path)

top_tbl <- diag_tbl_all %>%
  mutate(shift_pct = ifelse(is.na(shift_pct), NA, round(shift_pct, 1)),
         beta = round(beta, 3), se = round(se, 3),
         se_ratio = round(se_ratio, 2)) %>%
  select(type, var, beta, se, shift_pct, sign_flip, se_ratio, nobs) %>%
  slice_head(n = 12)

tex_path <- file.path(outdir, "table_control_shifts.tex")
cat("\\begin{table}[!ht]\\centering\n\\caption{Incremental control diagnostics: shift in $\\hat\\beta$ (Exposure $\\times$ Post)}\n",
    "\\begin{tabular}{llrrrrrr}\n\\toprule\n",
    "Type & Control & $\\hat\\beta$ & SE & $\\Delta\\%$ & Flip & SE/SE$_0$ & $N$\\\\\\midrule\n",
    file = tex_path)
apply(top_tbl, 1, function(r) {
  cat(sprintf("%s & %s & %s & %s & %s & %s & %s & %s\\\\\n",
              r[["type"]], r[["var"]], r[["beta"]], r[["se"]],
              ifelse(is.na(as.numeric(r[["shift_pct"]])), "NA", r[["shift_pct"]]),
              ifelse(as.logical(r[["sign_flip"]]), "Yes", "No"),
              r[["se_ratio"]], r[["nobs"]]),
      file = tex_path, append = TRUE)
})
cat("\\bottomrule\n\\end{tabular}\n\\label{tab:control_shifts}\n\\end{table}\n", file = tex_path, append = TRUE)

message("Wrote: ", csv_path, " and ", tex_path)

# ========= Heuristic recommendations =========
rec_keep   <- diag_tbl %>% filter(type == "safe", flag_big) %>% pull(var)
rec_avoid  <- diag_tbl %>% filter(type == "mediator", flag_big) %>% pull(var)
message("Recommend INCLUDE (safe & big shift): ", paste(rec_keep, collapse = ", "))
message("Recommend AVOID (mediator & big shift): ", paste(rec_avoid, collapse = ", "))

# ========= Optional: preferred model with recommended safe controls =========
if (length(rec_keep) > 0) {
  rhs_pref <- paste0("exposure_std:post + ", paste(unique(rec_keep), collapse = " + "))
  m_pref   <- run_fe("ln_gdp_pc", rhs_pref, panel)
  if (!is.null(m_pref)) {
    print(m_pref)
    texreg(list(`Preferred (+selected safe)` = m_pref),
           custom.coef.map = list(`exposure_std:post` = "Exposure $\\times$ Post"),
           stars = c(0.10, 0.05, 0.01),
           caption = "Preferred GDPpc model with selected safe controls.",
           label   = "tab:ssa_preferred",
           booktabs = TRUE, use.packages = FALSE,
           file = file.path(outdir, "table_preferred.tex"))
  }
}

# ========= Final console summary =========
cat("\nBaseline beta:", round(as.numeric(b0), 4), "  SE:", round(as.numeric(se0), 4), "\n")
cat("Top movers saved in:", tex_path, "\n")
