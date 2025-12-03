############################################################
# Summary stats for FINAL empirical spec
# Variables: ln_gdp_pc, agri_gdp, polstab_l1, exposure_std, post
# Output: table_summary_final.tex (booktabs LaTeX)
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(readr)
})

# ---------- Paths ----------
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

# ---------- Map series → variables (minimal set for final spec) ----------
map_key <- tibble(
  pattern = c(
    "^GDP per capita", "^GDP \\(current US\\$\\)", "^Population, total",
    "^Incidence of malaria", "^Infection Prevalence",
    "^Agriculture, forestry, and fishing", "^Political Stability and Absence"
  ),
  var = c(
    "gdp_pc_cur","gdp_current_usd","pop",
    "mal_incidence","mal_prev",
    "agri_gdp","polstab"
  )
)
match_var <- function(sn) {
  hit <- map_key$var[str_detect(sn, regex(map_key$pattern, ignore_case = TRUE))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# ---------- Long → wide panel ----------
panel <- raw_wide %>%
  pivot_longer(all_of(year_cols), names_to = "year_raw", values_to = "val_raw") %>%
  mutate(year  = readr::parse_number(year_raw),
         var   = vapply(`Series Name`, match_var, character(1)),
         value = readr::parse_number(as.character(val_raw))) %>%
  filter(!is.na(var)) %>%
  transmute(country = `Country Name`, year, var, value) %>%
  group_by(country, year, var) %>% summarise(value = dplyr::first(value[!is.na(value)]), .groups = "drop") %>%
  pivot_wider(id_cols = c(country, year), names_from = var, values_from = value) %>%
  arrange(country, year)

# ---------- Build final-spec variables ----------
panel <- panel %>%
  mutate(
    gdp_pc_cur = coalesce(gdp_pc_cur, gdp_current_usd / pop),
    ln_gdp_pc  = ifelse(is.finite(gdp_pc_cur) & gdp_pc_cur > 0, log(gdp_pc_cur), NA_real_),
    post       = as.integer(year >= 2005)
  )

# Exposure: avg 2000–2002 (incidence; fallback PfPR), standardized
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
  filter(year >= 2000) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(polstab_l1 = if ("polstab" %in% names(cur_data_all())) lag(polstab, 1) else polstab_l1) %>%
  ungroup()

# ---------- Preferred-spec sample (complete cases) ----------
vars_final <- c("ln_gdp_pc","agri_gdp","polstab_l1","exposure_std","post")
dd <- panel %>% select(country, year, all_of(vars_final)) %>% drop_na()
N  <- nrow(dd)

stopifnot("No rows left after filtering to final-spec variables." = N > 0)

# ---------- Summary stats ----------
summarise_var <- function(x) {
  c(mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE))
}

stats <- map_dfr(vars_final, ~{
  v <- .x
  tibble(
    Variable = c(
      "ln(GDP per capita)",
      "Agriculture share of GDP (%)",
      "Political Stability (lag, 0–100)",
      "Malaria incidence (per 1,000 at risk, baseline SD units)",
      "Post indicator (t ≥ 2005)"
    )[match(v, vars_final)],
    Mean = summarise_var(dd[[v]])["mean"],
    SD   = summarise_var(dd[[v]])["sd"],
    Min  = summarise_var(dd[[v]])["min"],
    Max  = summarise_var(dd[[v]])["max"]
  )
})

# Nicely rounded
stats <- stats %>%
  mutate(across(c(Mean, SD, Min, Max), ~ ifelse(is.finite(.x), round(.x, 2), NA)))

# ---------- Write LaTeX (booktabs) ----------
tex_path <- file.path(outdir, "table_summary_final.tex")
con <- file(tex_path, open = "wt", encoding = "UTF-8")

writeLines(c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Summary Statistics: Variables Used in Final Empirical Framework (SSA, 2000--2023)}",
  "\\label{tab:summary_final}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Variable & Mean & SD & Min & Max \\\\",
  "\\midrule"
), con)

apply(stats, 1, function(r) {
  line <- sprintf("%s & %s & %s & %s & %s \\\\",
                  r[["Variable"]], r[["Mean"]], r[["SD"]], r[["Min"]], r[["Max"]])
  writeLines(line, con)
})

writeLines(c(
  sprintf("Observations & \\multicolumn{4}{c}{%d} \\\\", N),
  "\\bottomrule",
  "\\multicolumn{5}{l}{\\footnotesize{Notes: Statistics computed on the preferred-spec sample (complete cases).}} \\\\",
  "\\end{tabular}",
  "\\end{table}"
), con)

close(con)
cat("Saved LaTeX summary table to:", tex_path, "\n")
