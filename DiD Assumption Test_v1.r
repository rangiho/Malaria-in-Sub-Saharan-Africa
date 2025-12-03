# ========= Setup =========
library(fixest); library(dplyr); library(broom); library(texreg)
outdir <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"  # change if needed

controls_rhs <- ~ u5_mort_l5 + fert_l1 + hiv_l1 + urban_pct + precip_mm +
  health_gdp + agri_gdp + polstab_l1

# Make relative year (ref = 2004) for SSA panel
panel_ssa_es <- panel_ssa %>% mutate(rel = year - 2004)

# --- GDP per capita event-study (SSA) ---
es_gdp <- feols(
  ln_gdp_pc ~ i(rel, exposure_std, ref = -1) + ..controls_rhs | country + year,
  data = panel_ssa_es, cluster = ~ country
)

# Joint pre-trend test: all leads k <= -2
pre_test_gdp <- wald(es_gdp, keep = "^rel::(-[2-9]|-1[0-9]+)", cluster = "country")

# Export regression (full dynamic) as a compact LaTeX table
texreg(list(es_gdp),
       custom.model.names = "ES: ln GDPpc (SSA)",
       file = file.path(outdir, "es_lnGDPpc_ssa.tex"),
       caption = "Event study for ln(GDP per capita) in SSA (ref = 2004).",
       label = "tab:es_gdp_ssa")

# Save the joint pre-trend test as a tiny LaTeX table
data.frame(
  Test = "Joint pre-trends (k<=-2)",
  F_stat = pre_test_gdp$stat[1],
  df1 = pre_test_gdp$df[1,1], df2 = pre_test_gdp$df[1,2],
  p_value = pre_test_gdp$p[1]
) |>
  xtable::xtable(caption = "Pre-trend joint test for ln(GDP per capita), SSA.",
                 label = "tab:pretrend_gdp_ssa") |>
  print(file = file.path(outdir, "pretrend_gdp_ssa.tex"),
        include.rownames = FALSE)
