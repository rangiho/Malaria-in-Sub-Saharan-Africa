############################################################
# Malaria controls → GDP & Education (z)
# - Auto-sheet selection
# - Nearest-year match (±3y) for sparse coverage
# - HC3-robust effects for Mosquito nets vs Antimalarial treatment
# - Heterogeneity by baseline income (continuous & groups)
# Outputs (in outdir):
#   - treat_effects_nearest.csv / .tex
#   - heterogeneity_continuous.csv
#   - heterogeneity_groups.csv / .tex
############################################################

suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(janitor); library(readr); library(stringr)
  library(sandwich); library(lmtest); library(broom)
})

# ---------- Paths ----------
outdir    <- "/Users/rangiho/Desktop/2025/Y4S1/malaria"
xlsx_path <- file.path(outdir, "Working Data_v3 (1).xlsx")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Parameters ----------
max_gap <- 3   # nearest-year tolerance (years)

# ---------- Helpers ----------
na_tokens <- c("", "NA", "NaN", "..", "…", "—", "-")
norm <- function(x) x %>% str_replace_all("\\s+"," ") %>% str_trim() %>% tolower()
zfun <- function(v) if (all(is.na(v))) NA_real_ else as.numeric(scale(v))
pctl <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE))

# ---------- Series labels (tolerant matching) ----------
label_nets_exact  <- "Use of insecticide-treated bed nets (% of under-5 population)"
label_treat_exact <- "Children with fever receiving antimalarial drugs (% of children under age 5 with fever)"

# ---------- Sheet scanner → returns panel for best sheet ----------
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
  
  id_like   <- c("Series Name","Series Code","Country Name","Country Code")
  year_cols <- setdiff(names(df)[str_detect(names(df), "^(\\d{4})(\\s*\\[YR\\1\\])?$")], id_like)
  if (!length(year_cols)) return(NULL)
  
  sN <- norm(df$`Series Name`)
  
  is_nets  <- sN == norm(label_nets_exact) |
    str_detect(sN, "insecticide.*bed net|treated.*net.*under.?5")
  is_treat <- sN == norm(label_treat_exact) |
    str_detect(sN, "fever.*antimalarial|receiv.*antimalarial.*fever")
  map <- rep(NA_character_, length(sN))
  map[sN == norm("GDP per capita (current US$)")] <- "gdp_pc_cur"
  map[sN == norm("GDP (current US$)")]            <- "gdp_current_usd"
  map[sN == norm("Population, total")]            <- "pop"
  map[sN == norm("Agriculture, forestry, and fishing, value added (% of GDP)")] <- "agri_gdp"
  map[sN == norm("Political Stability and Absence of Violence/Terrorism: Percentile Rank")] <- "polstab"
  map[sN == norm("School enrollment, primary (% gross)")] <- "enroll_primary_gross"
  map[sN == norm("Primary completion rate, total (% of relevant age group)")]   <- "primary_completion"
  map[sN == norm("Educational attainment, at least completed primary, population 25+ years, total (%) (cumulative)")] <- "attain_primary_25p"
  map[is_nets]  <- "nets_u5"
  map[is_treat] <- "treat_u5"
  df$var_compact <- map
  
  panel <- df %>%
    filter(!is.na(var_compact)) %>%
    pivot_longer(all_of(year_cols), names_to="year_raw", values_to="val_raw") %>%
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
    summarise(value = dplyr::first(value[!is.na(value)]), .groups="drop") %>%
    pivot_wider(id_cols = c(country, code, year),
                names_from = var, values_from = value) %>%
    arrange(country, year) %>%
    mutate(
      gdp_pc_cur = coalesce(gdp_pc_cur, gdp_current_usd / pop),
      ln_gdp_pc  = ifelse(is.finite(gdp_pc_cur) & gdp_pc_cur > 0, log(gdp_pc_cur), NA_real_)
    )
  
  # Coverage score for sheet selection
  cov <- tibble(
    nets  = panel %>% summarise(sum(!is.na(nets_u5)))  %>% pull(),
    treat = panel %>% summarise(sum(!is.na(treat_u5))) %>% pull(),
    gdp   = panel %>% summarise(sum(!is.na(ln_gdp_pc)))%>% pull(),
    edu   = panel %>% summarise(sum(!is.na(enroll_primary_gross) |
                                      !is.na(primary_completion) |
                                      !is.na(attain_primary_25p))) %>% pull()
  )
  score <- 2*min(cov$nets > 0, cov$treat > 0) + 2*min(cov$gdp > 0, cov$edu > 0) +
    cov$nets + cov$treat + cov$gdp + cov$edu
  
  list(sheet = sheet_name, panel = panel,
       coverage = cov, score = as.numeric(score))
}

sheets <- excel_sheets(xlsx_path)
scanned <- purrr::map(sheets, map_and_make_panel) %>% purrr::compact()
stopifnot("No usable sheets found." = length(scanned) > 0)

cov_df <- tibble(
  sheet = sapply(scanned, `[[`, "sheet"),
  nets  = sapply(scanned, function(x) x$coverage$nets),
  treat = sapply(scanned, function(x) x$coverage$treat),
  gdp   = sapply(scanned, function(x) x$coverage$gdp),
  edu   = sapply(scanned, function(x) x$coverage$edu),
  score = sapply(scanned, `[[`, "score")
) %>% arrange(desc(score))

cat("\n--- SHEET COVERAGE ---\n"); print(cov_df, n=Inf)
best <- scanned[[which.max(sapply(scanned, `[[`, "score"))]]
panel <- best$panel
cat("\nChosen sheet:", best$sheet, "\n")

# ---------- Nearest-year matching (by country; robust to missing ISO codes) ----------
by_keys <- "country"

nets_long <- panel %>% select(country, year, nets_u5)  %>% drop_na(nets_u5)  %>% rename(year_r = year)
trt_long  <- panel %>% select(country, year, treat_u5) %>% drop_na(treat_u5) %>% rename(year_r = year)

out_gdp <- panel %>% select(country, year, y = ln_gdp_pc)           %>% drop_na(y) %>% rename(year_l = year)
out_enr <- panel %>% select(country, year, y = enroll_primary_gross) %>% drop_na(y) %>% rename(year_l = year)
out_cpl <- panel %>% select(country, year, y = primary_completion)   %>% drop_na(y) %>% rename(year_l = year)
out_att <- panel %>% select(country, year, y = attain_primary_25p)   %>% drop_na(y) %>% rename(year_l = year)

nearest_for <- function(out_df, right_df, right_val_name) {
  out_df %>%
    inner_join(right_df, by = by_keys, relationship = "many-to-many") %>% # many-to-many is expected
    mutate(gap = abs(year_l - year_r)) %>%
    group_by(country, year_l) %>%
    slice_min(gap, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(gap <= max_gap) %>%
    transmute(country, year = year_l, y, !!right_val_name := .data[[right_val_name]])
}

gdp_nets <- nearest_for(out_gdp, nets_long, "nets_u5")
gdp_trt  <- nearest_for(out_gdp, trt_long,  "treat_u5")

edu_nets <- bind_rows(
  nearest_for(out_enr, nets_long, "nets_u5"),
  nearest_for(out_cpl, nets_long, "nets_u5"),
  nearest_for(out_att, nets_long, "nets_u5")
) %>%
  group_by(country, year) %>%
  summarise(y = mean(y, na.rm = TRUE),
            nets_u5 = mean(nets_u5, na.rm = TRUE), .groups = "drop")

edu_trt <- bind_rows(
  nearest_for(out_enr, trt_long, "treat_u5"),
  nearest_for(out_cpl, trt_long, "treat_u5"),
  nearest_for(out_att, trt_long, "treat_u5")
) %>%
  group_by(country, year) %>%
  summarise(y = mean(y, na.rm = TRUE),
            treat_u5 = mean(treat_u5, na.rm = TRUE), .groups = "drop")

gdp_s <- full_join(gdp_nets, gdp_trt, by = c("country","year","y")) %>%
  mutate(
    z_y     = zfun(y),
    z_nets  = if (all(is.na(nets_u5))) NA_real_ else zfun(nets_u5),
    z_treat = if (all(is.na(treat_u5))) NA_real_ else zfun(treat_u5)
  ) %>%
  filter(!is.na(z_y)) %>% filter(!(is.na(z_nets) & is.na(z_treat)))

edu_s <- full_join(edu_nets, edu_trt, by = c("country","year","y")) %>%
  mutate(
    z_y     = zfun(y),
    z_nets  = if (all(is.na(nets_u5))) NA_real_ else zfun(nets_u5),
    z_treat = if (all(is.na(treat_u5))) NA_real_ else zfun(treat_u5)
  ) %>%
  filter(!is.na(z_y)) %>% filter(!(is.na(z_nets) & is.na(z_treat)))

cat("Matched samples (rows = country×year): GDP =", nrow(gdp_s), " | EDU =", nrow(edu_s), "\n")

# ---------- Minimal OLS (HC3), fallback later if needed ----------
safe_fit <- function(df){
  if (nrow(df) < 3) return(NULL)
  keep <- c("z_nets","z_treat")[colSums(!is.na(df[,c("z_nets","z_treat")])) > 0]
  if (!length(keep)) return(NULL)
  fml <- as.formula(paste0("z_y ~ ", paste(keep, collapse = " + ")))
  m   <- try(lm(fml, data = df %>% drop_na(z_y)), silent = TRUE)
  if (inherits(m,"try-error")) return(NULL)
  vc  <- try(sandwich::vcovHC(m, type = "HC3"), silent = TRUE)
  if (inherits(vc,"try-error")) return(NULL)
  list(model = m)
}

res_gdp <- safe_fit(gdp_s)
res_edu <- safe_fit(edu_s)

# ---------- Extract effects (broom tidy, HC3) ----------
collect_rows <- function(res, outcome_label) {
  if (is.null(res)) return(tibble())
  vc <- sandwich::vcovHC(res$model, type = "HC3")
  broom::tidy(res$model, vcov = vc, conf.int = FALSE) %>%
    filter(term %in% c("z_nets", "z_treat")) %>%
    transmute(
      Outcome      = outcome_label,
      Intervention = recode(term, z_nets="Mosquito nets", z_treat="Antimalarial treatment"),
      Effect       = round(estimate, 3),
      SE           = round(std.error, 3),
      t            = round(statistic, 2),
      p            = signif(p.value, 3),
      Countries    = nobs(res$model),
      Method       = sprintf("OLS (HC3), |Δyear|≤%d", max_gap)
    )
}

results <- bind_rows(
  collect_rows(res_gdp, "GDP (z)"),
  collect_rows(res_edu, "Education Index (z)")
)

# If empty (shouldn't be, but guard): fall back to Spearman
spearman_rows <- function(df, lab){
  out <- list()
  if (sum(complete.cases(df[,c("z_nets","z_y")]))>=3) {
    ct <- suppressWarnings(cor.test(df$z_nets, df$z_y, method="spearman"))
    out[[1]] <- tibble(Outcome=lab, Intervention="Mosquito nets",
                       Effect=round(unname(ct$estimate),3), SE=NA_real_, t=NA_real_,
                       p=signif(ct$p.value,3), Countries=sum(complete.cases(df[,c("z_nets","z_y")])),
                       Method=sprintf("Spearman ρ, |Δyear|≤%d", max_gap))
  }
  if (sum(complete.cases(df[,c("z_treat","z_y")]))>=3) {
    ct <- suppressWarnings(cor.test(df$z_treat, df$z_y, method="spearman"))
    out[[2]] <- tibble(Outcome=lab, Intervention="Antimalarial treatment",
                       Effect=round(unname(ct$estimate),3), SE=NA_real_, t=NA_real_,
                       p=signif(ct$p.value,3), Countries=sum(complete.cases(df[,c("z_treat","z_y")])),
                       Method=sprintf("Spearman ρ, |Δyear|≤%d", max_gap))
  }
  bind_rows(out)
}
if (nrow(results) == 0) {
  results <- bind_rows(
    spearman_rows(gdp_s, "GDP (z)"),
    spearman_rows(edu_s, "Education Index (z)")
  )
}

# ---------- Optional: Wald test H0: nets = treatment ----------
wald_equal <- function(res) {
  if (is.null(res)) return(NA_real_)
  b  <- coef(res$model)
  if (!all(c("z_nets","z_treat") %in% names(b))) return(NA_real_)
  V  <- sandwich::vcovHC(res$model, type = "HC3")
  R  <- c(z_nets = 1, z_treat = -1)
  diff <- sum(R[names(b)] * b)
  var  <- as.numeric(t(R[names(b)]) %*% V[names(b), names(b)] %*% R[names(b)])
  if (!is.finite(var) || var <= 0) return(NA_real_)
  2 * pnorm(-abs(diff / sqrt(var)))
}
p_wald_gdp <- wald_equal(res_gdp); p_wald_edu <- wald_equal(res_edu)

# ---------- Save main effects ----------
csv_main <- file.path(outdir, "treat_effects_nearest.csv")
readr::write_csv(results, csv_main)

tex_main <- file.path(outdir, "treat_effects_nearest.tex")
con <- file(tex_main, open="wt", encoding="UTF-8")
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  sprintf("\\caption{Nets vs Antimalarial Treatment: nearest-year match (|$\\Delta$year|≤%d), standardized}", max_gap),
  "\\label{tab:treat_nearest}",
  "\\small",
  "\\begin{tabular}{@{}lccccccc@{}}\\toprule",
  "Outcome & Intervention & Effect & SE & $t$ & $p$ & Countries & Method \\\\ \\midrule"
), con)
if (nrow(results) == 0) {
  writeLines(" \\multicolumn{8}{c}{No estimable effects or correlations.} \\\\", con)
} else {
  apply(results, 1, function(r) {
    writeLines(sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                       r[["Outcome"]], r[["Intervention"]], r[["Effect"]],
                       ifelse(is.na(r[["SE"]]),"",r[["SE"]]),
                       ifelse(is.na(r[["t"]]),"",r[["t"]]),
                       ifelse(is.na(r[["p"]]),"",r[["p"]]),
                       r[["Countries"]], r[["Method"]]), con)
  })
}
wald_line <- sprintf("Wald test (nets = treatment): GDP $p$ = %s; EDU $p$ = %s.",
                     ifelse(is.na(p_wald_gdp),"NA",signif(p_wald_gdp,3)),
                     ifelse(is.na(p_wald_edu),"NA",signif(p_wald_edu,3)))
writeLines(c(
  "\\bottomrule",
  paste0("\\multicolumn{8}{p{0.95\\linewidth}}{\\footnotesize HC3 robust SEs. Outcome years matched to nearest treatment year. ",
         wald_line, "}"),
  "\\end{tabular}\\end{table}"
), con)
close(con)

cat("\nSaved main effects:\n- ", csv_main, "\n- ", tex_main, "\n", sep = "")

# =========================
# HETEROGENEITY BY INCOME
# =========================

# Baseline income: pre-2005 mean ln(GDPpc)
baseline <- panel %>%
  filter(year %in% 2000:2004) %>%
  group_by(country) %>%
  summarise(baseline_lngdp = if (all(is.na(ln_gdp_pc))) NA_real_ else mean(ln_gdp_pc, na.rm = TRUE),
            .groups = "drop")

# Merge into matched samples
gdp_h <- gdp_s %>% left_join(baseline, by = "country") %>% filter(!is.na(baseline_lngdp))
edu_h <- edu_s %>% left_join(baseline, by = "country") %>% filter(!is.na(baseline_lngdp))

# ---------- Continuous moderation ----------
gdp_c <- gdp_h %>% mutate(base_z = zfun(baseline_lngdp))
edu_c <- edu_h %>% mutate(base_z = zfun(baseline_lngdp))

fit_cont <- function(df) {
  keep <- c("z_nets","z_treat")
  have <- keep[colSums(!is.na(df[, keep])) > 0]
  if (!length(have) || nrow(df) < 5) return(NULL)
  rhs  <- paste(have, collapse = " + ")
  inter<- paste(paste0(have, ":base_z"), collapse = " + ")
  fml  <- as.formula(paste0("z_y ~ ", rhs, " + base_z + ", inter))
  m    <- lm(fml, data = df)
  list(m = m, df = df)
}

cont_gdp <- fit_cont(gdp_c)
cont_edu <- fit_cont(edu_c)

marg_effects_cont <- function(obj, label) {
  if (is.null(obj)) return(tibble())
  b  <- coef(obj$m)
  V  <- sandwich::vcovCL(obj$m, cluster = obj$df$country, type = "HC1")
  at_low  <- pctl(obj$df$base_z, 0.10)
  at_high <- pctl(obj$df$base_z, 0.90)
  
  one <- function(var) {
    if (!(var %in% names(b))) return(NULL)
    var_int <- paste0(var, ":base_z")
    lin <- function(x){
      R <- setNames(rep(0, length(b)), names(b)); R[var] <- 1
      if (var_int %in% names(b)) R[var_int] <- x
      est <- sum(R*b); se <- sqrt(drop(t(R)%*%V%*%R))
      tibble(Effect = est, SE = se, at = x)
    }
    rows <- bind_rows(lin(at_low), lin(at_high)) %>%
      mutate(Intervention = ifelse(var=="z_nets","Mosquito nets","Antimalarial treatment"),
             Outcome = label, Method = "OLS (clustered), continuous moderation")
    rows
  }
  bind_rows(one("z_nets"), one("z_treat")) %>%
    mutate(t = Effect/SE, p = 2*pnorm(-abs(t)),
           Effect = round(Effect,3), SE = round(SE,3),
           t = round(t,2), p = signif(p,3),
           at = ifelse(at==at_low,"Low (10th pct)","High (90th pct)")) %>%
    select(Outcome, Intervention, at, Effect, SE, t, p, Method)
}

cont_out <- bind_rows(
  marg_effects_cont(cont_gdp, "GDP (z)"),
  marg_effects_cont(cont_edu, "Education Index (z)")
)

# ---------- Group moderation (poor vs non-poor: bottom tercile) ----------
cut_terc <- pctl(baseline$baseline_lngdp, 1/3)
gdp_g <- gdp_h %>% mutate(poor = as.integer(baseline_lngdp <= cut_terc))
edu_g <- edu_h %>% mutate(poor = as.integer(baseline_lngdp <= cut_terc))

fit_group <- function(df) {
  keep <- c("z_nets","z_treat")
  have <- keep[colSums(!is.na(df[, keep])) > 0]
  if (!length(have) || nrow(df) < 5) return(NULL)
  inter <- paste(paste0(have, ":poor"), collapse = " + ")
  fml   <- as.formula(paste0("z_y ~ ", paste(have, collapse = " + "), " + poor + ", inter))
  m     <- lm(fml, data = df)
  list(m = m, df = df)
}

grp_gdp <- fit_group(gdp_g)
grp_edu <- fit_group(edu_g)

marg_effects_group <- function(obj, label) {
  if (is.null(obj)) return(tibble())
  b  <- coef(obj$m)
  V  <- sandwich::vcovCL(obj$m, cluster = obj$df$country, type = "HC1")
  
  one <- function(var) {
    if (!(var %in% names(b))) return(NULL)
    int <- paste0(var, ":poor")
    
    # Non-poor
    R0 <- setNames(rep(0, length(b)), names(b)); R0[var] <- 1
    est0 <- sum(R0*b); se0 <- sqrt(drop(t(R0)%*%V%*%R0))
    # Poor
    R1 <- R0; if (int %in% names(b)) R1[int] <- 1
    est1 <- sum(R1*b); se1 <- sqrt(drop(t(R1)%*%V%*%R1))
    # Wald (difference)
    Rd <- R1 - R0
    pd <- try({
      zd <- sum(Rd*b)/sqrt(drop(t(Rd)%*%V%*%Rd))
      2*pnorm(-abs(zd))
    }, silent = TRUE)
    if (inherits(pd,"try-error")) pd <- NA_real_
    
    tibble(
      Outcome = label,
      Intervention = ifelse(var=="z_nets","Mosquito nets","Antimalarial treatment"),
      Group = c("Non-poor","Poor (bottom tercile)"),
      Effect = round(c(est0, est1),3),
      SE     = round(c(se0, se1),3),
      t      = round(c(est0/se0, est1/se1),2),
      p      = signif(c(2*pnorm(-abs(est0/se0)), 2*pnorm(-abs(est1/se1))),3),
      Wald_p_diff = c(NA, ifelse(is.na(pd),"", signif(pd,3))),
      Method = "OLS (clustered), poor vs non-poor"
    )
  }
  bind_rows(one("z_nets"), one("z_treat"))
}

grp_out <- bind_rows(
  marg_effects_group(grp_gdp, "GDP (z)"),
  marg_effects_group(grp_edu, "Education Index (z)")
)

# ---------- Save heterogeneity outputs ----------
csv_cont <- file.path(outdir, "heterogeneity_continuous.csv")
csv_grp  <- file.path(outdir, "heterogeneity_groups.csv")
readr::write_csv(cont_out, csv_cont)
readr::write_csv(grp_out,  csv_grp)

tex_grp <- file.path(outdir, "heterogeneity_groups.tex")
con <- file(tex_grp, "wt", encoding = "UTF-8")
writeLines(c(
  "\\begin{table}[htbp]\\centering",
  "\\caption{Heterogeneous effects by income: Poor vs non-poor (baseline 2000--2004), nearest-year matched, clustered SEs}",
  "\\label{tab:hetero_groups}",
  "\\small",
  "\\begin{tabular}{@{}l l l c c c c c@{}}\\toprule",
  "Outcome & Intervention & Group & Effect & SE & $t$ & $p$ & Wald $p$ (diff) \\\\ \\midrule"
), con)
if (nrow(grp_out) > 0) {
  apply(grp_out, 1, function(r) {
    writeLines(sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                       r[["Outcome"]], r[["Intervention"]], r[["Group"]],
                       r[["Effect"]], r[["SE"]], r[["t"]], r[["p"]],
                       ifelse(is.na(r[["Wald_p_diff"]]) || r[["Wald_p_diff"]]=="" , "", r[["Wald_p_diff"]])), con)
  })
} else {
  writeLines("\\multicolumn{8}{c}{No estimable group effects.}", con)
}
writeLines(c(
  "\\bottomrule",
  "\\multicolumn{8}{p{0.95\\linewidth}}{\\footnotesize Baseline income is mean ln(GDPpc) over 2000--2004.",
  " SEs are clustered by country (HC1). Wald $p$ tests equality of effects between groups.}",
  "\\end{tabular}\\end{table}"
), con)
close(con)

cat("\nSaved heterogeneity:\n- ", csv_cont, "\n- ", csv_grp, "\n- ", tex_grp, "\n", sep = "")
