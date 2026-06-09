# robustness + exploratory models
# run after main script: dM, A_M, G_M, priors_main in memory

library(posterior)

# merge focal year + complexity components (by SCCS.Name)
extra <- read_csv("magicdolls_data_corrected.csv", show_col_types = FALSE) %>%
  select(SCCS.Name, focal_year,
         cc_writing, cc_fixity, cc_agriculture, cc_urbanization, cc_techspec,
         cc_landtransport, cc_money, cc_popdensity, cc_polinteg, cc_socstrat)
dM <- dM %>% left_join(extra, by = "SCCS.Name")

# standardise by 2 sd
z2 <- function(x) as.numeric((x - mean(x)) / (2*sd(x)))
cc <- c("cc_writing","cc_fixity","cc_agriculture","cc_urbanization","cc_techspec",
        "cc_landtransport","cc_money","cc_popdensity","cc_polinteg","cc_socstrat")
dM <- dM %>%
  mutate(Paragraphs_z = z2(log1p(Paragraphs)),
         Focal_year_z = z2(focal_year),
         across(all_of(cc), z2, .names = "{.col}_z"))

# summary helper (OR, 95% CrI, pd)
summ <- function(fit, term){
  x <- as_draws_df(fit)[[paste0("b_", term)]]
  data.frame(term = term,
             OR = round(median(exp(x)), 2),
             lo = round(quantile(exp(x), .025), 2),
             hi = round(quantile(exp(x), .975), 2),
             pd = round(max(mean(x > 0), mean(x < 0)), 3))
}

# fit helper (phylo + space, same as M2)
fit_ps <- function(form) brm(
  form, data = dM, data2 = list(A_M = A_M, G_M = G_M),
  family = bernoulli(), prior = priors_main,
  chains=4, cores=4, iter=4000,
  control=list(adapt_delta=0.99, max_treedepth=12),
  seed=123
)


# coverage control
fit_cov <- fit_ps(Magic ~ 1 + Complexity_z + Paragraphs_z +
                    (1 | gr(phylo, cov = A_M)) + (1 | gr(space, cov = G_M)))

# focal year control
fit_year <- fit_ps(Magic ~ 1 + Complexity_z + Focal_year_z +
                     (1 | gr(phylo, cov = A_M)) + (1 | gr(space, cov = G_M)))

# both
fit_both <- fit_ps(Magic ~ 1 + Complexity_z + Paragraphs_z + Focal_year_z +
                     (1 | gr(phylo, cov = A_M)) + (1 | gr(space, cov = G_M)))

# robustness table
rbind(cbind(model = "coverage", summ(fit_cov,  "Complexity_z")),
      cbind(model = "coverage", summ(fit_cov,  "Paragraphs_z")),
      cbind(model = "year",     summ(fit_year, "Complexity_z")),
      cbind(model = "year",     summ(fit_year, "Focal_year_z")),
      cbind(model = "both",     summ(fit_both, "Complexity_z")),
      cbind(model = "both",     summ(fit_both, "Paragraphs_z")),
      cbind(model = "both",     summ(fit_both, "Focal_year_z")))


# subtype models
subs <- c("Fertility","Protection","Punishment","Spirits")
fsub <- function(y) as.formula(paste0(
  y, " ~ 1 + Complexity_z + (1 | gr(phylo, cov = A_M)) + (1 | gr(space, cov = G_M))"))
fit_sub <- lapply(subs, function(y) fit_ps(fsub(y)))
names(fit_sub) <- subs

# subtype table
do.call(rbind, lapply(subs, function(y)
  cbind(subtype = y, summ(fit_sub[[y]], "Complexity_z"))))


# complexity component models (collinear, marginal not independent)
fcc <- function(v) as.formula(paste0(
  "Magic ~ 1 + ", v, "_z + (1 | gr(phylo, cov = A_M)) + (1 | gr(space, cov = G_M))"))
fit_cc <- lapply(cc, function(v) fit_ps(fcc(v)))
names(fit_cc) <- cc

# component table
do.call(rbind, lapply(cc, function(v)
  cbind(component = v, summ(fit_cc[[v]], paste0(v, "_z")))))


# phylo + spatial sds
sd_summ <- function(fit){
  dr <- as_draws_df(fit)
  data.frame(phylo = round(median(dr$sd_phylo__Intercept), 2),
             space = round(median(dr$sd_space__Intercept), 2))
}
rbind(magic = sd_summ(fit_magic_M2), play = sd_summ(fit_play_M2))

# play vs magic prevalence (pp difference)
pp_play  <- rowMeans(posterior_epred(fit_play_M2,  re_formula = NULL))
pp_magic <- rowMeans(posterior_epred(fit_magic_M2, re_formula = NULL))
diff_pp  <- 100 * (pp_play - pp_magic)
c(prob = mean(diff_pp > 0), med = median(diff_pp),
  lo = quantile(diff_pp, .025), hi = quantile(diff_pp, .975))
