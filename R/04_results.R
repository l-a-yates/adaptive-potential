## 04_results.R: mean fitness, V_A(w) and its decomposition, diet contrasts,
## r_A and the margin contrasts.

source(here::here("R", "00_functions.R"))

fitdat <- readRDS(out("fitdat.rds"))
draws  <- readRDS(FIT_PATH)$draws

d    <- decompose_va_w(draws)
Wbar <- wbar_from_draws(draws, fitdat, d)

saveRDS(list(d = d, Wbar = Wbar), out("va.rds"))


## ---- mean fitness -----------------------------------------------------------
## Offspring per day of laying, with sire, dam and zi effects marginalised and
## assay sets averaged; raw_per_day is the unmodelled sample mean.
wbar_tab <- bind_rows(lapply(TREATS, function(tt) {
  raw <- mean(fitdat$TOTAL[fitdat$TREAT == tt]) / LAYDAYS[[tt]]
  cbind(data.frame(treat = tt, temp = TEMP_OF[[tt]], diet = DIET_LAB[[DIET_OF[[tt]]]]),
        smry(Wbar[, tt]),
        data.frame(raw_per_day = round(raw, 3),
                   estimate = fmt(Wbar[, tt]),
                   log_estimate = fmt(log(Wbar[, tt]))))
}))
print(wbar_tab)
write.csv(wbar_tab, out("tab_mean_fitness.csv"), row.names = FALSE)


## ---- V_A(w) -----------------------------------------------------------------
s_mu <- par_mat(draws, "sd_mu_s"); s_zi <- par_mat(draws, "sd_zi_s")
va_tab <- bind_rows(lapply(seq_along(TREATS), function(i) {
  tt <- TREATS[i]
  data.frame(treat = tt, temp = TEMP_OF[[tt]], diet = DIET_LAB[[DIET_OF[[tt]]]],
             q = fmt(d$q[, i]), alpha = fmt(d$alpha[, i]),
             sd_mu = fmt(s_mu[, i]), sd_zi = fmt(s_zi[, i]),
             count = fmt(d$count[, i]), zi = fmt(d$zi[, i]), VA_w = fmt(d$VA_w[, i]),
             ## ratio of posterior means; the per-draw ratio blows up when zi is near 0
             zi_pct = sprintf("%.0f%%", 100 * mean(d$zi[, i]) / mean(d$VA_w[, i])))
}))
print(va_tab)
write.csv(va_tab, out("tab_va_per_treatment.csv"), row.names = FALSE)


## ---- diet effects and temperature x diet interactions on V_A(w) -------------
con_tab <- summarise_draws_col(va_contrasts(d$VA_w)) %>%
  transmute(quantity, VA_w = estimate, p_gt0 = round(p_gt0, 3))
print(con_tab)
write.csv(con_tab, out("tab_contrasts.csv"), row.names = FALSE)


## ---- genetic correlations in fitness ----------------------------------------
R_w <- cor_from_G(d$G_w)
Rbar <- apply(R_w, 2:3, mean)
print(round(Rbar, 3))
write.csv(round(Rbar, 4), out("tab_rA_matrix.csv"))

mc_tab <- summarise_draws_col(margin_contrasts(R_w)) %>%
  transmute(margin, contrast, definition, r_A = estimate, p_gt0 = round(p_gt0, 3))
print(mc_tab)
write.csv(mc_tab, out("tab_margin_contrasts.csv"), row.names = FALSE)

