## 03_diagnostics.R: convergence (Rhat, ESS, sampler), posterior predictive checks,
## trace plots and PSIS-LOO.

library(ggplot2)

source(here::here("R", "00_functions.R"))
theme_set(theme_classic())

fitdat <- readRDS(out("fitdat.rds"))
r <- readRDS(FIT_PATH)


## ---- convergence ------------------------------------------------------------
## Rhat and ESS for every parameter
pars <- setdiff(names(r$draws), c(".chain", ".iteration", ".draw"))
## keep only the upper triangle of Omega
dup  <- grepl("^Omega\\[", pars) &
        vapply(strsplit(gsub("^Omega\\[|\\]$", "", pars), ","), function(ij)
          length(ij) == 2 && as.integer(ij[1]) >= as.integer(ij[2]), logical(1))
pars <- pars[!dup]

cv_tab <- bind_rows(lapply(pars, function(p) {
  m <- as_chains(r$draws[[p]])
  data.frame(variable = p, rhat = posterior::rhat(m),
             ess_bulk = posterior::ess_bulk(m), ess_tail = posterior::ess_tail(m))
}))

conv <- data.frame(
  n_monitored        = nrow(cv_tab),
  max_rhat           = round(max(cv_tab$rhat, na.rm = TRUE), 4),
  n_rhat_above_101   = sum(cv_tab$rhat > 1.01, na.rm = TRUE),
  min_ess_bulk       = round(min(cv_tab$ess_bulk, na.rm = TRUE)),
  divergences        = sum(r$diag$num_divergent),
  max_treedepth_hits = sum(r$diag$num_max_treedepth),
  min_ebfmi          = round(min(r$diag$ebfmi), 3),
  draws              = nrow(r$draws))

cat("\n=== convergence ===\n")
print(conv, row.names = FALSE)
write.csv(conv, out("tab_convergence.csv"), row.names = FALSE)
write.csv(cv_tab %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
          out("tab_convergence_by_parameter.csv"), row.names = FALSE)

cat("\n  worst five by Rhat:\n")
print(head(cv_tab[order(-cv_tab$rhat), ], 5), digits = 4, row.names = FALSE)
cat("\n  worst five by bulk ESS:\n")
print(head(cv_tab[order(cv_tab$ess_bulk), ], 5), digits = 4, row.names = FALSE)

cat("\n=== the zero-inflation SDs ===\n")
print(as.data.frame(r$summary[grep("^sd_zi", r$summary$variable),
                              c("variable", "median", "q5", "q95", "rhat", "ess_bulk")]),
      digits = 3, row.names = FALSE)


## ---- posterior predictive checks, by treatment ------------------------------
obs <- function(f) vapply(TREATS, function(tt) f(fitdat$TOTAL[fitdat$TREAT == tt]), numeric(1))
o_zero <- obs(function(y) mean(y == 0))
o_mean <- obs(mean)
nT     <- as.vector(table(fitdat$TREAT)[TREATS])

rep_zfrac <- sweep(r$ppc$rep_zeros, 2, nT, "/")
ppc_p <- function(R, o) vapply(seq_along(o), function(i) mean(R[, i] >= o[i]), numeric(1))

ppc_tab <- data.frame(
  treat    = TREATS,
  zero_obs = round(o_zero, 3), zero_rep = round(colMeans(rep_zfrac), 3),
  zero_p   = round(ppc_p(rep_zfrac, o_zero), 2),
  mean_obs = round(o_mean, 1), mean_rep = round(colMeans(r$ppc$rep_mean), 1),
  mean_p   = round(ppc_p(r$ppc$rep_mean, o_mean), 2))

cat("\n=== posterior predictive checks (p = Pr(replicate >= observed)) ===\n")
print(ppc_tab, row.names = FALSE)
write.csv(ppc_tab, out("tab_ppc.csv"), row.names = FALSE)


## ---- per-sire zero counts ---------------------------------------------------
sz    <- r$ppc$sire_zeros
ij    <- do.call(rbind, lapply(strsplit(gsub("^sire_zeros\\[|\\]$", "", colnames(sz)), ","),
                               as.integer))
sires <- levels(fitdat$SIRE)

sire_tab <- bind_rows(lapply(seq_along(TREATS), function(t) {
  cols <- which(ij[, 2] == t)
  cols <- cols[order(ij[cols, 1])]
  x    <- fitdat[fitdat$TREAT == TREATS[t], ]
  o    <- tapply(x$TOTAL == 0, factor(as.character(x$SIRE), levels = sires), sum)
  o[is.na(o)] <- 0
  rep_sd <- apply(sz[, cols, drop = FALSE], 1, sd)
  data.frame(treat = TREATS[t],
             obs_total = sum(o),   rep_total = round(mean(rowSums(sz[, cols, drop = FALSE])), 1),
             obs_sd    = round(sd(o), 3), rep_sd = round(mean(rep_sd), 3),
             p_sd      = round(mean(rep_sd >= sd(o)), 3))
}))

cat("\n=== per-sire zero counts: is failure clustered by sire as observed? ===\n")
cat("  obs_sd / rep_sd : SD of the per-sire zero count, observed vs replicated\n")
cat("  p_sd near 0 means the model UNDER-disperses failure across sires\n\n")
print(sire_tab, row.names = FALSE)
write.csv(sire_tab, out("tab_ppc_sire_zeros.csv"), row.names = FALSE)


## ---- figures ----------------------------------------------------------------
pdf(out("fig_ppc.pdf"), width = 10, height = 7)
checks <- list(`zero fraction`     = list(R = rep_zfrac,      o = o_zero),
               `mean productivity` = list(R = r$ppc$rep_mean, o = o_mean))
for (nm in names(checks)) {
  R <- checks[[nm]]$R; o <- checks[[nm]]$o
  dl <- tibble(treat = rep(TREATS, each = nrow(R)), value = as.vector(R))
  print(ggplot(dl, aes(value)) +
          geom_histogram(bins = 40, fill = "grey60") +
          geom_vline(data = tibble(treat = TREATS, o = o), aes(xintercept = o),
                     colour = "#D95F02", linewidth = 0.8) +
          facet_wrap(~ factor(treat, TREATS), scales = "free") +
          labs(x = nm, y = NULL,
               title = paste0("Replicated ", nm, "; orange = observed")))
}
dev.off()

pdf(out("fig_trace.pdf"), width = 11, height = 8)
for (grp in c("sd_mu_s", "sd_zi_s", "sd_mu_d", "sd_zi_d")) {
  cols <- sprintf("%s[%d]", grp, seq_along(TREATS))
  dl <- bind_rows(lapply(seq_along(cols), function(i) {
    m <- as_chains(r$draws[[cols[i]]])
    tibble(treat = TREATS[i], iter = rep(seq_len(nrow(m)), N_CHAINS),
           chain = factor(rep(seq_len(N_CHAINS), each = nrow(m))),
           value = as.vector(m))
  }))
  print(ggplot(dl, aes(iter, value, colour = chain)) +
          geom_line(linewidth = 0.2, alpha = 0.7) +
          facet_wrap(~ factor(treat, TREATS), scales = "free_y") +
          labs(title = grp, x = "iteration", y = NULL))
}
dev.off()


## ---- PSIS-LOO ---------------------------------------------------------------
k <- r$loo$diagnostics$pareto_k
loo_tab <- data.frame(
  elpd_loo = round(r$loo$estimates["elpd_loo", 1], 1),
  se       = round(r$loo$estimates["elpd_loo", 2], 1),
  p_loo    = round(r$loo$estimates["p_loo", 1]),
  n_obs    = length(k),
  k_gt_0.7 = sum(k > 0.7, na.rm = TRUE),
  k_gt_1   = sum(k > 1, na.rm = TRUE),
  max_k    = round(max(k, na.rm = TRUE), 2))
cat("\n=== PSIS-LOO ===\n")
print(loo_tab, row.names = FALSE)
write.csv(loo_tab, out("tab_loo.csv"), row.names = FALSE)
