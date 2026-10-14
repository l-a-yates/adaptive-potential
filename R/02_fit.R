## 02_fit.R: fit model k0 with CmdStan and save results/data/k0.rds (skipped if it exists).

source(here::here("R", "00_functions.R"))

if (file.exists(FIT_PATH)) {
  message("using cached fit ", FIT_PATH, "; delete to resample")
} else {

  library(cmdstanr)   # from the renv library; CmdStan itself is per-machine (~/.cmdstan)

  fitdat <- readRDS(out("fitdat.rds"))

  ## Stan CSVs
  CSV_DIR <- path.expand("~/.cache/adaptive-potential-stan/k0")
  EXE_DIR <- file.path(CSV_DIR, "exe")
  dir.create(EXE_DIR, recursive = TRUE, showWarnings = FALSE)

  ## Priors: student_t(3, 0, 3.4) on every SD (brms default), lkj(1) on Omega;
  ## logistic(0, 1) on the zi intercepts and flat cell means/shapes are in the Stan file.
  sdat <- list(
    N = nrow(fitdat), K = length(TREATS), Ncell = nlevels(fitdat$SET),
    S = nlevels(fitdat$SIRE), F = nlevels(fitdat$DAM),
    y    = fitdat$TOTAL,
    tr   = as.integer(fitdat$TREAT),
    cell = as.integer(fitdat$SET),
    sire = as.integer(fitdat$SIRE),
    dam  = as.integer(fitdat$DAM),
    lkj_eta = 1,
    pr_sd_mu_s = 3.4, pr_sd_zi_s = 3.4, pr_sd_mu_d = 3.4, pr_sd_zi_d = 3.4)

  mod <- cmdstan_model(here::here("stan", "k0_zinb.stan"), dir = EXE_DIR)
  f <- mod$sample(
    data = sdat, chains = N_CHAINS, parallel_chains = N_CHAINS, seed = SEED,
    output_dir = CSV_DIR, iter_warmup = N_WARMUP, iter_sampling = N_SAMPLES,
    adapt_delta = 0.95, max_treedepth = 12, refresh = 250)

  ## Hyperparameter draws, chains merged in order
  vars <- c("beta", "psi", "log_phi", "sd_mu_s", "sd_zi_s", "sd_mu_d", "sd_zi_d", "Omega")
  d <- as.data.frame(posterior::merge_chains(f$draws(variables = vars, format = "draws_df")))

  ## log_lik (~0.9 GB) is used for loo and then dropped.
  ll <- f$draws("log_lik")
  lo <- suppressWarnings(loo::loo(ll, cores = 2))
  rm(ll)

  ppc_vars <- c("rep_zeros", "rep_mean", "sire_zeros")
  ppc <- setNames(lapply(ppc_vars, function(v) posterior::as_draws_matrix(f$draws(v))), ppc_vars)

  ## Sire breeding values for relative fitness (Fig. 3): a_mu[t,s] + alpha_t * a_zi[t,s],
  ## per draw, summarised per sire and treatment.
  alpha <- decompose_va_w(d)$alpha
  a_mu  <- as.matrix(posterior::as_draws_matrix(posterior::merge_chains(f$draws("a_mu_s"))))
  a_zi  <- as.matrix(posterior::as_draws_matrix(posterior::merge_chains(f$draws("a_zi_s"))))
  sires <- levels(fitdat$SIRE)
  obs   <- paste(fitdat$SIRE, fitdat$TREAT)
  sire_bv <- bind_rows(lapply(seq_along(TREATS), function(k) {
    cm <- sprintf("a_mu_s[%d,%d]", k, seq_along(sires))
    cz <- sprintf("a_zi_s[%d,%d]", k, seq_along(sires))
    bv <- a_mu[, cm] + alpha[, k] * a_zi[, cz]
    data.frame(SIRE = sires, TREAT = TREATS[k],
               bv_mean = colMeans(bv), bv_sd = apply(bv, 2, sd),
               mu_mean = colMeans(a_mu[, cm]), zi_mean = colMeans(a_zi[, cz]),
               observed = paste(sires, TREATS[k]) %in% obs)
  }))

  res <- list(
    loo     = lo,
    ppc     = ppc,
    draws   = d,
    sire_bv = sire_bv,
    diag    = f$diagnostic_summary(),
    summary = f$summary(c("psi", "log_phi", "sd_mu_s", "sd_zi_s", "sd_mu_d", "sd_zi_d")))
  saveRDS(res, FIT_PATH)

  cat(sprintf("saved %s | divergences %d | max treedepth hits %d\n", FIT_PATH,
              sum(res$diag$num_divergent), sum(res$diag$num_max_treedepth)))
}
