## 00_functions.R: constants, paths, posterior summaries and the V_A(w) calculations

library(dplyr)

TREATS  <- c("T14DA", "T14DB", "T25DA", "T25DB", "T28DA", "T28DB")
TEMP_OF <- c(T14DA = 14, T14DB = 14, T25DA = 25, T25DB = 25, T28DA = 28, T28DB = 28)
DIET_OF <- c(T14DA = "DA", T14DB = "DB", T25DA = "DA", T25DB = "DB",
             T28DA = "DA", T28DB = "DB")
DIET_LAB <- c(DA = "Control diet", DB = "Low-calorie diet")
LAYDAYS <- c(T14DA = 11, T14DB = 15, T25DA = 6, T25DB = 8, T28DA = 6, T28DB = 9)

MARGINS <- c(Cold = 14, Hot = 28)

N_CHAINS  <- 4
N_WARMUP  <- 2000
N_SAMPLES <- 4000
SEED      <- 20261001

FIT_PATH <- here::here("results", "data", "k0.rds")

out <- function(name) {
  ext  <- tolower(tools::file_ext(name))
  kind <- if (ext %in% c("pdf", "png")) "plots" else if (ext == "rds") "data" else "tables"
  dir.create(here::here("results", kind), recursive = TRUE, showWarnings = FALSE)
  here::here("results", kind, name)
}

## Reshape a stacked vector of draws to
## iterations x chains for rhat() and ess_*().
as_chains <- function(v) matrix(v, ncol = N_CHAINS)

## D x k matrix of a vector parameter, e.g. par_mat(df, "psi").
par_mat <- function(df, par, k = length(TREATS))
  as.matrix(df[, sprintf("%s[%d]", par, seq_len(k))])

## Apply f(i, j) over the k x k treatment pairs, giving a D x k x k array.
pairs_array <- function(D, f) {
  k <- length(TREATS); A <- array(0, c(D, k, k))
  for (i in seq_len(k)) for (j in seq_len(k)) A[, i, j] <- f(i, j)
  A
}


## alpha_t = -E[pi(1-pi)] / E[1-pi] and q_t = E[pi], with pi = logistic(psi + z),
## z ~ N(0, sigma_z^2). Quadrature on an equally spaced grid over +/- 8 SD,
## weighted by the normal density
zi_gradient <- function(psi, sigma_z, n = 129L, lim = 8) {
  x <- seq(-lim, lim, length.out = n)
  w <- dnorm(x); w <- w / sum(w)
  p      <- plogis(psi + outer(sigma_z, x))        # D x n
  E_1mp  <- as.vector((1 - p) %*% w)
  E_p1mp <- as.vector((p * (1 - p)) %*% w)
  list(alpha = -E_p1mp / E_1mp, q = 1 - E_1mp)
}


## V_A(w) per draw. With v_mu = 2 s_mu and v_zi = 2 alpha s_zi (alpha integrates over
## sire plus dam zi variance), G_w = (v_mu v_mu' + v_zi v_zi') * Omega and its
## diagonal splits into count = v_mu^2 and zi = v_zi^2. All D x 6 except G_w (D x 6 x 6).
decompose_va_w <- function(df) {
  s_mu <- par_mat(df, "sd_mu_s"); s_zi <- par_mat(df, "sd_zi_s")
  g <- lapply(seq_along(TREATS), function(i)
    zi_gradient(df[[sprintf("psi[%d]", i)]], sqrt(s_zi[, i]^2 + df[[sprintf("sd_zi_d[%d]", i)]]^2)))
  alpha <- sapply(g, `[[`, "alpha"); q <- sapply(g, `[[`, "q")
  v_mu <- 2 * s_mu; v_zi <- 2 * alpha * s_zi
  G_w <- pairs_array(nrow(df), function(i, j)
    (v_mu[, i] * v_mu[, j] + v_zi[, i] * v_zi[, j]) * df[[sprintf("Omega[%d,%d]", i, j)]])
  nm <- function(m) { colnames(m) <- TREATS; m }
  list(alpha = nm(alpha), q = nm(q), count = nm(v_mu^2), zi = nm(v_zi^2),
       VA_w = nm(v_mu^2 + v_zi^2), G_w = G_w)
}

## D x k x k correlations from a D x k x k covariance array.
cor_from_G <- function(G) {
  s <- sapply(seq_along(TREATS), function(i) sqrt(G[, i, i]))
  R <- pairs_array(dim(G)[1], function(i, j) G[, i, j] / (s[, i] * s[, j]))
  dimnames(R) <- list(NULL, TREATS, TREATS)
  R
}


## Mean fitness per day of laying, random effects marginalised.
wbar_from_draws <- function(df, fitdat, d) {
  cells  <- levels(fitdat$SET)
  cell_n <- as.vector(table(fitdat$SET))
  cell_t <- as.character(fitdat$TREAT[match(cells, fitdat$SET)])
  mns    <- exp(par_mat(df, "beta", length(cells)))
  m_t <- sapply(TREATS, function(tt) {
    w <- cell_n * (cell_t == tt)
    as.vector(mns %*% (w / sum(w)))
  })
  sd2 <- par_mat(df, "sd_mu_s")^2 + par_mat(df, "sd_mu_d")^2
  Wbar <- (1 - d$q) * m_t * exp(0.5 * sd2) / rep(LAYDAYS[TREATS], each = nrow(df))
  colnames(Wbar) <- TREATS
  Wbar
}


## Diet effect on V_A(w) at each temperature, and the temperature x diet
## interaction at each thermal margin relative to 25 C, per draw.
va_contrasts <- function(VA) {
  de <- lapply(c(14, 25, 28), function(tm) VA[, sprintf("T%dDB", tm)] - VA[, sprintf("T%dDA", tm)])
  tibble(
    quantity = c("diet effect at 14C", "diet effect at 25C", "diet effect at 28C",
                 "interaction cold (14C - 25C)", "interaction hot (28C - 25C)"),
    draws = list(de[[1]], de[[2]], de[[3]], de[[1]] - de[[2]], de[[3]] - de[[2]]))
}

## Margin contrasts on a D x 6 x 6 correlation array R:
##   P1 = rho_T(control) - rho_T(low-cal),  P2 = rho_D(focal) - rho_D(25 C).
margin_contrasts <- function(R) {
  r <- function(a, b) R[, a, b]
  bind_rows(lapply(names(MARGINS), function(m) {
    f <- MARGINS[[m]]; t <- function(temp, diet) sprintf("T%d%s", temp, diet)
    tibble(
      margin   = m,
      contrast = c("P1", "P2"),
      definition = c(sprintf("rho_T(control) - rho_T(low-cal); 25 vs %d C", f),
                     sprintf("rho_D(%d C) - rho_D(25 C)", f)),
      draws = list(r(t(25, "DA"), t(f, "DA")) - r(t(25, "DB"), t(f, "DB")),
                   r(t(f, "DA"), t(f, "DB")) - r(t(25, "DA"), t(25, "DB"))))
  }))
}


## ---- posterior summaries: mean and 90% HPD ----------------------------------

hpd <- function(x) {
  h <- coda::HPDinterval(coda::as.mcmc(x), prob = 0.90)
  c(lower = h[1], upper = h[2])
}

smry <- function(x) {
  h <- hpd(x)
  tibble(mean = mean(x), lower = unname(h[1]), upper = unname(h[2]), p_gt0 = mean(x > 0))
}

fmt <- function(x) {
  h <- hpd(x)
  sprintf("%.3f (%.3f, %.3f)", mean(x), h[1], h[2])
}

## Summarise a list-column of draws, dropping the draws.
summarise_draws_col <- function(tab) {
  bind_cols(select(tab, -draws), bind_rows(lapply(tab$draws, smry))) %>%
    mutate(estimate = vapply(tab$draws, fmt, character(1)))
}
