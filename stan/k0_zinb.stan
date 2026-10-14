// k0_zinb.stan

data {
  int<lower=1> N;                             // females
  int<lower=1> K;                             // treatments
  int<lower=1> Ncell;                         // treatment-by-set cells
  int<lower=1> S;                             // sires
  int<lower=1> F;                             // dams (full-sib families)
  array[N] int<lower=0> y;
  array[N] int<lower=1, upper=K>     tr;
  array[N] int<lower=1, upper=Ncell> cell;
  array[N] int<lower=1, upper=S>     sire;
  array[N] int<lower=1, upper=F>     dam;

  real<lower=0> lkj_eta;                      // LKJ shape on Omega
  real<lower=0> pr_sd_mu_s;                   // half-t(3, 0, .) scales
  real<lower=0> pr_sd_zi_s;
  real<lower=0> pr_sd_mu_d;
  real<lower=0> pr_sd_zi_d;
}
transformed data {
  array[K] int nT = rep_array(0, K);          // females per treatment
  for (n in 1:N) nT[tr[n]] += 1;

  // Split the zeros so positive counts can be vectorised
  array[N] int is0;
  int n0 = 0;
  for (n in 1:N) { is0[n] = (y[n] == 0); n0 += is0[n]; }
  array[n0]   int idx0;
  array[N-n0] int idx1;
  { int a = 1; int b = 1;
    for (n in 1:N) if (is0[n]) { idx0[a] = n; a += 1; } else { idx1[b] = n; b += 1; } }
}
parameters {
  vector[Ncell] beta;                         // log mu, cell means
  vector[K]     psi;                          // logit zi intercepts
  vector[K]     log_phi;                      // log NB shape
  vector<lower=0>[K] sd_mu_s, sd_zi_s, sd_mu_d, sd_zi_d;
  cholesky_factor_corr[K] L_Om;               // between-treatment, shared
  matrix[K, S] Z1_s, Z2_s;
  matrix[K, F] Z1_d, Z2_d;
}
transformed parameters {
  matrix[K, S] a_mu_s = diag_pre_multiply(sd_mu_s, L_Om * Z1_s);
  matrix[K, S] a_zi_s = diag_pre_multiply(sd_zi_s, L_Om * Z2_s);
  matrix[K, F] a_mu_d = diag_pre_multiply(sd_mu_d, Z1_d);
  matrix[K, F] a_zi_d = diag_pre_multiply(sd_zi_d, Z2_d);
}
model {
  vector[N] eta;                              // log mu
  vector[N] zeta;                             // logit zi
  vector[N] phi;
  for (n in 1:N) {
    int t = tr[n];
    eta[n]  = beta[cell[n]] + a_mu_s[t, sire[n]] + a_mu_d[t, dam[n]];
    zeta[n] = psi[t]        + a_zi_s[t, sire[n]] + a_zi_d[t, dam[n]];
    phi[n]  = exp(log_phi[t]);
  }

  psi     ~ logistic(0, 1);
  sd_mu_s ~ student_t(3, 0, pr_sd_mu_s);
  sd_zi_s ~ student_t(3, 0, pr_sd_zi_s);
  sd_mu_d ~ student_t(3, 0, pr_sd_mu_d);
  sd_zi_d ~ student_t(3, 0, pr_sd_zi_d);
  L_Om    ~ lkj_corr_cholesky(lkj_eta);
  to_vector(Z1_s) ~ std_normal();  to_vector(Z2_s) ~ std_normal();
  to_vector(Z1_d) ~ std_normal();  to_vector(Z2_d) ~ std_normal();
  // beta and log_phi: flat.

  // positive counts
  target += sum(log1m_inv_logit(zeta[idx1]))
          + neg_binomial_2_log_lpmf(y[idx1] | eta[idx1], phi[idx1]);
  // zeros: structural, or a sampling zero from the negative binomial
  for (i in 1:n0) {
    int n = idx0[i];
    target += log_sum_exp(log_inv_logit(zeta[n]),
                          log1m_inv_logit(zeta[n])
                            + neg_binomial_2_log_lpmf(0 | eta[n], phi[n]));
  }
}
generated quantities {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(L_Om);
  vector[N] log_lik;

  // Compact posterior predictive summaries
  array[K] int rep_zeros = rep_array(0, K);
  vector[K]    rep_mean  = rep_vector(0.0, K);
  array[S, K]  int sire_zeros = rep_array(0, S, K);

  for (n in 1:N) {
    int t = tr[n];
    real eta  = beta[cell[n]] + a_mu_s[t, sire[n]] + a_mu_d[t, dam[n]];
    real zeta = psi[t]        + a_zi_s[t, sire[n]] + a_zi_d[t, dam[n]];
    real phi  = exp(log_phi[t]);
    log_lik[n] = y[n] == 0
      ? log_sum_exp(log_inv_logit(zeta),
                    log1m_inv_logit(zeta) + neg_binomial_2_log_lpmf(0 | eta, phi))
      : log1m_inv_logit(zeta) + neg_binomial_2_log_lpmf(y[n] | eta, phi);

    // Replicate. eta is clamped with fmin only for the RNG
    int yr = bernoulli_logit_rng(zeta) == 1
               ? 0 : neg_binomial_2_log_rng(fmin(eta, 12.0), phi);
    if (yr == 0) { rep_zeros[t] += 1; sire_zeros[sire[n], t] += 1; }
    rep_mean[t] += yr;
  }
  for (t in 1:K) rep_mean[t] /= nT[t];
}
