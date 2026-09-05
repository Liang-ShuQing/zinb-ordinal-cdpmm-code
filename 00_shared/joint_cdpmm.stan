// Joint ZINB + ordinal with truncated stick-breaking CDPMM (no PG).
// Matches paper links: logit zero-inflation, logit NB mean, cumulative logit ordinal.
// Scenario 2: pz=4, pc=5, po=3; C=5; G=8; q=3.
// Cluster labels marginalized via log_sum_exp.

data {
  int<lower=1> N;
  int<lower=1> Ntot;
  int<lower=2> C;
  int<lower=2> G;
  int<lower=1> pz;
  int<lower=1> pc;
  int<lower=1> po;
  array[N] int<lower=1> nis;
  array[Ntot] int<lower=1, upper=N> id;
  array[Ntot] int<lower=0> y1;
  array[Ntot] int<lower=1, upper=C> y2;
  matrix[Ntot, pz] Xz;
  matrix[Ntot, pc] Xc;
  matrix[Ntot, po] Xo;
  real<lower=0> prior_sd;
  real<lower=0> tau_a1;
  real<lower=0> tau_a2;
}

parameters {
  vector[pz] alpha;
  vector[pc] beta;
  vector[po] gamma;
  real<lower=0> r;
  positive_ordered[C - 2] delta_raw;

  real<lower=0> tau;
  vector<lower=0, upper=1>[G - 1] nu;
  array[G] vector[3] mu_star;
  array[G] cholesky_factor_corr[3] L_Omega;
  array[G] vector<lower=0>[3] sigma_Omega;
  array[N] vector[3] b;
}

transformed parameters {
  vector[C - 1] delta;
  vector[G] pi_w;
  array[G] vector[3] mu;
  vector[3] mu_bar;
  delta[1] = 0;
  for (c in 2:(C - 1)) delta[c] = delta_raw[c - 1];

  {
    real rest = 1;
    for (g in 1:(G - 1)) {
      pi_w[g] = rest * nu[g];
      rest *= (1 - nu[g]);
    }
    pi_w[G] = rest;
  }
  mu_bar = rep_vector(0, 3);
  for (g in 1:G) mu_bar += pi_w[g] * mu_star[g];
  for (g in 1:G) mu[g] = mu_star[g] - mu_bar;
}

model {
  alpha ~ normal(0, prior_sd);
  beta ~ normal(0, prior_sd);
  gamma ~ normal(0, prior_sd);
  r ~ gamma(0.01, 0.01);
  delta_raw ~ normal(0, 2);
  tau ~ gamma(tau_a1, tau_a2);
  for (g in 1:(G - 1)) nu[g] ~ beta(1, tau);
  for (g in 1:G) {
    mu_star[g] ~ normal(0, 2);
    L_Omega[g] ~ lkj_corr_cholesky(2);
    sigma_Omega[g] ~ student_t(3, 0, 1);
  }

  for (i in 1:N) {
    vector[G] lps;
    for (g in 1:G) {
      matrix[3, 3] Lcov = diag_pre_multiply(sigma_Omega[g], L_Omega[g]);
      lps[g] = log(pi_w[g]) + multi_normal_cholesky_lpdf(b[i] | mu[g], Lcov);
    }
    target += log_sum_exp(lps);
  }

  for (n in 1:Ntot) {
    int i = id[n];
    real eta1 = Xz[n] * alpha + b[i][1];
    real eta2 = Xc[n] * beta + b[i][2];
    real eta3 = Xo[n] * gamma + b[i][3];
    real pi_at = inv_logit(eta1);
    real phi = inv_logit(eta2);
    real mu_nb = r * phi / fmax(1 - phi, 1e-12);

    if (y1[n] == 0) {
      target += log_sum_exp(
        log1m(pi_at),
        log(pi_at) + neg_binomial_2_lpmf(0 | mu_nb, r)
      );
    } else {
      target += log(pi_at) + neg_binomial_2_lpmf(y1[n] | mu_nb, r);
    }

    y2[n] ~ ordered_logistic(eta3, delta);
  }
}

generated quantities {
  matrix[3, 3] Sigma;
  vector[3] rho;
  Sigma = rep_matrix(0, 3, 3);
  for (g in 1:G) {
    matrix[3, 3] Om = multiply_lower_tri_self_transpose(
      diag_pre_multiply(sigma_Omega[g], L_Omega[g])
    );
    Sigma += pi_w[g] * (Om + mu[g] * mu[g]');
  }
  rho[1] = Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]);
  rho[2] = Sigma[1, 3] / sqrt(Sigma[1, 1] * Sigma[3, 3]);
  rho[3] = Sigma[2, 3] / sqrt(Sigma[2, 2] * Sigma[3, 3]);
}
