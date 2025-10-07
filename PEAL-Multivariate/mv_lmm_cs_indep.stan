data {
  int<lower=1> N;              // observations
  int<lower=1> H;              // sites
  int<lower=1> P;              // patients
  int<lower=1> R;              // outcomes
  int<lower=1> p;              // predictors
  matrix[N, p] X;
  matrix[N, R] y;
  int<lower=1, upper=H> site_id[N];
  int<lower=1, upper=P> pat_id[N];
  int<lower=1, upper=H> pat_site[P]; // site index for each patient
}

parameters {
  matrix[p, R] beta;                 // fixed effects
  vector[H] u_raw;                   // site random intercepts (raw)
  vector[P] v_raw;                   // patient random intercepts (raw)

  real<lower=0> sigma_u;             // site RI SD
  vector<lower=0>[H] sigma_v_h;      // patient RI SD, one per site

  real<lower=0> sigma_e;             // residual SD (shared over outcomes)
}

transformed parameters {
  vector[H] u = sigma_u * u_raw;
  vector[P] v;
  for (p_i in 1:P) {
    v[p_i] = sigma_v_h[ pat_site[p_i] ] * v_raw[p_i];
  }
}

model {
  // Priors
  to_vector(beta) ~ normal(0, 2);
  u_raw ~ normal(0, 1);
  v_raw ~ normal(0, 1);

  sigma_u  ~ student_t(3, 0, 1);
  sigma_v_h ~ student_t(3, 0, 1);
  sigma_e  ~ student_t(3, 0, 1);
  // weakly informative prior on rho
  // (implicitly uniform on the constrained interval)

  // Build exchangeable covariance once per leapfrog step
  matrix[R, R] J = rep_matrix(1.0, R, R);
  matrix[R, R] I = diag_matrix(rep_vector(1.0, R));
  matrix[R, R] Sigma = square(sigma_e) * ((1 - 0) * I + 0 * J);

  // Likelihood
  for (n in 1:N) {
    row_vector[R] mu_n = (row(X, n) * beta) + rep_row_vector(u[site_id[n]] + v[pat_id[n]], R);
    y[n] ~ multi_normal(mu_n', Sigma);
  }
}

generated quantities {
  // For convenience: convert to variances and correlation
  real sigma2_u = square(sigma_u);
  vector[H] sigma2_v_h;
  for (h in 1:H) sigma2_v_h[h] = square(sigma_v_h[h]);
  real sigma2_e = square(sigma_e);
}
