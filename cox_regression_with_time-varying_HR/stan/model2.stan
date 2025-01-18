functions {
  vector gp(array[] real x, vector mu, vector eta,
            real a, real rho, real s2) {
    int N = size(x);
    matrix[N,N] k_addI = add_diag(gp_exp_quad_cov(x, a, rho), s2);
    vector[N] out = mu + cholesky_decompose(k_addI)*eta;
    return out;
  }
  
  vector gp_pred_rng(array[] real x, array[] real xp,
                     vector mu, vector mup, vector out,
                     real a, real rho, real s2, real sp2) {
    int N  = size(x);
    int Np = size(xp);
    matrix[N, N ] k_addI = add_diag(gp_exp_quad_cov(x, a, rho), s2);
    matrix[Np,Np] kpp_addI = add_diag(gp_exp_quad_cov(xp, a, rho), sp2);
    matrix[N, Np] kp = gp_exp_quad_cov(x, xp, a, rho);
    matrix[Np,N ] coef = kp' / k_addI;
    vector[Np] mup_cond = mup + coef * (out - mu);
    matrix[Np,Np] covp_cond = kpp_addI - coef * kp;
    return multi_normal_rng(mup_cond, covp_cond);
  }
}

data {
  int N;
  int Np;
  int D;
  matrix[N, D] X;
  array[N] real Time;
  array[Np] real Timep;
  vector[N] Event;
}

transformed data {
  vector[N] Mu = rep_vector(0, N);
  vector[Np] Mup = rep_vector(0, Np);
}

parameters {
  vector[N] eta;
  real<lower=0> a;
  real<lower=0> rho;
}

transformed parameters {
  matrix[N, D] b;
  for (d in 1:D) {
    b[:,d] = gp(Time, Mu, eta, a, rho, 1e-8);  
  }
}

model {
  vector[N] Xb = rows_dot_product(X, b);
  vector[N] log_cumsum = log(reverse(cumulative_sum(reverse(exp(Xb)))));
  vector[N] pl = Event .* (Xb - log_cumsum);
  target += sum(pl);
  
  eta ~ std_normal();
  a   ~ std_normal();
  rho ~ normal(0.1, 0.1);
  to_vector(b) ~ std_normal();
}

generated quantities {
  matrix[Np, D] bp;
  for (d in 1:D) {
    bp[:,d] = gp_pred_rng(Time, Timep, Mu, Mup, b[:,d], a, rho, 1e-8, 1e-8);  
  }
}
