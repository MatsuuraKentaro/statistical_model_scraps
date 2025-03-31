functions {
  real binormal_cdf_func(real z1, real z2, real rho) {
    if (z1 == 0 && z2 == 0) {
      return 0.25 + asin(rho) / (2 * pi());
    }
    real z1p = z1 == 0 ? z1 + 1E-20 : z1;
    real z2p = z2 == 0 ? z2 + 1E-20 : z2;
    real a1 = (z2 / z1p - rho) / sqrt(1 - rho^2);
    real a2 = (z1 / z2p - rho) / sqrt(1 - rho^2);
    real delta = z1 * z2 < 0 || (z1 * z2 == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2);
  }
}

data {
  int N_om; // the number of subjects in group om (Y1, Y2) = (observed, missing)
  int N_mo; // the number of subjects in group mo (Y1, Y2) = (missing, observed)
  int N_oo; // the number of subjects in group oo (Y1, Y2) = (observed, observed)
  int Y_om_1x; // the number of (Y1, Y2) = (1, missing) in group om
  int Y_mo_x1; // the number of (Y1, Y2) = (missing, 1) in group mo
  array[4] int Y_oo;  // the numbers of (Y1,Y2) = (0,0),(1,0),(0,1),(1,1) in group oo
}

parameters {
  vector[2] x;
  real<lower=-1, upper=1> rho;
}

transformed parameters {
  real p_om_1x = 1 - binormal_cdf_func(x[1], 100, rho);
  real p_mo_x1 = 1 - binormal_cdf_func(100, x[2], rho);
  real p_oo_00 = binormal_cdf_func(x[1], x[2], rho);
  real p_oo_10 = (1 - p_mo_x1) - p_oo_00;
  real p_oo_01 = (1 - p_om_1x) - p_oo_00;
  real p_oo_11 = 1 - (p_oo_00 + p_oo_10 + p_oo_01);
}

model {
  Y_om_1x ~ binomial(N_om, p_om_1x);
  Y_mo_x1 ~ binomial(N_mo, p_mo_x1);
  Y_oo ~ multinomial([p_oo_00, p_oo_10, p_oo_01, p_oo_11]');
}
