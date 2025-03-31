library(mvtnorm)
library(cmdstanr)

set.seed(1234)

N_om <- 300
N_mo <- 300
N_oo <- 300
x <- c(0.25, 0.15)
rho <- 0.7
Omega <- matrix(c(1, rho, rho, 1), 2, 2)
p_om_1x_true <- pmvnorm(lower = c(x[1], -Inf), upper = c(Inf, Inf), corr = Omega)[1]
p_mo_x1_true <- pmvnorm(lower = c(-Inf, x[2]), upper = c(Inf, Inf), corr = Omega)[1]
p_oo_00_true <- pmvnorm(upper = x, corr = Omega)[1]
p_oo_10_true <- pmvnorm(lower = c(x[1], -Inf), upper = c(Inf, x[2]), corr = Omega)[1]
p_oo_01_true <- pmvnorm(lower = c(-Inf, x[2]), upper = c(x[1], Inf), corr = Omega)[1]
p_oo_11_true <- pmvnorm(lower = c(x[1], x[2]), upper = c(Inf, Inf), corr = Omega)[1]

Y_om_1x <- rbinom(n = 1, size = N_om, prob = p_om_1x_true)
Y_mo_x1 <- rbinom(n = 1, size = N_mo, prob = p_mo_x1_true)
Y_oo <- rmultinom(
  n = 1, size = N_oo, 
  prob = c(p_oo_00_true, p_oo_10_true, p_oo_01_true, p_oo_11_true)) |> 
  as.vector()

## estimation
data <- list(N_om = N_om, N_mo = N_mo, N_oo = N_oo,
             Y_om_1x = Y_om_1x, Y_mo_x1 = Y_mo_x1, Y_oo = Y_oo)
model <- cmdstan_model("stan/model1.stan")
fit <- model$sample(data = data, seed = 1234)

