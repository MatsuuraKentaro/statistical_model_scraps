library(tidyverse)
library(cmdstanr)

d <- read_csv("data/BRAIN.csv") |> 
  arrange(weeks)

X <- d |> select(resect75)
D <- ncol(X)
data <- list(N = nrow(d), D = D, X = X, Event = d$event)

## Estimation
model_a <- cmdstan_model(stan_file = "stan/model1a.stan")
fit_a <- model_a$sample(data = data, seed = 123, parallel_chains = 4)

model_b <- cmdstan_model(stan_file = "stan/model1b.stan")
fit_b <- model_b$sample(data = data, seed = 123, parallel_chains = 4)
