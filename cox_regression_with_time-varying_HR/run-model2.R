library(tidyverse)
library(cmdstanr)

d <- read_csv("data/BRAIN.csv") |> 
  arrange(weeks)
  
X <- d |> select(resect75)
D <- ncol(X)
Np <- 101
Timep <- seq(0, 1, len = Np)
max_week <- 200
data <- list(N = nrow(d), Np = Np, D = D, X = X, 
             Time = d$weeks/max_week, Timep = Timep, Event = d$event)

## Estimation
model <- cmdstan_model(stan_file = "stan/model2.stan")
fit <- model$sample(data = data, seed = 123, parallel_chains = 4)


## Visualization
d_est <- fit$draws("bp", format = "matrix") |>
  apply(2, quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95)) |> t() |>
  data.frame(weeks = Timep * max_week, check.names = FALSE)

p <- ggplot(data = d_est) +
  theme(text = element_text(size = 20)) +
  geom_ribbon(aes(x = weeks, ymin = `5%`, ymax = `95%`), fill = "black", alpha = 1/6) +
  geom_ribbon(aes(x = weeks, ymin = `25%`, ymax = `75%`), fill = "black", alpha = 2/6) +
  geom_line(aes(x = weeks, y = `50%`)) +
  geom_hline(yintercept = 0, linetype = "31") +
  labs(x = "Weeks", y = "log(HR)")
ggsave(filename = "fig/fig2.png", plot = p, width = 6, h = 5)
