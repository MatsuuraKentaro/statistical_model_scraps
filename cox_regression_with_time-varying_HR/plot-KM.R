library(survival)
library(ggsurvfit)
library(tidyverse)

d <- read_csv("data/BRAIN.csv")

fit <- survfit2(Surv(weeks, event) ~ resect75, data = d)
p <- ggsurvfit(fit) +
  theme(text = element_text(size = 20)) +
  add_confidence_interval()
ggsave(filename = "fig/fig1.png", plot = p, width = 6, h = 5)
