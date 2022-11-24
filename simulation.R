# packages
library(tidyverse)
library(broom)

# set seed
set.seed(1116)

# define a expit (logistic sigmoid) function
expit<-function(x){1/(1+exp(-(x)))}

n<-100000 # number of patients
nsim = 1000 # number of simulation
g = c(0.001, 0.01, 0.1, 0, 1, 2) # gene
e = c(0.001, 0.01, 0.1, 0, 1, 2) # environment
interaction = c(0) # interaction term
results = expand.grid(g = g,
                      e = e,
                      interaction = interaction)
# creating a simulation dataset

for(simInd in 1:nrow(results)){
  g = results[simInd, "g"]
  e = results[simInd, "e"]
  interaction = results[simInd, "interaction"]
  eval <- c()
  
  for (i in 1:nsim) {
    set.seed(i)
    simdata = tibble::tibble(
      gene = rbinom(n, 1, .25),
      environment = rbinom(n, 1, .25),
      outcome = rbinom(n,
                       1,
                       expit(g*gene+e*environment+interaction*gene*environment)))
    fit <- glm(outcome ~ gene*environment, 
               family = binomial(link = "logit"),
               data = simdata)
    sum <- tidy(fit, exponentiate = FALSE, conf.int = TRUE)
    l.ci <- sum[4, "conf.low"]
    u.ci <- sum[4, "conf.high"]
    eval <- c(eval, if_else(l.ci < interaction & interaction < u.ci, 1, 0))
  }
  results[simInd, "eval"] = sum(eval)/nsim
}
