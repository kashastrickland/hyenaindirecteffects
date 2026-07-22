
## simulate from the posterior of the full IGE model to check statistical power to detect paramters
## should be run after "prepare_standata.R" and "IGE_models_SRI.R"

library(future)
library(future.apply)
library(MCMCglmm)
library(nadiv)
library(corpcor)
library(rstan)
library(ggplot2)
library(tidybayes)
library(posterior)

setwd("~/Documents/GitHub/hyenaindirecteffects")

rstan_options(auto_write = TRUE)
set.seed(12340)

load("./models/full_IGE_models_SRI.RData")
source("./R/sims_func.R")

plan(multisession, workers = 4)
dir.create("sim_fits_full", showWarnings = FALSE)

##draw from posterior
# Extract all post-warmup draws as a dataframe
draws_df <- as.data.frame(mod_stan_IGE_NB_VA_fullIGE)
colnames(draws_df)[1:6]<-colnames(X)

# Extract 100 random draw
random_draw <- draws_df[sample(nrow(draws_df), 100,replace = FALSE), ]

###
### Simulation ####
###

IGE_sim_files <- future_lapply(
  seq_len(nrow(random_draw)),
  run_one_sim_and_fit,
  future.seed = TRUE,
  future.packages = c("rstan", "MCMCglmm")
)
