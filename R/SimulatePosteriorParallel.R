
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

###
#### ------ Process simulations
###

IGE_sim_files<-list.files("sim_fits_full")

pars<-c("(Intercept)","sexmale","mother_rank_ord","clan_size","f","X","var_mother",
        "sd_by","sd_clan","sd_indirectE","sd_hro","Sigma_A[1,1]","Sigma_A[2,1]","Sigma_A[1,2]","Sigma_A[2,2]")

CI_sims<-list()
for (i in 1:length(IGE_sim_files)){
  sim_list<-readRDS(IGE_sim_files[[i]])
  sim_fit<-sim_list$fit
  names(sim_fit)[1:6]<-pars[1:6]
  CI <- as.data.frame(t(sapply(pars, function(x) {
    quantile(rstan::extract(sim_fit, pars = x)[[1]],probs = c(0.05,0.95),na.rm = TRUE)
  })
  ))
  
  CI$parameter <- rownames(CI)
  rownames(CI) <- NULL
  CI$true <- as.numeric(sim_list$true_parameters[CI$parameter])
  CI <- CI[match(pars, CI$parameter), ]
  CI$ci_95<-ifelse(CI$true>=CI$`5%`&CI$true<=CI$`95%`,TRUE,FALSE)
  CI_sims[[i]]<-CI[c("parameter","ci_95")]
  cat(i)
}

CI_sims_df<-do.call("rbind",CI_sims)
CI_sims_df_sp<-split(CI_sims_df,CI_sims_df$parameter)

CI95<-lapply(CI_sims_df_sp, function(x) sum(x$ci_95 == TRUE, na.rm = TRUE))




