
## ------------------------------------------------------------
## Analyses for "Social genetic effects on fitness, Strickland et al."; in prep
## Kasha Strickland & Joel Pick, 2026
## ------------------------------------------------------------

## ------------------------------------------------------------
## Load packages and set directory
## ------------------------------------------------------------

library(MCMCglmm)
library(nadiv)
library(corpcor)
library(rstan)
library(ggplot2)
library(tidybayes)

setwd("/Users/kashastrickland/Documents/GitHub/hyenaindirecteffects")

#load helper functions
source("./R/functions.R")

##run data prep code
source("./R/prepare_standata.R")

## ------------------------------------------------------------
## Run models
## ------------------------------------------------------------

## ------------------------------------------------------------
## VA only model
## ------------------------------------------------------------

stan_model_NB_Va <- stan_model(file = "R/stan/VA_hyenas.stan")
##run model sampling
mod_stan_NB_Va <- sampling(stan_model_NB_Va, data = stan_dat_IGE_NB_HRO, cores=3,
                           chains=3, 
                           iter = 30000, warmup = 2000, 
                           pars=c("b","shape","var_mother","var_by","var_clan","var_animal",
                                  "sd_mother","sd_by","sd_clan","sd_animal"))

summary(mod_stan_NB_Va)$summary[,c(1,4,8:10)]
pairs(mod_stan_NB_Va)
rstan::traceplot(mod_stan_NB_Va,
                 pars=c("b","shape","var_mother","var_by","var_clan","var_animal",
                        "sd_mother","sd_by","sd_clan","sd_animal"))

## ------------------------------------------------------------
## VA + HRO as covariance matrix###
## ------------------------------------------------------------

stan_model_NB_Vahro <- stan_model(file = "R/stan/VA_HRO_hyenas.stan")

##run model sampling
mod_stan_NB_Vahro<- sampling(stan_model_NB_Vahro, data = stan_dat_IGE_NB_HRO, cores=3,
                                chains=3, 
                                iter = 10000, warmup = 2000,
                                pars=c("b","shape","var_mother","var_by","var_clan","var_animal","var_hro",
                                       "sd_mother","sd_by","sd_clan","sd_animal","sd_hro"))

summary(mod_stan_NB_Vahro)$summary[,c(1,4,8:10)]
pairs(mod_stan_NB_Vahro)
rstan::traceplot(mod_stan_NB_Vahro,
                 pars=c("b","shape","var_mother","var_by","var_clan","var_animal","var_hro"))

## ------------------------------------------------------------
## VA + soc similarity + HRO as covariance matrix
## ------------------------------------------------------------

stan_model_NB_Vahrosocsim <- stan_model(file = "R/stan/VA_HROSocSim_hyenas.stan")

##run model sampling
mod_stan_NB_Vahrosocsim<- sampling(stan_model_NB_Vahrosocsim, data = stan_dat_IGE_NB_HRO, cores=3,
                             chains=3, 
                             iter = 10000, warmup = 2000, 
                             pars=c("b","shape","var_mother","var_by","var_clan","var_animal","var_hro","var_ss",
                                    "sd_mother","sd_by","sd_clan","sd_animal","sd_hro","sd_ss"))

summary(mod_stan_NB_Vahrosocsim)$summary[,c(1,4,8:10)]
pairs(mod_stan_NB_Vahrosocsim)
rstan::traceplot(mod_stan_NB_Vahrosocsim,
                 pars=c("b","shape","var_mother","var_by","var_clan","var_animal","var_hro","var_ss"))

## ------------------------------------------------------------
## set network as design matrix
## ------------------------------------------------------------

## ------------------------------------------------------------
## VA + Indirect effects (network as design matrix) + HRO
## ------------------------------------------------------------

stan_model_NB_VaIEHRO <- stan_model(file = "R/stan/VA_IE_hyenas.stan")

##sampling
mod_stan_NB_VaIEHRO <- sampling(stan_model_NB_VaIEHRO, data = stan_dat_IGE_NB_HRO, cores=3,
                                chains=3, 
                                iter = 10000, warmup = 2000, 
                                pars=c("b","shape","var_mother","var_by","var_clan","var_hro","var_indirect","var_animal",
                                       "sd_mother","sd_by","sd_clan","sd_hro","sd_indirect","sd_animal"))

summary(mod_stan_NB_VaIEHRO)$summary[,c(1,4,8:10)]
pairs(mod_stan_NB_VaIEHRO)
rstan::traceplot(mod_stan_NB_VaIEHRO,
                 pars=c("b","shape","var_mother","var_by","var_clan","var_hro","var_animal"))

## ------------------------------------------------------------
## VA + Indirect effects (G & E) + HRO
## ------------------------------------------------------------

stan_model_IGE_NB_VA_fullIGE <- stan_model(file = "R/stan/VA_IGE_IEE_hyenas.stan")

## va and hro and ige
mod_stan_IGE_NB_VA_fullIGE <- sampling(stan_model_IGE_NB_VA_fullIGE, data = stan_dat_IGE_NB_HRO, cores=3,
                            chains=3, 
                            iter = 10000, warmup = 2000, 
                            pars=c("b","shape","var_mother","var_by","var_clan","var_indirectE","var_hro",
                                   "sd_by","sd_clan","sd_indirectE","sd_hro",
                                   "Sigma_A","corr_A"))

summary(mod_stan_IGE_NB_VA_fullIGE)$summary[,c(1,4,8:10)]
pairs(mod_stan_IGE_NB_VA_fullIGE)
rstan::traceplot(mod_stan_IGE_NB_VA_fullIGE,
                 pars=c("b","shape","var_mother","var_by","var_clan","var_indirectE","var_hro",
                        "Sigma_A","corr_A"))

## ------------------------------------------------------------
## Save all model outputs
## ------------------------------------------------------------
save.image(file="./models/full_IGE_models_SRI.RData")


## ------------------------------------------------------------
## Plot variances
## ------------------------------------------------------------

##extract everything
models_list <- ls()
models_list <- models_list[sapply(models_list, function(x) inherits(get(x), "stanfit"))]
models_list <- setNames(lapply(models_list, get), models_list)

var_dfs<-list()

for(i in 1:length(models_list)){
md_posterior <- rstan::extract(models_list[[i]])
md_vars <- as.data.frame(models_list[[i]], 
              pars = names(md_posterior)[(grep("var",names(md_posterior)))])
if(grepl("IGE",names(models_list[i]))){
##G matrix
  IGE_Gmat <- as.data.frame(models_list[[i]],pars = c("Sigma_A"))
  md_vars$var_animal<-IGE_Gmat$`Sigma_A[1,1]`
  md_vars$var_IGE<-IGE_Gmat$`Sigma_A[2,2]`
}
else{
}
md_vars_long<-reshape2::melt(md_vars)
md_vars_long$model<-names(models_list[i])
var_dfs[[i]] <- md_vars_long
}

allvars<-do.call("rbind",var_dfs)
allvars$sd<-sqrt(allvars$value)

subset_mods<-allvars[allvars$model==c("mod_stan_NB_Va","mod_stan_NB_Vahrosocsim",
                                      "mod_stan_NB_VaIEHRO","mod_stan_IGE_NB_VA_fullIGE"),]

subset_mods$variable<-factor(subset_mods$variable,levels=c("var_IGE","var_indirectE", "var_indirect","var_ss",
                                                           "var_hro","var_by","var_clan","var_mother","var_animal"))

subset_mods$model<-factor(subset_mods$model,levels=c("mod_stan_IGE_NB_VA_fullIGE","mod_stan_NB_VaIEHRO",
                                                        "mod_stan_NB_Vahrosocsim","mod_stan_NB_Va"))

####plot everything
p_colors<-RColorBrewer::brewer.pal(8,"Dark2")[3:8]

pdf("./tables_figures/Fig2.pdf",width=11,height=9)
ggplot(subset_mods,aes(y=variable,x = sd,colour =model)) +
  geom_vline(xintercept=0,linetype=2,alpha=1/3,col="black",key_glyph = "path")+
  stat_pointinterval(position = position_dodge(0.7))+
  theme_classic()+
  #facet_grid(~model)+
  scale_color_manual(name="Model",breaks = c("mod_stan_NB_Va","mod_stan_NB_Vahrosocsim",
                                             "mod_stan_NB_VaIEHRO","mod_stan_IGE_NB_VA_fullIGE"),
                   labels=c("Additive genetic","Social similarity","IE","Full IGE"),values=p_colors)+
  scale_y_discrete(labels=c("Indirect genetic","Indirect environment","Indirect effects (total)","Social covariance","Shared environment",
                            "Cohort","Clan","Maternal","Additive genetic"))+
  ylab("")+
  xlab(expression(V^~~~LRS))+
  theme(text=element_text(size=16),
        legend.justification = c(1, 1), legend.position = c(0.99, 0.99),
        legend.title.align=0.5)
dev.off()



