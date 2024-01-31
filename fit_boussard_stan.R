## Use stan to fit the RM model to the Boussard data

## Libraries --------------------------------------------------------
library("shinystan") 
library(data.table)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(tidyverse)
theme_set(new = theme_classic())
here()

## Read the data 

boussard_data <- read_excel(here("data","boussard_2020",
                                 "SRL(all_data).xlsx"))


boussard_data <-boussard_data %>% as.data.table()

# boussard_data <-boussard_data[tankID<5 &
#                                 reversal<2]

# get number of individuals
Nind <- boussard_data[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data[,interaction(trial,reversal)] %>% 
              unique() %>% length()

setorder(boussard_data,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data[is.na(success),success:=-1]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

# using cmdstanr

library(cmdstanr)
# Set cmdstan path 
## Erase to run in cluster
stan_path <- here("..","cmdstan","cmdstan-2.32.2")
set_cmdstan_path(stan_path)

# Compile stan model with random alpha
boussard_RW_cmd<-cmdstan_model("boussard_RW.stan")

dim(block_r)

# sample from posterior
fit_boussard_RW_cmd <- boussard_RW_cmd$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                 block_r=block_r,
                                 treat_ID=treat_Inds,
                                 y=boussard_wide),
                                 parallel_chains = getOption("mc.cores", 5),
                                 chains = 5)

# Save samples to file
fit_boussard_RW_cmd$save_object(file = "fit_boussard_stan.RDS")

fit_boussard_RW_cmd<-readRDS("fit_boussard_stan.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_boussard_RW_cmd)


# Compile stan model with random alpha and tau
boussard_RW_2<-cmdstan_model("boussard_RW_2.stan")

dim(block_r)

# sample from posterior
fit_boussard_RW_2 <- boussard_RW_2$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                                   block_r=block_r,
                                                   treat_ID=treat_Inds,
                                                   y=boussard_wide),
                                              parallel_chains = getOption("mc.cores", 5),
                                              chains = 5)

# Save samples to file
fit_boussard_RW_2$save_object(file = "fit_boussard_stan_2.RDS")

fit_boussard_RW_2<-readRDS("fit_boussard_stan_2.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_boussard_RW_2)


# fit_boussard_RW_cmd<-readRDS("MCMC_alphas.rda")

# print(fit_boussard_RW_cmd)
# launch_shinystan(fit_boussard_RW_cmd)

# fit_boussard_RW_cmd <- readRDS(temp_rds_file)


# # using rstan
# library("rstan")
# 
# boussard_RW<-stan_model("boussard_RW.stan")
# 
# fit_boussard_RW <-sampling(boussard_RW,list(N=Nind,B=Ntreat,Tr=Ntrials,
#                                            block_r=block_r,
#                                            treat_ID=treat_Inds,
#                                            y=boussard_wide))

## Shortening the data-set -----------------------------------------------------
## Let's try to find out why the model with flexible learning rates 
# does not fit the Boussard data set so well.

boussard_data_short <- boussard_data[reversal==0,]

# get number of individuals
Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>% 
  unique() %>% length()

setorder(boussard_data_short,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data_short[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data_short[is.na(success),success:=0]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

boussard_RW_cmd<-cmdstan_model(here("stanModels","boussard_RW.stan"))

dim(block_r)

# sample from posterior
fit_boussard_RW_short <- boussard_RW_cmd$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                                   block_r=block_r,
                                                   treat_ID=treat_Inds,
                                                   y=boussard_wide),
                                              parallel_chains = getOption("mc.cores", 5),
                                              chains = 5)

# Save samples to file
# fit_boussard_RW_short$save_object(file = "fit_boussard_stan_short.RDS")
# 
# fit_boussard_RW_short<-readRDS("fit_boussard_stan_short.RDS")

# Use shinystan to evaluate the performance of the model
# launch_shinystan(fit_boussard_RW_short)

pars2plot <- c("tau", "mu_alpha", "alphasT[1]", "alphasT[2]",
                "sigma_a")

posteriors<-fit_boussard_RW_short$draws(
  variables = pars2plot)


# png("boussard_mcmc_short.png")
mcmc_intervals(posteriors)
# dev.off()

preds <- fit_boussard_RW_short$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                   measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]


# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]

rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= as.integer((trial-1)/30)) %>%
  mutate(RTrial=((trial-1) %% 30)+1)

# mean_ind[,reversal:=reversal+1]

colnames(mean_ind) <- c("chain","iteration","brainsize","Ttrial","success",
                        "reversal","trial")


# png("boussard_ppchecks.png")
boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=brainsize))+
    stat_summary(fun = mean,geom = "point")+
    stat_summary(fun = mean,geom = "line")+
    # geom_point()+
    theme(legend.position = c(0.5,0.6),
          legend.direction = "horizontal",
          strip.text.y = element_blank())+
    scale_x_continuous(breaks=c(1,15,30))+
    guides(fill=guide_legend(title="Brain size"))+
    ggtitle("Repeated reversal vs brainsize")+
    facet_grid(brainsize~reversal)+
    stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.2,fun.max = function(x){
      quantile(x,0.95)},
      fun.min = function(x){
      quantile(x,0.05)},colour=NA)+
    stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.5,fun.max = function(x){
      quantile(x,0.75)},
      fun.min = function(x){
      quantile(x,0.25)},colour=NA)+
    facet_grid(brainsize~reversal)
# dev.off()


## Shortened data-set with flexible learning rates -----------------------------

boussard_data_short <- boussard_data[reversal==0,]

# get number of individuals
Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>%
  unique() %>% length()

Nrev <- boussard_data_short[,reversal] %>%
  unique() %>% length()

NtriRev <- Ntrials/Nrev

setorder(boussard_data_short,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]

# set the reversal structure
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data_short[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data_short[is.na(success),success:=0]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)


# Compile stan model with random alpha
boussard_RW_rev<-cmdstan_model(here("stanModels","boussard_RW_rev.stan"))


# sample from posterior
fit_boussard_RW_rev_short <- boussard_RW_rev$sample(list(N=Nind,B=Ntreat,Tr=NtriRev,
                                                   Rev=Nrev,TotTr=Ntrials,
                                                   block_r=block_r,treat_ID=treat_Inds,
                                                   y=boussard_wide),
                                              parallel_chains = getOption("mc.cores", 5),
                                              chains = 5)

# Save samples to file
fit_boussard_RW_rev_short$save_object(file = "fit_boussard_rev_short.RDS")

fit_boussard_RW_rev_short<-readRDS("fit_boussard_rev_short.RDS")


pars2plot <- c("tau", "mu_alpha", "alphasT[1]", "alphasT[2]",
                "sigma_a")

posteriors<-fit_boussard_RW_rev_short$draws(
  variables = pars2plot)

posteriorsRev <- fit_boussard_RW_rev_short$draws(variables = "alphasRev")

# png("boussard_mcmc_interv_rev.png")
cowplot::plot_grid(mcmc_intervals(posteriors),
                   mcmc_intervals(posteriorsRev))
# dev.off()

preds <- fit_boussard_RW_short$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                             measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]


# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]

rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= as.integer((trial-1)/30)) %>%
  mutate(RTrial=((trial-1) %% 30)+1)

colnames(mean_ind) <- c("chain","iteration","brainsize","Ttrial","success",
                        "reversal","trial")


# png("boussard_ppchecks.png")
boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=brainsize))+
  stat_summary(fun = mean,geom = "point")+
  stat_summary(fun = mean,geom = "line")+
  # geom_point()+
  theme(legend.position = c(0.3,0.6),
        legend.direction = "horizontal",
        strip.text.y = element_blank())+
  scale_x_continuous(breaks=c(1,15,30))+
  guides(fill=guide_legend(title="Brain size"))+
  ggtitle("Repeated reversal vs brainsize")+
  facet_grid(brainsize~reversal)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.2,fun.max = function(x){
                 quantile(x,0.95)},
               fun.min = function(x){
                 quantile(x,0.05)},colour=NA)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.5,fun.max = function(x){
                 quantile(x,0.75)},
               fun.min = function(x){
                 quantile(x,0.25)},colour=NA)+
  facet_grid(brainsize~reversal)

## shortened data-set with tau as a fixed effect ------------------------------

boussard_data_short <- boussard_data[reversal<1,]

# get number of individuals
Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>% 
  unique() %>% length()

setorder(boussard_data_short,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data_short[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data_short[is.na(success),success:=0]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

boussard_RW_tau<-cmdstan_model(here("stanModels","boussard_RW_tau.stan"))

dim(block_r)

# sample from posterior
fit_boussard_RW_short_tau <- boussard_RW_tau$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                                     block_r=block_r,
                                                     treat_ID=treat_Inds,
                                                     y=boussard_wide),
                                                parallel_chains = getOption("mc.cores", 5),
                                                chains = 5)

# Save samples to file
fit_boussard_RW_short_tau$save_object(file = "fit_boussard_stan_short_tau.RDS")

fit_boussard_RW_short_tau<-readRDS("fit_boussard_stan_short_tau.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_boussard_RW_short_tau)

pars2plot <- c("alpha", "mu_tau", "tausT[1]", "tausT[2]",
               "sigma_t")

posteriors<-fit_boussard_RW_short_tau$draws(
  variables = pars2plot)


# png("boussard_mcmc_short.png")
mcmc_intervals(posteriors)
# dev.off()

preds <- fit_boussard_RW_short_tau$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                             measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]


# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]

rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= as.integer((trial-1)/30)) %>%
  mutate(RTrial=((trial-1) %% 30)+1)

colnames(mean_ind) <- c("chain","iteration","brainsize","Ttrial","success",
                        "reversal","trial")


# png("boussard_ppchecks.png")
boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=brainsize))+
  stat_summary(fun = mean,geom = "point")+
  stat_summary(fun = mean,geom = "line")+
  # geom_point()+
  theme(legend.position = c(0.5,0.6),
        legend.direction = "horizontal",
        strip.text.y = element_blank())+
  scale_x_continuous(breaks=c(1,15,30))+
  guides(fill=guide_legend(title="Brain size"))+
  ggtitle("Repeated reversal vs brainsize")+
  facet_grid(brainsize~reversal)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.2,fun.max = function(x){
                 quantile(x,0.95)},
               fun.min = function(x){
                 quantile(x,0.05)},colour=NA)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.5,fun.max = function(x){
                 quantile(x,0.75)},
               fun.min = function(x){
                 quantile(x,0.25)},colour=NA)+
  facet_grid(brainsize~reversal)
# dev.off()



## Let´s fit a different tau for each reversal block ---------------------------

boussard_data_short <- boussard_data[reversal<2,]

# get number of individuals
Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>%
  unique() %>% length()

Nrev <- boussard_data_short[,reversal] %>%
  unique() %>% length()

NtriRev <- Ntrials/Nrev

setorder(boussard_data_short,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data_short[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data_short[is.na(success),success:=0]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

boussard_RW_tau_rev<-cmdstan_model(here("stanModels","boussard_RW_tau_rev.stan"))

dim(block_r)

# sample from posterior
fit_boussard_RW_short_tau_rev <- boussard_RW_tau_rev$sample(list(N=Nind,B=Ntreat,Tr=NtriRev,
                                                         Rev=Nrev,TotTr=Ntrials,
                                                         block_r=block_r,treat_ID=treat_Inds,
                                                         y=boussard_wide),
                                                    parallel_chains = getOption("mc.cores", 5),
                                                    chains = 5)

# Save samples to file
fit_boussard_RW_short_tau_rev$save_object(file = "fit_boussard_stan_short_tau_rev.RDS")

fit_boussard_RW_short_tau_rev<-readRDS("fit_boussard_stan_short_tau_rev.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_boussard_RW_short_tau_rev)

pars2plot <- c("alpha", "mu_tau", "tausT[1]", "tausT[2]",
               "sigma_t")

posteriors <- fit_boussard_RW_short_tau_rev$draws(variables = pars2plot)
posteriorsRev <- fit_boussard_RW_short_tau_rev$draws(variables = "tausRev")

# png("boussard_mcmc_interv_rev.png")
cowplot::plot_grid(mcmc_intervals(posteriors),
                   mcmc_intervals(posteriorsRev))
# dev.off()

preds <- fit_boussard_RW_short_tau_rev$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                             measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]


# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]

rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= as.integer((trial-1)/30)) %>%
  mutate(RTrial=((trial-1) %% 30)+1)

colnames(mean_ind) <- c("chain","iteration","brainsize","Ttrial","success",
                        "reversal","trial")


# png("boussard_ppchecks.png")
boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=brainsize))+
  stat_summary(fun = mean,geom = "point")+
  stat_summary(fun = mean,geom = "line")+
  # geom_point()+
  theme(legend.position = c(0.5,0.6),
        legend.direction = "horizontal")+
  scale_x_continuous(breaks=c(1,15,30))+
  guides(fill=guide_legend(title="Brain size"))+
  ggtitle("Repeated reversal vs brainsize")+
  facet_grid(brainsize~reversal)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.2,fun.max = function(x){
                 quantile(x,0.95)},
               fun.min = function(x){
                 quantile(x,0.05)},colour=NA)+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.5,fun.max = function(x){
                 quantile(x,0.75)},
               fun.min = function(x){
                 quantile(x,0.25)},colour=NA)+
  facet_grid(brainsize~reversal)
# dev.off()


## Let's fit a different tau and alpha for each reversal block -----------------


boussard_data_short <- boussard_data[reversal<2,]

# get number of individuals
Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>%
  unique() %>% length()

Nrev <- boussard_data_short[,reversal] %>%
  unique() %>% length()

NtriRev <- Ntrials/Nrev

setorder(boussard_data_short,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data_short[tankID==1])[1]]))


# Get a unique ID for trials along the reversal blocks
boussard_data_short[,trial.long:=interaction(reversal,trial)]

# Set NA as failure to choose the rewarding option
boussard_data_short[is.na(success),success:=0]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

boussard_RW_tau_alpha_rev<-cmdstan_model(here("stanModels","boussard_RW_tau_alpha_rev.stan"))

dim(block_r)

# sample from posterior
fit_boussard_RW_short_tau_alpha_rev <- boussard_RW_tau_alpha_rev$sample(
  list(N=Nind,B=Ntreat,Tr=NtriRev,
           Rev=Nrev,TotTr=Ntrials,
           block_r=block_r,treat_ID=treat_Inds,
           y=boussard_wide),
      parallel_chains = getOption("mc.cores", 5),
      chains = 5)

# Save samples to file
# fit_boussard_RW_short_tau_alpha_rev$save_object(file = 
  # "fit_boussard_stan_short_tau_alpha_rev.RDS")

# fit_boussard_RW_short_tau_alpha_rev<-readRDS("fit_boussard_stan_short_tau_alpha_rev.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_boussard_RW_short_tau_alpha_rev)

pars2plot <- c("tausT[1]", "tausT[2]")

posteriors <- fit_boussard_RW_short_tau_alpha_rev$draws(variables = pars2plot)
posteriorsTRev <- fit_boussard_RW_short_tau_alpha_rev$draws(variables = "tausRev")
posteriorsARev <- fit_boussard_RW_short_tau_alpha_rev$draws(variables = "alphasRev")
posteriorsSigRev <- 
  fit_boussard_RW_short_tau_alpha_rev$draws(variables = c("alphaSigmasRev",
                                                          "tauSigmasRev"))

# png("boussard_mcmc_interv_rev.png")
interv_01<-cowplot::plot_grid(mcmc_intervals(posteriors),
                   mcmc_intervals(posteriorsTRev),
                   mcmc_intervals(posteriorsARev),
                   mcmc_intervals(posteriorsSigRev))
# dev.off()

cowplot::plot_grid(interv_0,interv_1,interv_01)

preds <- fit_boussard_RW_short_tau_alpha_rev$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                             measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]


# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]

rm(list="preds_long")


mean_ind_01 <-mean_ind %>%
  mutate(reversal= as.integer((trial-1)/30)) %>%
  mutate(RTrial=((trial-1) %% 30)+1)

# mean_ind_1[,reversal:=reversal+1]

colnames(mean_ind_01) <- c("chain","iteration","brainsize","Ttrial","success",
                        "reversal","trial")


# png("boussard_ppchecks.png")
cowplot::plot_grid(nrow = 2,
  boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
    ggplot(aes(y=success,x=trial,col=brainsize))+
    stat_summary(fun = mean,geom = "point")+
    stat_summary(fun = mean,geom = "line")+
    geom_hline(yintercept = c(0.5,1))+
    theme(legend.position = c(0.3,0.6),
          legend.direction = "horizontal")+
    scale_x_continuous(breaks=c(1,15,30))+
    guides(fill=guide_legend(title="Brain size"))+
    ggtitle("Repeated reversal vs brainsize")+
    facet_grid(brainsize~reversal)+
    stat_summary(data=mean_ind_01,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.2,fun.max = function(x){
                   quantile(x,0.95)},
                 fun.min = function(x){
                   quantile(x,0.05)},colour=NA)+
    stat_summary(data=mean_ind_01,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.5,fun.max = function(x){
                   quantile(x,0.75)},
                 fun.min = function(x){
                   quantile(x,0.25)},colour=NA)+
    facet_grid(brainsize~reversal),
    cowplot::plot_grid(boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
    filter(reversal==0) %>% 
    ggplot(aes(y=success,x=trial,col=brainsize))+
    stat_summary(fun = mean,geom = "point")+
    stat_summary(fun = mean,geom = "line")+
    geom_hline(yintercept = c(0.5,1))+
    theme(legend.position = c(0.3,0.6),
            legend.direction = "horizontal")+
      scale_x_continuous(breaks=c(1,15,30))+
    stat_summary(data=mean_ind_0,aes(x=trial,y=success,col=as.factor(brainsize)),
                        geom="ribbon",alpha = 0.2,fun.max = function(x){
                          quantile(x,0.95)},
                        fun.min = function(x){
                          quantile(x,0.05)},colour=NA)+
    stat_summary(data=mean_ind_0,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.5,fun.max = function(x){
                   quantile(x,0.75)},
                 fun.min = function(x){
                   quantile(x,0.25)},colour=NA)+
    facet_grid(brainsize~reversal),
    boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
    filter(reversal==1) %>% 
    ggplot(aes(y=success,x=trial,col=brainsize))+
    stat_summary(fun = mean,geom = "point")+
    stat_summary(fun = mean,geom = "line")+
    geom_hline(yintercept = c(0.5,1))+
    scale_x_continuous(breaks=c(1,15,30))+
      theme(legend.position = c(0.3,0.6),
            legend.direction = "horizontal")+
      
    stat_summary(data=mean_ind_1,aes(x=trial,y=success,col=as.factor(brainsize)),
               geom="ribbon",alpha = 0.2,fun.max = function(x){
                 quantile(x,0.95)},
               fun.min = function(x){
                 quantile(x,0.05)},colour=NA)+
    stat_summary(data=mean_ind_1,aes(x=trial,y=success,col=as.factor(brainsize)),
                 geom="ribbon",alpha = 0.5,fun.max = function(x){
                   quantile(x,0.75)},
                 fun.min = function(x){
                   quantile(x,0.25)},colour=NA)+
    facet_grid(brainsize~reversal))
)
# dev.off()

## Fit each individual block separately ---------------------------------------

## Shortened data-set with flexible learning rates 

# get number of individuals
pdf(here('images','allBlockIndi.pdf'))
for(i in 0:10){
  # i <- 1
  print(i)
  boussard_data_short <- boussard_data[reversal==i,]
  
  init_prop <- boussard_data_short[trial==1,mean(success,na.rm=TRUE)]
  
  E_1 <- 0-log(init_prop/(1-init_prop))/2
  
  Nind <- boussard_data_short[,tankID] %>% unique() %>% length()

# get number of treatment groups
  Ntreat <- boussard_data_short[,brainsize] %>% unique() %>% length()
  
  # get total number of trials (including all reversals)
  Ntrials <- boussard_data_short[,interaction(trial,reversal)] %>%
    unique() %>% length()
  
  Nrev <- boussard_data_short[,reversal] %>%
    unique() %>% length()
  
  NtriRev <- Ntrials/Nrev
  
  setorder(boussard_data_short,tankID,reversal,trial)
  
  # get the treatment for all individuals
  treat_Inds <- boussard_data_short[,unique(brainsize),by=tankID][,V1]
  
  # set the reversal structure
  block_r <-cbind(c(rep(c(0,1),each=30,
                        times=6)[1:dim(boussard_data_short[tankID==1])[1]]),
                  c(rep(c(1,0),each=30,
                        times=6)[1:dim(boussard_data_short[tankID==1])[1]]))
  
  
  # Get a unique ID for trials along the reversal blocks
  boussard_data_short[,trial.long:=interaction(reversal,trial)]
  
  # Set NA as failure to choose the rewarding option
  boussard_data_short[is.na(success),success:=0]
  
  # Transform the success variable, to fit stan model
  boussard_wide<-dcast(boussard_data_short[,.(trial.long,tankID,success)],
                       trial.long~tankID,value.var = "success")
  
  boussard_wide[,trial.long:=NULL]
  
  boussard_wide<-t(boussard_wide)
  
  
  # Compile stan model with random alpha
  # boussard_RW_rev<-cmdstan_model(here("stanModels","boussard_RW_rev.stan"))
  
  
  # sample from posterior
  fit_boussard_RW_rev_short <- boussard_RW_rev$sample(list(N=Nind,B=Ntreat,Tr=NtriRev,
                                                     Rev=Nrev,TotTr=Ntrials,
                                                     initProp=init_prop,
                                                     block_r=block_r,treat_ID=treat_Inds,
                                                     y=boussard_wide),
                                                parallel_chains = getOption("mc.cores", 5),
                                                chains = 1,iter_sampling  = 500)
  
  # Save samples to file
  # fit_boussard_RW_rev_short$save_object(file = "fit_boussard_rev_short.RDS")
  
  # fit_boussard_RW_rev_short<-readRDS("fit_boussard_rev_short.RDS")
  
  
  pars2plot <- c("tau", "mu_alpha", "alphasT[1]", "alphasT[2]",
                  "sigma_a")
  
  posteriors<-fit_boussard_RW_rev_short$draws(
    variables = pars2plot)
  
  posteriorsRev <- fit_boussard_RW_rev_short$draws(variables = "alphasRev")
  
  # png("boussard_mcmc_interv_rev.png")
  interv <- cowplot::plot_grid(mcmc_intervals(posteriors),
                     mcmc_intervals(posteriorsRev),nrow = 2)
  # dev.off()
  
  preds <- fit_boussard_RW_rev_short$draws(variables = "y_pred")
  
  preds_df <- posterior::as_draws_df(preds)
  
  preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                               measure.vars=grep("y_pred",colnames(preds_df)))
  
  rm(list=c("preds","preds_df"))
  
  preds_long <- as.data.table(preds_long)
  
  preds_long[,c("individual","trial"):=tstrsplit(variable,",")]
  
  preds_long[,variable:=NULL]
  
  preds_long[,`:=`(individual=parse_number(individual),
                   trial=parse_number(trial))]
  
  preds_long[,treatment:=as.factor(ifelse(individual<50,0,1))]
  
  
  # Average over de MCMC samples
  mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,treatment,trial)]
  
  rm(list="preds_long")
  
  
  mean_ind <-mean_ind %>%
    mutate(reversal= as.integer((trial-1)/30)) %>%
    mutate(RTrial=((trial-1) %% 30)+1)
  
  colnames(mean_ind) <- c("chain","iteration","brainsize","Ttrial","success",
                          "reversal","trial")
  
  mean_ind[,reversal:=reversal+i]
  
  # png(here("images","boussard_short_0_rev.png"))
  print(
    cowplot::plot_grid(ncol = 2,labels = paste0("    Reversal: ",i),
     interv,
      boussard_data_short %>% mutate(brainsize=as.factor(brainsize)) %>%
        ggplot(aes(y=success,x=trial,col=brainsize))+
        stat_summary(fun = mean,geom = "point")+
        stat_summary(fun = mean,geom = "line")+
        # geom_point()+
        theme(legend.position = c(0.6,0.6),
              legend.direction = "vertical",
              strip.text.y = element_blank())+
        scale_x_continuous(breaks=c(1,15,30))+
        guides(fill=guide_legend(title="Brain size \n"))+
        ggtitle("Repeated reversal vs brainsize")+
        facet_grid(brainsize~reversal)+
        stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
                     geom="ribbon",alpha = 0.2,fun.max = function(x){
                       quantile(x,0.95)},
                     fun.min = function(x){
                       quantile(x,0.05)},colour=NA)+
        stat_summary(data=mean_ind,aes(x=trial,y=success,col=as.factor(brainsize)),
                     geom="ribbon",alpha = 0.5,fun.max = function(x){
                       quantile(x,0.75)},
                     fun.min = function(x){
                       quantile(x,0.25)},colour=NA)+
        facet_grid(brainsize~reversal)
    )
  )
}

dev.off()


