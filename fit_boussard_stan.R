## Use stan to fit the RM model to the Boussard data


library("shinystan") 
library(data.table)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
here()

## Read the data 

boussard_data <- read_excel(here("data","doi_10.5061_dryad.5mkkwh72s__v2",
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
boussard_data[is.na(success),success:=0]

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

