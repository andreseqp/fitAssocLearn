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

# Set parameters for the simulated data

pars<-

# set number of individuals
Nind <- 20

# get number of treatment groups
Ntreat <- 2

# get total number of trials (including all reversals)
Ntrials <- 50

setorder(boussard_data,tankID,reversal,trial)

# set the treatment for all individuals
treat_Inds <- rep(x=c(0,1),each=Nind/2)
  # boussard_data[,unique(brainsize),by=tankID][,V1]  

# set the reversal structure 
block_r <- cbind(rep(0,Ntrials),rep(1,Ntrials))
  

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

# Compile stan model
boussard_RW_cmd<-cmdstan_model("boussard_RW.stan")

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
