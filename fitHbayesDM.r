## Let´s see how well our model can be fit to HBayesDM datasets
library(here)
library(tidyverse)
library("hBayesDM")
library('data.table')
library('cmdstanr')
# Set cmdstan path 
## Erase to run in cluster
stan_path <- here("..","cmdstan","cmdstan-2.32.2")
set_cmdstan_path(stan_path)

# Let´s choose one of the data sets
prl_example <- fread(here('data',"prl_exampleData.txt")) %>% as.data.table()
str(prl_example)

prl_ewa <- cmdstan_model("prl_ewa.stan")

Nind<-prl_example[,unique(subjID)] %>% length()
Trials<-prl_example[,max(trial)] 
TrialSub<-prl_example[,max(trial),by=subjID]$V1
choice<-dcast(prl_example,trial~subjID,
              value.var = "choice")
choice[,trial:=NULL]
choice<-t(choice)

reward<-dcast(prl_example,trial~subjID,
              value.var = "outcome")
reward[,trial:=NULL]
reward<-t(reward)

fit_prl_ewa <- prl_ewa$sample(data=list(
  N=Nind,T=Trials,Tsubj=TrialSub,choice=choice,outcome=reward),
  iter_sampling = 2000,parallel_chains =  4,
  iter_warmup=1000, chains=4
)

# Save samples to file
fit_prl_ewa$save_object(file = "fit_prl_ewa.RDS")
# read MCMC files
fit_prl_ewa<-readRDS("fit_prl_ewa.RDS")
# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_prl_ewa)

## Use boussard inspired model to fit hBayesDM dataset

# get number of individuals
Nind <- prl_example[,subjID] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- prl_example[,trial] %>% 
  unique() %>% length()

# set the reversal structure 

# library(abind)
# block_r <- data.matrix(prl_example[,.(rew_1=ifelse(choice==1,
#                                   outcome,-outcome),
#                              rew_2=ifelse(choice==2,
#                                   outcome,-outcome))]
# )
# library(abind)
# 
# block_r_3d <-  abind(lapply(split(seq_len(nrow(block_r)), 
#                             (seq_len(nrow(block_r))-1) %/% Ntrials + 1), 
#                       function(x) block_r[x, ]), along = 3)



# Transform the choice variable, to fit stan model
prl_choice_wide<-dcast(prl_example[,.(trial,subjID,choice)],
                     trial~subjID,value.var = "choice")
prl_choice_wide[,trial:=NULL]
prl_choice_wide<-t(prl_choice_wide)
str(prl_choice_wide)

# Transform the choice variable, to fit stan model
prl_reward_wide<-dcast(prl_example[,.(trial,subjID,outcome)],
                       trial~subjID,value.var = "outcome")
prl_reward_wide[,trial:=NULL]
prl_reward_wide<-t(prl_reward_wide)
str(prl_reward_wide)


# Compile stan model
prl_RW<-cmdstan_model("prl_RW.stan")

# sample from posterior
fit_prl_Data <- prl_RW$sample(list(N=Nind,
                             Tr=Ntrials,
                             reward=prl_reward_wide,
                             choice=prl_choice_wide),
                        parallel_chains = getOption("mc.cores", 5),
                        chains = 5)


# Save samples to file
fit_prl_Data$save_object(file = "fit_prl_Data.RDS")

fit_prl_Data<-readRDS("fit_prl_Data.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_prl_Data)



