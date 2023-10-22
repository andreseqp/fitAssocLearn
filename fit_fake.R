## Use stan to fit the RM model to the Boussard data

library("shinystan") 
library(data.table)
library(readxl)
library(here)
library(tidyverse)
library(Rcpp)
here()

## Read the data 

boussard_data <- read_excel(here("data","doi_10.5061_dryad.5mkkwh72s__v2",
                                 "SRL(all_data).xlsx"))


boussard_data <-boussard_data %>% as.data.table()

# Set parameters for the simulated data

pars.gen<-list(tau=1,mu_alpha=0.2,
               alphasT=c(-1,1),sigma_a=0.5)

# set number of individuals
Nind <- 100

# set number of treatment groups
Ntreat <- 2

# set number of trials per reversal
NtriRev <- 200

Nrev <-1

# get total number of trials (including all reversals)
Ntrials <- NtriRev*Nrev

pars.gen$alphasID<-rnorm(Nind,pars.gen$mu_alpha,sd = pars.gen$sigma_a)

pars.gen

# set the treatment for all individuals
treat_Inds <- rep(x=c(0,1),each=Nind/2)

# set the reversal structure 
#block_r <- cbind(rep(0,Ntrials),rep(1,Ntrials))

block_r <-cbind(c(rep(c(0,1),each=NtriRev,
                      times=6)[1:Ntrials]),
                c(rep(c(1,0),each=NtriRev,
                      times=6)[1:Ntrials]))

# Simulate learning for all individuals
prediction.ind <- data.frame(
  val.1 = rep(0,dim(block_r)[1]),
  val.2 = rep(0,dim(block_r)[1]),
  choice = rep(0,dim(block_r)[1]),
  rew.1 = block_r[,1],
  rew.2 = block_r[,2],
  success = rep(0,dim(block_r)[1])
)

# using cmdstanr
library(cmdstanr)
# Set cmdstan path 
## Erase to run in cluster
stan_path <- here("..","cmdstan","cmdstan-2.32.2")
set_cmdstan_path(stan_path)

# Compile stan model
boussard_RW_cmd<-cmdstan_model("boussard_RW.stan")

zeros<-matrix(data=rep(0,Nind*Ntrials),nrow = Nind,ncol = Ntrials)

# Simulate data using stan
sim_data_stan <- boussard_RW_cmd$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                             block_r=block_r,
                                             treat_ID=treat_Inds,
                                             y=zeros),
                                        fixed_param = TRUE,chains = 1,
                                        iter_sampling = 1,
                                        init=list(pars.gen))



summSim<-sim_data_stan$summary()

predictions <- summSim %>% filter(grepl("y_pred",variable)) %>% 
  select(c("variable","mean")) %>% 
  separate(variable,into=c("Individual","Trial"),sep = ",") %>% 
  mutate(Individual=parse_number(Individual),Trial=parse_number(Trial))

predictions <-rename(predictions,success=mean)

pred_wide<-pivot_wider(predictions,names_from = Trial,
                       values_from = "success") %>%
  mutate(Individual=NULL)
  

# sample from posterior
fit_simulated_RW_cmd <- boussard_RW_cmd$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                           block_r=block_r,
                           treat_ID=treat_Inds,
                           y=as.matrix(pred_wide)),
                           parallel_chains = getOption("mc.cores", 5),
                           chains = 5,iter_sampling = 1000)

# Save samples to file
fit_simulated_RW_cmd$save_object(file = "fit_sim_stan_trials.RDS")


# Save samples to file
fit_boussard_RW_cmd$save_object(file = "fit_boussard_stan.RDS")

fit_simulated_RW_cmd<-readRDS("fit_sim_stan_trials.RDS")

# Use shinystan to evaluate the performance of the model
launch_shinystan(fit_simulated_RW_cmd)
