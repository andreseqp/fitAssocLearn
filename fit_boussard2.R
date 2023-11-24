## Fitting a RW model to the data of Boussard et al. 2020
##  Boussard, A., S. D. Buechel, M. Amcoff, A. Kotrschal, and N. Kolm. 2020. 
#   Brain size does not predict learning strategies in a serial reversal 
#   learning test. The Journal of Experimental Biology 223:jeb224741.

## Import important packages

library(data.table)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)

## Read the data 

boussard2_data <- read_excel(here("data","boussard2",
                                 "Cognitive_ageing_original_data.xlsx"),
                             sheet = 'FL orig. data',
                             col_types = c('text','numeric','numeric',
                                           'numeric','text','text','text',
                                           'text','numeric','numeric'))
boussard2_data <- boussard2_data %>% mutate(reversal=0)

boussard2_data_b <- read_excel(here("data","boussard2",
                                  "Cognitive_ageing_original_data.xlsx"),
                             sheet = 'RL orig. data',
                             col_types = c('text','numeric','numeric',
                                           'numeric','text','text','text',
                                           'text','numeric','numeric'))

boussard2_data_b <- boussard2_data_b %>% mutate(reversal=1)

boussard2_data[,unique(tankID)] %>% length()

boussard2_data <- rbind(boussard2_data,boussard2_data_b) 

boussard2_data$tankID %>% unique() %>% length()
boussard2_data_b$tankID %>% unique() %>% length()


# tankID coincides with the number of fish supposedly tested in the experimental 
# set up

interaction(boussard2_data$tankID,boussard2_data$fishID) %>% 
  unique() %>% length()
table(boussard2_data[,c("fishID","replicate")])
table(boussard2_data[,c("tankID","fishID")])
table(boussard2_data[,c("tankID","replicate")])
table(boussard2_data[,c("fishID","replicate","tankID")])

# get number of individuals
Nind <- boussard2_data[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- boussard2_data[,brainsize] %>% unique() %>% length()

# get total number of trials (including all reversals)
Ntrials <- boussard2_data[,interaction(trial,reversal)] %>%
  unique() %>% length()

Nrev <- boussard2_data[,reversal] %>%
  unique() %>% length()

setorder(boussard2_data,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard2_data[,unique(brainsize),by=tankID][,V1]

# set the reversal structure
block_r <-  do.call(rbind,lapply(1:Nind,function(x){
              rbind(matrix(rep(c(0,1),each=24),nrow = 24),
                    matrix(rep(c(1,0),each=42),ncol = 2))
            }))        
# length of the reward block does not fit the dataset
dim(block_r)[1]-
dim(boussard2_data)[1]

# tankID 57 does not have reversal. 
# What shall we do about this??


boussard2_data[,interaction(trial.long,tankID)] %>%
  unique() %>% length()

boussard2_data[,unique(trial)]
boussard2_data[,max(trial),by=.(tankID,reversal)]

table(boussard2_data[,c("tankID","trial.long")])

# Get a unique ID for trials along the reversal blocks
boussard2_data[,trial.long:=interaction(reversal,trial)]
boussard2_data[,trial.long.num:]

# Set NA as failure to choose the rewarding option
boussard2_data[,sum(is.na(success))]

# Transform the success variable, to fit stan model
boussard_wide<-dcast(boussard_data[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)


# using cmdstanr

library(cmdstanr)
# Set cmdstan path
## Erase to run in cluster
stan_path <- here("..","cmdstan","cmdstan-2.32.2")
set_cmdstan_path(stan_path)

# Compile stan model with different random alpha for each reversal
boussard_RW_rev<-cmdstan_model("boussard_RW_rev.stan")


# sample from posterior
fit_boussard_RW_rev <- boussard_RW_rev$sample(list(N=Nind,B=Ntreat,Tr=Ntrials/Nrev,
                                                   Rev=Nrev,TotTr=Ntrials,
                                                   block_r=block_r,
                                                   treat_ID=treat_Inds,
                                                   y=boussard_wide),
                                              parallel_chains = getOption("mc.cores", 5),
                                              chains = 3)

# Save samples to file
fit_boussard_RW_cmd$save_object(file = "fit_boussard_stan.RDS")

fit_boussard_RW_cmd<-readRDS("fit_boussard_stan.RDS")
