library(Rcpp)
library(here)
library("BayesianTools")
library(data.table)
source("RLModel.R")
library(readxl)
library(dplyr)
library(parallel)

## Parameters:
## 1-43 -> alphas brainsize 0
## 44-96 -> alphas brainsize 1
## 97 -> temp
## 98 -> overall alpha brainsize 0
## 99 -> overall alpha brainsize 1
## 100 -> sig random

## Read the data 

boussard_data <- read_excel(here("data","doi_10.5061_dryad.5mkkwh72s__v2",
                                 "SRL(all_data).xlsx"))

boussard_data <-boussard_data %>% as.data.table()

boussard_data.tankID.1 <- boussard_data %>% filter(tankID==1)
boussard_data.tankID.1

prediction.ind <- data.frame(
  val.1 = rep(0,dim(boussard_data.tankID.1)[1]+1),
  val.2 = rep(0,dim(boussard_data.tankID.1)[1]+1),
  choice = rep(0,dim(boussard_data.tankID.1)[1]+1),
  reversal = c(0,boussard_data.tankID.1$reversal),
  rew.1 = c(0,rep(c(0,1),each=30,times=6)[1:dim(boussard_data.tankID.1)[1]]),
  rew.2 = c(1,rep(c(1,0),each=30,times=6)[1:dim(boussard_data.tankID.1)[1]]),
  success = rep(0,dim(boussard_data.tankID.1)[1]+1),
  trial_long = 0:dim(boussard_data.tankID.1)[1],
  trial_reversal = c(0,rep(1:30,11))
)


# define bayesian set-up
bayesianSetup <- createBayesianSetup(likelihood = likelihood, 
                                     lower = c(rep(-1,96),0,-1,-1,0),
                                     upper = c(rep(1,96),100,1,1,50),
                                     best = c(rep(0.1,96),10,0.1,0.11,0.5),
                                     names = c(paste0("alpha",1:96),"temp",
                                               "alpha_SB","alpha_LB","RE_sd")
                                       )

settings <- list(iterations = 10000, nrChains = 1)
res <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

nChains<-5
cl <- parallel::makeCluster(nChains,outfile="out")
parallel::clusterExport(cl,c("boussard_data",
                             "prediction.ind")
)
parallel::clusterEvalQ(cl, {
   library(BayesianTools)
  library(Rcpp)
  library(here)
  library(data.table)
  source("RLModel.R")
  defaultPars<-foc.param
}
)

## calculate parallel n chains, for each chain the likelihood will be calculated on one core
MCMC.alphas <- parallel::parLapply(cl, 1:nChains, 
                                fun = function(X, bayesianSetup, settings) 
                                  runMCMC(bayesianSetup, settings, sampler = "DEzs"), 
                                bayesianSetup, settings)


## Combine the chains
MCMC.alphas <- createMcmcSamplerList(MCMC.FAA)

# Save files for future analysis
saveRDS(MCMC.alphas, file= here(paste0(scenario,"_"),"MCMC_alphas.rda"))
