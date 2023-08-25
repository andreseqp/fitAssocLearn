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
                                     lower = c(rep(-20,96),-20,-20,-20,0),
                                     upper = c(rep(20,96),20,20,20,50),
                                     best = c(rep(0,96),0,0,0,0.5),
                                     names = c(paste0("alpha",1:96),"temp",
                                               "alpha_SB","alpha_LB","RE_sd")
                                       )

settings <- list(iterations = 1000, nrChains = 1, message = TRUE)

likelihood(c(rep(0.7,96),10,0.5,0.5,00.001))
likelihood_plain(c(1,0.5,0.5))

res <- runMCMC(bayesianSetup = bayesianSetup, settings = settings,
               sampler = "DREAM")

pnorm(-5)
par(mfrow=c(1,1))
plot(res)
summary(res)

res$chain[[1]]

likelihoods<-sapply(1:300, function(x){likelihood(res$codaChain[[1]][x,])})
likelihood(as.array(res$codaChain[[3]][1,]))
par(mfrow=c(1,1))
plot(likelihoods)

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
}
)

## calculate parallel n chains, for each chain the likelihood will be calculated on one core
MCMC.alphas <- parallel::parLapply(cl, 1:nChains, 
                                fun = function(X, bayesianSetup, settings) 
                                  runMCMC(bayesianSetup, settings, sampler = "DEzs"), 
                                bayesianSetup, settings)


## Combine the chains
MCMC.alphas <- createMcmcSamplerList(MCMC.alphas)

# Save files for future analysis
saveRDS(MCMC.alphas, file= here("MCMC_alphas.rda"))


# Plot the MCMC chains ----------------------------------------------------------

par()
plot(MCMC.alphas)
summary(MCMC.alphas)
marginalPlot(MCMC.alphas)
gelmanDiagnostics(MCMC.FAA.loaded)

## Let's try to fit a model without the hierarchical structure -----------------
 
# define bayesian set-up
bayesianSetup <- createBayesianSetup(likelihood = likelihood_plain, 
                                     lower = c(-20,-20,-20),
                                     upper = c(20,20,20),
                                     best = c(0,0,0),
                                     names = c("temp",
                                               "alpha_SB","alpha_LB")
)

settings <- list(iterations = 1000, nrChains = 1)

likelihood_plain(c(10,0.5,0.5))

res <- runMCMC(bayesianSetup = bayesianSetup, settings = settings,
               sampler = "SMC")

plot(res)
summary(res)

likelihoods<-sapply(1:334, function(x){likelihood_plain(res$codaChain[[1]][x,])})
sum(likelihoods!=res$codaChain[[1]][,5])
likelihood(as.array(res$codaChain[[3]][1,]))
par(mfrow=c(1,1))
plot(likelihoods)

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
}
)

## calculate parallel n chains, for each chain the likelihood will be calculated on one core
MCMC.alphas <- parallel::parLapply(cl, 1:nChains, 
                                   fun = function(X, bayesianSetup, settings) 
                                     runMCMC(bayesianSetup, settings, sampler = "DEzs"), 
                                   bayesianSetup, settings)


## Combine the chains
MCMC.alphas <- createMcmcSamplerList(MCMC.alphas)

# Save files for future analysis
saveRDS(MCMC.alphas, file= here("MCMC_alphas.rda"))


# Plot the MCMC chains ----------------------------------------------------------

par()
plot(MCMC.alphas)
summary(MCMC.alphas)
marginalPlot(MCMC.alphas)
gelmanDiagnostics(MCMC.FAA.loaded)