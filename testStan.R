## Test the use of stan
library(here)
library("rstan")
library("shinystan") #install.packages("shinystan")

N<-100
realMean<-1.2
realSigma<-0.4
y<-rnorm(N,realMean,realSigma)

model.1<-stan_model("myfirstStan.stan")

fit.1<- sampling(model.1,list(N=N,y=y),iter=1000,chains=4)

print(fit.1)
launch_shinystan(fit.1)

logist_model<-stan_model("logisticHierar.stan")

boussard_RW<-stan_model("boussard_RW.stan")

## Test the use of cmdstanr
library(here)
library(cmdstanr)
set_cmdstan_path("M:\\Projects\\cmdstan\\cmdstan-2.32.2")


N<-100
realMean<-1.2
realSigma<-0.4
y<-rnorm(N,realMean,realSigma)

model.1<-cmdstan_model("myfirstStan.stan")

fit.1<- model.1$sample(list(N=N,y=y))

print(fit.1)
launch_shinystan(fit.1)

fit.1$save_object(file = "myfirststan.RDS")

fit.2<-readRDS("myfirststan.RDS")

logist_model<-stan_model("logisticHierar.stan")

boussard_RW<-stan_model("boussard_RW.stan")
