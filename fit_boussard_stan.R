## Use stan to fit the RM model to the Boussard data


library("shinystan") #install.packages("shinystan")
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

Nind <- boussard_data[,tankID] %>% unique() %>% length()

Ntreat <- boussard_data[,brainsize] %>% unique() %>% length()

Ntrials <- boussard_data[,interaction(trial,reversal)] %>% 
              unique() %>% length()

setorder(boussard_data,tankID,reversal,trial)

treat_Inds <- boussard_data[,unique(brainsize),by=tankID][,V1]  

block_r <-cbind(c(rep(c(0,1),each=30,
                      times=6)[1:dim(boussard_data[tankID==1])[1]]),
                c(rep(c(1,0),each=30,
                      times=6)[1:dim(boussard_data[tankID==1])[1]]))

# block_r <-cbind(c(rep(c(0,1),each=30,
#                       times=2)[1:dim(boussard_data[tankID==1])[1]]),
#                 c(rep(c(1,0),each=30,
#                       times=2)[1:dim(boussard_data[tankID==1])[1]]))

boussard_data[,trial.long:=interaction(reversal,trial)]

boussard_data[is.na(success),success:=0]

boussard_wide<-dcast(boussard_data[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

boussard_wide[,trial.long:=NULL]

boussard_wide<-t(boussard_wide)

str(boussard_wide)

# using cmdstanr

library(cmdstanr)
#set_cmdstan_path("M:\\Projects\\cmdstan\\cmdstan-2.32.2")


boussard_RW_cmd<-cmdstan_model("boussard_RW.stan")

fit_boussard_RW_cmd <- boussard_RW_cmd$sample(list(N=Nind,B=Ntreat,Tr=Ntrials,
                                 block_r=block_r,
                                 treat_ID=treat_Inds,
                                 y=boussard_wide),
                                 parallel_chains = getOption("mc.cores", 5),
                                 chains = 5)

temp_rds_file <- tempfile(fileext = ".RDS")
fit_boussard_RW_cmd$save_object(file = temp_rds_file)

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

