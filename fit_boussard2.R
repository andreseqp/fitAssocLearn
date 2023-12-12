## Fitting a RW model to the data of Boussard et al. 2021
#  Boussard, A., M. Amcoff, S. D. Buechel, A. Kotrschal, and N. Kolm. 2021. 
#  The link between relative brain size and cognitive ageing in female guppies 
# (Poecilia reticulata) artificially selected for variation in brain size. 
#  Experimental Gerontology 146:111218.


## Import important packages
library(data.table)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(shinystan)
library(stringi)
library(tidyverse)
theme_set(new = theme_classic())
## Read the data 

boussard2_data <- read_excel(here("data","boussard_2021",
                                 "Cognitive_ageing_original_data.xlsx"),
                             sheet = 'FL orig. data',
                             col_types = c('text','numeric','numeric',
                                           'numeric','text','text','text',
                                           'text','numeric','numeric'))
boussard2_data <- boussard2_data %>% mutate(reversal=0)

boussard2_data_b <- read_excel(here("data","boussard_2021",
                                  "Cognitive_ageing_original_data.xlsx"),
                             sheet = 'RL orig. data',
                             col_types = c('text','numeric','numeric',
                                           'numeric','text','text','text',
                                           'text','numeric','numeric'))

boussard2_data_b <- boussard2_data_b %>% mutate(reversal=1)

boussard2_data <- as.data.table(boussard2_data)

boussard2_data_b <- as.data.table(boussard2_data_b)

boussard2_data$tankID %>% unique() %>% length()

boussard2_data <-boussard2_data[tankID!="57"]

boussard2_data_b$tankID %>% unique() %>% length()

boussard2_data <- rbind(boussard2_data,boussard2_data_b) 

boussard2_data[,unique(brainsize)]
boussard2_data[,unique(age)]

boussard2_data[,age.int:=lapply(age,function(x){
  switch(x,young=1,middle=2,old=3)})]
boussard2_data[,brainsize:=brainsize+1]

# tankID coincides with the number of fish supposedly tested in the experimental 
# set up

interaction(boussard2_data$tankID,boussard2_data$fishID) %>% 
  unique() %>% length()
table(boussard2_data[,c("fishID","replicate")])
table(boussard2_data[,c("tankID","fishID")])
table(boussard2_data[,c("tankID","replicate")])
table(boussard2_data[,c("tankID","rewarded.colour")])
table(boussard2_data[,c("fishID","replicate","tankID")])

# get number of individuals
Nind <- boussard2_data[,tankID] %>% unique() %>% length()

# get number of treatment groups
Ntreat <- 2

# get the number of treatment groups for all treatments
NtreatEff<- boussard2_data[,lapply(.SD,function(x){
  unique(x) %>% length()}),.SDcols = c("brainsize","age")] %>% as.numeric()
  


# get total number of trials (including all reversals)
Ntrials <- boussard2_data[,interaction(trial,reversal)] %>%
  unique() %>% length()

Nrev <- boussard2_data[,reversal] %>%
  unique() %>% length()

setorder(boussard2_data,tankID,reversal,trial)

# get the treatment for all individuals
treat_Inds <- boussard2_data[,.(unique(brainsize),unique(age.int)),by=tankID] %>% 
  as.data.table()
treat_Inds[,tankID:=NULL]
setnames(treat_Inds,c("V1","V2"),c("brainsize","age"))
treat_Inds <- t(as.matrix(treat_Inds))

# set the reversal structure
block_r <-  rbind(matrix(rep(c(0,1),each=24),nrow = 24),
                    matrix(rep(c(1,0),each=42),ncol = 2))
            
# do.call(rbind,lapply(1:Nind,function(x){
#   rbind(matrix(rep(c(0,1),each=24),nrow = 24),
#         matrix(rep(c(1,0),each=42),ncol = 2))
# }))        


# tankID 57 does not have reversal. 
# What shall we do about this??
# For now I will delete individual (tankID) 57

boussard2_data[,trial.long:=reversal*24+trial]
  
# Get a unique ID for trials along the reversal blocks
boussard2_data[,interaction(trial.long,tankID)] %>%
  unique() %>% length()

boussard2_data[,unique(trial)]

# Set NA as failure to choose the rewarding option
boussard2_data[is.na(success),sum(is.na(success))]
boussard2_data[is.na(success),success:=0]

# There is a data point wrongly typed. 
# I will set it to 1
boussard2_data[success==9,]
boussard2_data[success==9,success:=0]


# Transform the success variable, to fit stan model
boussard2_wide<-dcast(boussard2_data[,.(trial.long,tankID,success)],
                     trial.long~tankID,value.var = "success")

treatAssign <- boussard2_data[,.(max(brainsize),max(age)),by=tankID]
treatAssign[,ind:=1:287]
setnames(treatAssign,c("tankID","V1","V2","ind"),c("tankID","brainsize","age","ind"))

boussard2_wide[,trial.long:=NULL]

boussard2_wide<-t(boussard2_wide)

# using cmdstanr

library(cmdstanr)
# Set cmdstan path
## Erase to run in cluster
stan_path <- here("..","cmdstan","cmdstan-2.32.2")
set_cmdstan_path(stan_path)

## Let's first fit a simple RW model
# Compile stan model with different random alpha for each reversal
boussard_RW_2treat<-cmdstan_model("stanModels/boussard_RW_2_treat.stan")


# sample from posterior
fit_boussard_RW_2treat <- boussard_RW_2treat$sample(list(N=Nind,B=Ntreat,
                                       BE= NtreatEff,
                                       Tr=Ntrials,
                                       block_r=block_r,
                                       treat_ID=treat_Inds,
                                       y=boussard2_wide),
                                    parallel_chains = getOption("mc.cores", 5),
                                    chains = 3)

# Save samples to file
fit_boussard_RW_2treat$save_object(file = "fit_boussard2_stan2treat.RDS")

fit_boussard_RW_2treat<-readRDS("fit_boussard2_stan2treat.RDS")

# Use shinystan to evaluate the performance of the model
# launch_shinystan(fit_boussard)


pars2plot <- c("tau", "mu_alpha", "alphasT[1]", "alphasT[2]",
               "alphasT[3]", "alphasT[4]","alphasT[5]",
                "sigma_a")

posteriors<-fit_boussard_RW_2treat$draws(
  variables = pars2plot)

dimnames(posteriors)$variable <- c(expression(tau),expression(mu[alpha]),
                                   "small brain size", "large brain size", 
                                        "young age", "middle age", "old age",
                                   expression(sigma[alpha]))


png("images/boussard2intervals.png")
mcmc_intervals(posteriors)  
dev.off()

preds <- fit_boussard_RW_2treat$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                   measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

# preds_long <- preds_long %>% filter(.chain<2)
# preds_long <- as.data.table(preds_long)

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long <- preds_long[.chain<3]
preds_long <- preds_long[treatAssign[,.(ind,brainsize,age)],on=.(individual=ind)]

# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),by=.(.chain,.iteration,brainsize,age,trial)]

rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= ifelse(trial>24,1,0)) %>%
  mutate(RTrial=reversal*(-24)+trial)

colnames(mean_ind) <- c("chain","iteration","brainsize","age","Ttrial","success",
                        "reversal","trial")


png("images/boussard2_ppchecks.png")
boussard2_data %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=age,fill=age))+
    stat_summary(fun = mean,geom = "point")+
    stat_summary(fun = mean,geom = "line")+
    # geom_point()+
    theme(legend.position = c(0.8,0.5),
          legend.direction = "horizontal",
          strip.text.y = element_blank())+
    guides(fill=guide_legend(title="age"))+
    ggtitle("Repeated reversal vs  brainsize")+
    stat_summary(data=mean_ind,aes(x=trial,y=success,col=age,fill=age),
                 geom="ribbon",alpha = 0.2,fun.max = function(x){
      quantile(x,0.95)},
      fun.min = function(x){
      quantile(x,0.05)})+
    stat_summary(data=mean_ind,aes(x=trial,y=success,col=age,fill=age),
                 geom="ribbon",alpha = 0.5,fun.max = function(x){
      quantile(x,0.75)},
      fun.min = function(x){
      quantile(x,0.25)})+
    facet_grid(brainsize~reversal,scales = "free_x")
dev.off()


# Include the effect of the color used for the association

treat_cols <- boussard2_data[,unique(rewarded.colour),by=.(tankID,reversal)] %>% 
  as.data.table()
treat_cols <- dcast(treat_cols,tankID~reversal,value.var = "V1")
treat_cols[,tankID:=NULL]
setnames(treat_cols,old = c("0","1"),new = c("Initial","Reversal"))
treat_cols[,`:=`(Initial=ifelse(Initial=="yellow",1,2),
                 Reversal=ifelse(Reversal=="Yellow",1,2))]
treat_cols <- as.matrix(treat_cols)

iniTR = boussard2_data[reversal==0,max(trial.long)]

## Let's first fit a simple RW model
# Compile stan model with different random alpha for each reversal
boussard_RW_2treat_col<-cmdstan_model("stanModels/boussard_RW_2_treat_col.stan")

# sample from posterior
fit_boussard_RW_2treat_col <- boussard_RW_2treat_col$sample(list(N=Nind,B=Ntreat,
                                                         BE= NtreatEff,
                                                         Tr=Ntrials,
                                                         iniTr=iniTR,
                                                         color_assn=treat_cols,
                                                         block_r=block_r,
                                                         treat_ID=treat_Inds,
                                                         y=boussard2_wide),
                                                    parallel_chains = getOption("mc.cores", 5),
                                                    chains = 3)

# Save samples to file
fit_boussard_RW_2treat_col$save_object(file = "fit_boussard2_stan2treat_col.RDS")

fit_boussard_RW_2treat_col<-readRDS("fit_boussard2_stan2treat_col.RDS")

# Use shinystan to evaluate the performance of the model
# launch_shinystan(fit_boussard_RW_2treat_col)


pars2plot <- c("tau", "mu_alpha", "alphasT[1]", "alphasT[2]",
               "alphasT[3]", "alphasT[4]","alphasT[5]","colAlphas[1]",
               "colAlphas[2]","sigma_a")

posteriors<-fit_boussard_RW_2treat_col$draws(
  variables = pars2plot)

dimnames(posteriors)$variable <- c(expression(tau),expression(mu[alpha]),
                                   "small brain size", "large brain size", 
                                   "young age", "middle age", "old age",
                                   "yellow", "red",
                                   expression(sigma[alpha]))


png("images/boussard2intervals_col.png")
mcmc_intervals(posteriors)  
dev.off()

preds <- fit_boussard_RW_2treat$draws(variables = "y_pred")

preds_df <- posterior::as_draws_df(preds)

preds_long <- reshape2::melt(preds_df,id=c('.chain','.iteration'),
                             measure.vars=grep("y_pred",colnames(preds_df)))

rm(list=c("preds","preds_df"))

# preds_long <- preds_long %>% filter(.chain<2)
# preds_long <- as.data.table(preds_long)

preds_long <- as.data.table(preds_long)

preds_long[,c("individual","trial"):=tstrsplit(variable,",")]

preds_long[,variable:=NULL]

preds_long[,`:=`(individual=parse_number(individual),
                 trial=parse_number(trial))]

preds_long <- preds_long[treatAssign[,.(ind,brainsize,age)],
                         on=.(individual=ind)]

treat_cols <- as.data.table(treat_cols) 
treat_cols$ind <- 1:287


preds_long <- preds_long[treat_cols,on=.(individual=ind)]

preds_long[,rewarded.colour:=ifelse(trial<=iniTR,
                                    Initial,Reversal)]

# Average over de MCMC samples
mean_ind <- preds_long[,mean(value),
                       by=.(.chain,.iteration,brainsize,age,trial,
                            rewarded.colour)]


rm(list="preds_long")


mean_ind <-mean_ind %>%
  mutate(reversal= ifelse(trial>24,1,0)) %>%
  mutate(RTrial=reversal*(-24)+trial)

colnames(mean_ind) <- c("chain","iteration","brainsize","age","Ttrial",
                        "rewarded.colour","success",
                        "reversal","trial")


png("images/boussard2_ppchecks_colour.png")
boussard2_data %>% mutate(brainsize=as.factor(brainsize)) %>%
  ggplot(aes(y=success,x=trial,col=age,fill=age))+
  stat_summary(fun = mean,geom = "point")+
  stat_summary(fun = mean,geom = "line")+
  # geom_point()+
  theme(legend.position = c(0.8,0.5),
        legend.direction = "horizontal",
        strip.text.y = element_blank())+
  guides(fill=guide_legend(title="age"))+
  ggtitle("Repeated reversal vs  brainsize")+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=age,fill=age),
               geom="ribbon",alpha = 0.2,fun.max = function(x){
                 quantile(x,0.95)},
               fun.min = function(x){
                 quantile(x,0.05)})+
  stat_summary(data=mean_ind,aes(x=trial,y=success,col=age,fill=age),
               geom="ribbon",alpha = 0.5,fun.max = function(x){
                 quantile(x,0.75)},
               fun.min = function(x){
                 quantile(x,0.25)})+
  facet_grid(brainsize~reversal+rewarded.colour,scales = "free_x")
dev.off()








