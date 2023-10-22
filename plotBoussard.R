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

boussard_data <- read_excel(here("data","doi_10.5061_dryad.5mkkwh72s__v2",
                           "SRL(all_data).xlsx"))

boussard_data <-boussard_data %>% as.data.table()

boussard_data$tankID %>% unique() %>% length()

# tankID coincides with the number of fish supposedly tested in the experimental 
# set up

interaction(boussard_data$tankID,boussard_data$fishID) %>% unique() %>% length()
table(boussard_data[,c("fishID","replicate")])
table(boussard_data[,c("tankID","fishID")])
table(boussard_data[,c("tankID","replicate")])
table(boussard_data[,c("fishID","replicate","tankID")])

# No idea what fishID and replicate are. They perfectly contained in fishID

table(boussard_data[,c("colour","brainsize")])
# Colour and brainsize are parfectly counterbalanced

table(boussard_data[,c("colour","brainsize","reversal")])
# however, there is something odd in the last reversal. Where,
# colour is not perfectly counterbalanced. 
# I guess there is mistake assigning the colour treatment. 
# Misbalanced seems to occur in reversal 10


boussard_data %>% filter(reversal>8) %>% select(c(colour,brainsize,reversal)) %>% 
  table()

boussard_data %>% select(c(colour,reversal)) %>% 
  table()

boussard_data %>% filter(reversal==10,brainsize==1)  %>% select(c(colour)) %>%
  table()



# episode_length<-3
# boussard_data <-boussard_data %>% mutate(episode = floor(trial/episode_length))
# 
# unique(boussard_data$episode)
# 
# 
# boussard_data$episode %>% unique
# 
# boussard_data %>% filter(tankID < 10,reversal<5) %>% 
#   ggplot(aes(y=success,x=episode,col=as.factor(tankID)))+
#     stat_summary(fun = mean,geom = "point")+
#     stat_summary(fun = mean,geom = "line")+
#     facet_grid(reversal~.)+
#     theme_classic()

## Let´s try to make a RL model to fit to the data








