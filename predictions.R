library(BayesianTools)
library("Rcpp")
source("RLModel.R")
library(cowplot)


boussard_data <- read_excel(here("data","doi_10.5061_dryad.5mkkwh72s__v2",
                                 "SRL(all_data).xlsx"))

boussard_data <-boussard_data %>% as.data.table()


setorder(boussard_data,reversal,brainsize,tankID)

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

learn.pars <- list(alpha=0.1,temp=2)

learned.ind <- learning(prediction_ind = prediction.ind,learnPar = learn.pars,
                        seed=2)

prediction.ind2 <- as.data.table(prediction.ind)

modelPArs<-c(rbeta(96,1,0.8),1,rbeta(2,1,0.8),0.2)

# str(prediction.ind)
# str(learned.ind)
# 
# start_time <- Sys.time()
# likelihood(modelPArs,seed = 2)
# end_time <- Sys.time()
# end_time - start_time
# 
# prediction.ind <-as.data.table(prediction.ind)
# start_time2 <- Sys.time()
# set.seed(2)
# likelihood2(modelPArs)
# end_time2 <- Sys.time()
# 
# end_time2 - start_time2

set.seed(2)
for(i in 1:(dim(boussard_data.tankID.1)[1])){
  tmp<-  update(as.numeric(prediction.ind2[trial_long==i-1,c("val.1","val.2")]),
                reward = 
                  as.numeric(prediction.ind2[trial_long==i,c("rew.1","rew.2")]),
                learnPars = learn.pars)

  prediction.ind2[trial_long==i,c('val.1','val.2','choice')] <- tmp
  prediction.ind2[trial_long==i,'success'] <- 
    with(prediction.ind2[trial_long==i,c("rew.1","rew.2","choice")],{
      ifelse(choice==1,rew.1,rew.2)
    })
}

loop <-prediction.ind2 %>% 
  mutate(prob_2= ChoiceSoftMax(val.2,val.1,learn.pars$temp,prob=TRUE)) %>% #1/(1+exp(-learnParam$temp*(val.2-val.1)))) %>% 
  as.data.table() %>% 
  melt(id.vars = c('trial_reversal','reversal'),
       measure.vars = c('val.1','val.2','prob_2')) %>% 
    ggplot(aes(y=value,x=trial_reversal,col=variable))+
     geom_line()+
     ylim(0,1)+
     facet_grid(~reversal)+
     theme_classic()

cpp <-learned.ind %>% as.data.table() %>%  
  mutate(prob_2= ChoiceSoftMax(val.2,val.1,learn.pars$temp,prob=TRUE)) %>% #1/(1+exp(-learnParam$temp*(val.2-val.1)))) %>% 
  melt(id.vars = c('trial_reversal','reversal'),
       measure.vars = c('val.1','val.2','prob_2')) %>% 
  ggplot(aes(y=value,x=trial_reversal,col=variable))+
  geom_line()+
  ylim(0,1)+
  facet_grid(~reversal)+
  theme_classic()

plot_grid(loop,cpp,nrow = 2) 
