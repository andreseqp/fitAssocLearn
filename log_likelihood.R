likelihood <- function(param,sum=TRUE){
  ## Parameters:
  ## 1-43 -> alphas brainsize 0
  ## 44-96 -> alphas brainsize 1
  ## 97 -> temp
  ## 98 -> overall alpha brainsize 0
  ## 99 -> overall alpha brainsize 1
  ## 100 -> sig random
  alpha.brain <- c(param[98],param[99])
  prediction<-
    do.call(rbind,lapply(boussard_data[,unique(tankID)],
                       FUN =  function(ind){
    ind.par<-list(alpha=
                    alpha.brain[boussard_data[tankID==ind,
                                                   unique(brainsize)+1]]+param[ind],
                   temp=param[97]) 
    prediction.ind <- learning(prediction_ind = prediction.ind,learnPar = ind.par)
     #  for(i in 1:(dim(prediction.ind)[1])){
     #    tmp<-  update(as.numeric(prediction.ind[trial_long==i-1,c("val.1","val.2")]),
     #                 reward = 
     #                   as.numeric(prediction.ind[trial_long==i,c("rew.1","rew.2")]),
     #                 learnPars = ind.par)
     #   
     #   prediction.ind[trial_long==i,`:=`(val.1=tmp$val.1,val.2=tmp$val.2,
     #                                     choice=tmp$choice)]
     # }
    prediction.ind <- as.data.table(prediction.ind)
     return(prediction.ind[trial_reversal>0])
                       }))
  prediction[,`:=`(success=ifelse(choice==1,rew.1,rew.2),
                   prob= ChoiceSoftMax(val.1,val.2,param[97],prob=TRUE))]
  llRandom <- dnorm(param[1:96],sd = param[100],log=T)
  llObservation <- dbinom(boussard_data[!is.na(success),success],size=1,
                          prob = prediction[!is.na(boussard_data$success),prob],
                          log = TRUE)
  ifelse(sum==TRUE,return(sum(c(llRandom,llObservation))),
         return(c(llRandom,llObservation)))
}

likelihood2 <- function(par,sum=TRUE){
  ## Parameters:
  ## 1-43 -> alphas brainsize 0
  ## 44-96 -> alphas brainsize 1
  ## 97 -> temp
  ## 98 -> overall alpha brainsize 0
  ## 99 -> overall alpha brainsize 1
  ## 100 -> sig random
  alpha.brain <- c(par[98],par[99])
  prediction<-
    do.call(rbind,lapply(boussard_data[,unique(tankID)],
     FUN = function(ind){
       ind.par<-list(alpha=alpha.brain[boussard_data[tankID==ind,
                                                     unique(brainsize)+1]]+par[ind],
                     temp=par[97]) 
       lapply(1:(dim(prediction.ind)[1]),
              FUN = function(i){
                # for(i in 1:(dim(prediction.ind)[1])){
                
                # i<-1
                # prediction.tankID.1[i,c('val.1','val.2','choice')]<-
                tmp<-  update(as.numeric(prediction.ind[trial_long==i-1,c("val.1","val.2")]),
                              reward = 
                                as.numeric(prediction.ind[trial_long==i,c("rew.1","rew.2")]),
                              learnPars = ind.par)
                
                prediction.ind[trial_long==i,c('val.1','val.2','choice')] <- tmp
                prediction.ind[trial_long==i,'success'] <- 
                  with(prediction.ind[trial_long==i,c("rew.1","rew.2","choice")],{
                    ifelse(choice==1,rew.1,rew.2)
                  })
              })
       return(prediction.ind[trial_reversal>0])
     }))
  prediction[,`:=`(success=ifelse(choice==1,rew.1,rew.2),
                   prob= ChoiceSoftMax(val.1,val.2,par[97],prob=TRUE))]
  llRandom <- dnorm(par[1:96],sd = par[100],log=T)
  llObservation <- dbinom(boussard_data[!is.na(success),success],size=1,
                          prob = prediction[!is.na(boussard_data$success),prob],
                          log = TRUE)
  ifelse(sum==TRUE,return(sum(c(llRandom,llObservation))),
         return(c(llRandom,llObservation)))
}