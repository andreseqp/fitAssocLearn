## Reinforcement learning model to fit data

# Soft max desicion making process
ChoiceSoftMax<-function(EV1,EV2,temp,prob=FALSE){
  p<-1/(1+exp(-temp*(EV2-EV1)))
  ifelse(prob,
         return(p),
         return(rbinom(1,1,p)+1)  
  )
}
# ChoiceSoftMax(c(1,1),10)

# function to update 
update<-function(EV0,reward,learnPars){
  # set: all the necessary inputs from previous time step
  with(learnPars,{
      EV1<-EV0
      choiceT0 <- ChoiceSoftMax(EV0[2],EV0[1],temp)
      predE <- reward[choiceT0]-EV0[choiceT0]
      EV1[choiceT0] <- EV0[choiceT0] + alpha*predE
      return(list(val.1=EV1[1],val.2=EV1[2],choice=choiceT0))
    }
  )
}

likelihood <- function(par,sum=TRUE,seed=1){
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
                         FUN =  function(ind){
                           ind.par<-list(alpha=alpha.brain[boussard_data[tankID==ind,
                                                                         unique(brainsize)+1]]+par[ind],
                                         temp=par[97]) 
                           prediction.ind <- learning(prediction_ind = prediction.ind,
                                                      learnPar = ind.par,seed = seed)
                           prediction.ind <- as.data.table(prediction.ind)
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

# for(i in 1:(dim(prediction.ind)[1])){
#   tmp<-  update(as.numeric(prediction.ind[trial_long==i-1,c("val.1","val.2")]),
#                 reward = 
#                   as.numeric(prediction.ind[trial_long==i,c("rew.1","rew.2")]),
#                 learnPars = ind.par)
#   
#   prediction.ind[trial_long==i,`:=`(val.1=tmp$val.1,val.2=tmp$val.2,
#                                     choice=tmp$choice)]

learning<-cppFunction(
  '
  Rcpp::DataFrame learning(Rcpp::DataFrame prediction_ind, Rcpp::List learnPar,
    int seed=1){
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);  
    double predE;
    Rcpp::NumericVector val_1 = prediction_ind["val.1"];
    Rcpp::NumericVector val_2 = prediction_ind["val.2"];
    Rcpp::NumericVector rew_1 = prediction_ind["rew.1"];
    Rcpp::NumericVector rew_2 = prediction_ind["rew.2"];
    Rcpp::NumericVector choice;
    int  choicet0; 
    double prob;
    for(int i = 1; i < val_1.size();++i){
      prob = 1/(1+exp(-double(learnPar["temp"])*(val_2[i-1]-val_1[i-1])));
      choicet0 = R::rbinom(1,prob)+1;
      choice.push_back(choicet0); 
      if(choicet0==2) {
        predE = rew_2[i-1] - val_2[i-1];
        val_2[i] = val_2[i-1] + double(learnPar["alpha"])*predE;
        val_1[i] = val_1[i-1];
      }
      else {
        predE = rew_1[i-1] - val_1[i-1];
        val_1[i] = val_1[i-1] + double(learnPar["alpha"])*predE;
        val_2[i] = val_2[i-1];
      }
    }
    choice.push_back(0);
    prediction_ind["choice"] = choice;
    prediction_ind["val.1"] = val_1;
    prediction_ind["val.2"] = val_2;
    prediction_ind["rew.1"] = rew_1;
    prediction_ind["rew.2"] = rew_2;
    return prediction_ind;
  }
  '
)

randombinom <- cppFunction(
  '
  Rcpp::LogicalVector randombinom (int size, int seed){
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);  
    Rcpp::LogicalVector nums;
    for(int i =0 ; i<size;i++){
      nums.push_back(R::rbinom(1,0.5));
    }
    return nums;
  }
  '
)

randombinom(10,3)
set.seed(3)
as.logical(rbinom(10,1,0.5))





