// Rescola-wagner model to fit the Boussard et al data set
// the model assumes the speed of learning (alpha) is given by a 
// deterministic effect (brain size treatment)
// and a random effect. The temperature (tau) is the same for all individuals
data {
  int<lower=0> N; 
  // number of individuals 
  int<lower=0> B; 
  // number of treatments
  int<lower=0> Tr; 
  // number of trials per reversal
  int<lower=0> Rev; 
  // number of reversal blocks
  int<lower=0> TotTr; 
  // number of total trials
  real <lower=0> initProp;
  // Intitial probability of success
  array[TotTr,2] int <lower=0,upper=1> block_r;    
  // Reward on each trial
  array[N] int treat_ID;    
  // treatment ID for the second option
  array[N,TotTr] int <lower=-1, upper=1> y;    
  // choice (response variable)
}

transformed data {
  // Default value for (re-)initializing parameter vectors
  // vector[2] initV;
  // initV = rep_vector(0.0, 2);
  real logodds_init = logit(initProp);
}

parameters {
  // Individual learning rates
  vector[N] alphasID;
  // Reversal specific effect on learning rates
  vector[Rev] alphasRev;
  // Temparature parameter
  real tau;
  // Treatment specific effect on learning rates
  vector[B] alphasT;
  // hyper-parameters of random effects:
  // mean effect on the speed of learning
  real mu_alpha;
  // variance of the individual effect on learning rates
  real <lower=0> sigma_a;
}

transformed parameters {
  // Transform subject-level raw parameters to the range required in the RL model
  // Learning rates between 0 and 1
  array[N,Rev] real <lower=0, upper=1>  alphasID_t;
  // Temperature above 0
  real tau_t <lower=0>;
  for (i in 1:N) {
    for(j in 1:Rev){
      alphasID_t[i,j]  = inv_logit(mu_alpha + alphasT[treat_ID[i]+1] 
                                      + sigma_a*alphasID[i]+alphasRev[j]); 
    }
  }
  tau_t = exp(tau);
}
model {
  // total trials
  int totTrial;
  // Prediction error
  real pred_error;
  // Probabilities for both options
  vector[2] probs;
  // Estimated values
  vector[2] est_values; 
  // reversal block
  
  
  // Hyperparameter priors
  mu_alpha  ~ normal(0, 10);
  sigma_a ~ normal(0, 10);
  alphasRev ~ normal(0,10);
  
  // Individual parameters
  alphasT  ~ normal(0, 10);
  alphasID  ~ normal(0, 10);
    // group parameters
  tau ~ normal(0, 10);

  for (ind in 1:N){
    // Initialize estimated values according to the initial prob
    est_values[2] = 0;
    est_values[1] = est_values[2] - logodds_init/tau_t;
    for (revBlock in 1:Rev){
      
      for (tr in 1:Tr){
        // id of the total trial
        totTrial = (revBlock-1)*Tr+tr;
        // Calculate probabilities based on values
        probs[2] = inv_logit((est_values[2]-est_values[1])*tau_t);
        probs[1] = 1-probs[2];
        
        if(block_r[totTrial,2]==1) { 
          // if the second option yields reward
          // use probability of second option to calculate loglikelihood of success
          if(y[ind,tr]==-1){ // if the data point is missing
            target += bernoulli_lpmf( 1 | probs[2] );
          }
          else{
            y[ind,totTrial] ~ bernoulli(probs[2]); 
            
            if(y[ind,totTrial]==1){ 
              // if choice was succesfull
              pred_error = block_r[totTrial,2]  - est_values[2];
              est_values[2] += alphasID_t[ind,revBlock]*pred_error;
              // update action values for option 2
            }
            else{ 
              // if choice was unsuccesfull
              pred_error =  block_r[totTrial,1]- est_values[1];
              est_values[1] += alphasID_t[ind,revBlock]*pred_error;
              // update action values for option 1
            }
          }
        }
        else { 
          // if the first option yields reward
          // use probability of first option to calculate loglikelihood of success
          if(y[ind,tr]==-1){ // if the data point is missing
          target += bernoulli_lpmf( 1 | probs[1] );
        }
        else{
            y[ind,totTrial] ~ bernoulli(probs[1]);
            
            if(y[ind,totTrial]==1){ // if choice was successful
              pred_error =  block_r[totTrial,1] - est_values[1];
              est_values[1] += alphasID_t[ind,revBlock]*pred_error;
              // update action values for option 1
            }
            else{ 
              // if choice was unsuccesfull
              pred_error = block_r[totTrial,2] - est_values[2];
              est_values[2] += alphasID_t[ind,revBlock]*pred_error;
              // update action values for option 2
            }
          }
        } 
      } // end of trials loop
    } // end of the reversals loop
  } // end of inddividuals loop
} // end of model block

generated quantities {
  // Block usted to run the model and get posterior predictive checks
  int totTrial_p;
  real pred_error_p;
  vector[2] probs_p;
  vector[2] est_values_p;
  array[N] real log_lik;
  // Set all posterior predictions to 0 (avoids NULL values)
  array[N, TotTr] int y_pred = -1;


  for (ind in 1:N){
    // est_values_p = initV;
    // Set the initial values according to the provided initial probability
    est_values_p[2] = 0;
    est_values_p[1] = est_values_p[2] - logodds_init/tau_t;
    log_lik[ind] = 0;
    for (revBlock in 1:Rev){ 
      // Reversal loop
      for (tr in 1:Tr){
        // trial loop
        totTrial_p = (revBlock-1)*Tr+tr;
        // calculate the total trial id
        probs_p[2] = inv_logit((est_values_p[2]-est_values_p[1])*tau_t);
        probs_p[1] = 1-probs_p[2];
        // Calculate the probabilities of ech option
        if(block_r[totTrial_p,2]==1) {
           // if the second option yields reward
          if(y[ind, tr]==-1){
          log_lik[ind] += bernoulli_lpmf(1 | probs_p[2]);
          }
          else{
            log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[2]);  
          }
          y_pred[ind,totTrial_p] = bernoulli_rng(probs_p[2]);
          // Draw choice according to the probbility
          if(y_pred[ind,totTrial_p]==1){
            // if choice was succesfull
            pred_error_p = block_r[totTrial_p,2]  - est_values_p[2];
            est_values_p[2] += alphasID_t[ind,revBlock]*pred_error_p;
            // update action values
          }
          else{
            pred_error_p =  block_r[totTrial_p,1]- est_values_p[1];
            est_values_p[1] += alphasID_t[ind,revBlock]*pred_error_p;
            // update action values
          }
        }
        else {
          // If the first option yields reward
          if(y[ind, tr]==-1){
          log_lik[ind] += bernoulli_lpmf(1 | probs_p[1]);
          }
          else{
            log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[1]);
          }
          y_pred[ind,totTrial_p] = bernoulli_rng(probs_p[1]);
          // Draw choice according to the respective probability
          if(y_pred[ind,totTrial_p]==1){
            // If the choice was succesfull
            pred_error_p =  block_r[totTrial_p,1] - est_values_p[1];
            est_values_p[1] += alphasID_t[ind,revBlock]*pred_error_p;
            // update action values
          }
          else{
            pred_error_p = block_r[totTrial_p,2] - est_values_p[2];
            est_values_p[2] += alphasID_t[ind,revBlock]*pred_error_p;
            // update action values
          }
        }
      } //end of reversal loop
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block
