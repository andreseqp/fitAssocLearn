// Rescola-wagner model to fit the Boussard et al data set
// the model assumes the speed of learning (alpha) is given by a 
// deterministic effect (brain size treatment)
// and a random effect. The temperature (tau) is the same for all individuals
data {
  int<lower=0> N; // number of individuals 
  int<lower=0> B; // number of treatments
  int<lower=0> Tr; // number of trials in total
  // int<lower=0> N_obs;
  // int<lower=0> N_mis;
  // array [2,N_obs] int <lower=1> ii_obs;
  // array [2,N_mis] int <lower=1> ii_mis;
  array[Tr,2] int <lower=0,upper=1> block_r;    // Reward on each trial
  array[N] int treat_ID;    // treatment ID for the second option
  array[N,Tr] int <lower=-1, upper=1> y;    // choice
  // array[N_obs] int <lower=0, upper=1> y_obs;    // choice
}

transformed data {
  // Default value for (re-)initializing parameter vectors
  vector[2] initV;
  initV = rep_vector(0.0, 2);
}

parameters {
  real tau;
  // treatment parameres
  vector[B] alphasT;
  // hyper-parameters
  real mu_alpha;
  real <lower=0> sigma_a;
  // Individual parameter
  vector[N] alphasID;
  // array [N_mis] int <lower=0, upper=1> y_mis;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector <lower=0, upper=1>[N]  alphasID_t;
  // array[N,Tr] real <lower=0, upper=1> y; 
  real tau_t;
  for (i in 1:N) {
    alphasID_t[i]  = inv_logit(mu_alpha + alphasT[treat_ID[i]+1] 
                                      + sigma_a*alphasID[i]);
  }
  tau_t = exp(tau);
  // for (i in 1:N_obs){
  //   y[ii_obs[1,i],ii_obs[2,i]] = y_obs[i];
  // }
  // for (j in 1:N_mis){
  //   y[ii_mis[1,j],ii_mis[2,j]] = y_mis[j];
  // }
  // // y[ii_obs[1,1:N_obs],ii_obs[2,1:N_obs]] = y_obs;
  // // y[ii_mis[1,1:N_mis],ii_mis[2,1:N_mis]] = y_mis;
}
model {
  // Prediction error
  real pred_error;
  // Probabilities for both options
  vector[2] probs;
  // Estimated values
  vector[2] est_values; 
  
  // Hyperparameters
  mu_alpha  ~ normal(0, 2);
  //print("target = ", target());
  sigma_a ~ normal(0, 2);
  //print("target = ", target());
  
  // Individual parameters
  // alphasT  ~ normal(mu_alpha, sigma_a);
  alphasT  ~ normal(0, 2);
  alphasID  ~ normal(0, 2);
    // group parameters
  tau ~ normal(0, 2);

  for (ind in 1:N){
    // Initialize with 0 estimated values
    est_values = initV;
    
    for (tr in 1:Tr){
      // Calculate probabilities based on values
      probs[2] = inv_logit((est_values[2]-est_values[1])*tau_t);
      probs[1] = 1-probs[2];
      
      if(block_r[tr,2]==1) { // if the second option yields reward
        
        if(y[ind,tr]==-1){ // if the data point is missing
          target += bernoulli_lpmf( 1 | probs[2] );
        }
        else{ 
          // use probability of second option to calculate loglikelihood of success
          y[ind,tr] ~ bernoulli(probs[2]); 
          
          if(y[ind,tr]==1){ // if choice was succesfull
            pred_error = block_r[tr,2]  - est_values[2];
            est_values[2] += alphasID_t[ind]*pred_error;
            // update action values for option 2
          }
          else{ // if choice was unsuccesfull
            pred_error =  block_r[tr,1]- est_values[1];
            est_values[1] += alphasID_t[ind]*pred_error;
            // update action values for option 1
          }  
        }
          
      }
      else { // if the first option yields reward
        if(y[ind,tr]==-1){ // if the data point is missing
          target += bernoulli_lpmf( 1 | probs[1] );
        }
        else{
        // use probability of first option to calculate loglikelihood of success
          y[ind,tr] ~ bernoulli(probs[1]);
          if(y[ind,tr]==1){ // if choice was successful
            pred_error =  block_r[tr,1] - est_values[1];
            est_values[1] += alphasID_t[ind]*pred_error;
            // update action values for option 1
          }
          else{ // if choice was unsuccesfull
            pred_error = block_r[tr,2] - est_values[2];
            est_values[2] += alphasID_t[ind]*pred_error;
            // update action values for option 2
          }
        }
      }
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block

generated quantities {

  real pred_error_p;
  vector[2] probs_p;
  vector[2] est_values_p;
  array[N] real log_lik;
  array[N, TotTr] int y_pred = -1;
  // Set all posterior predictions to 0 (avoids NULL values)


  for (ind in 1:N){
    // est_values_p = initV;
    // Set the initial values according to the provided initial probability
    est_values_p = initV;
    log_lik[ind] = 0;
    for (tr in 1:Tr){
      probs_p[2] = inv_logit((est_values_p[2]-est_values_p[1])*tau_t);
      probs_p[1] = 1-probs_p[2];
      if(block_r[tr,2]==1) {
        if(y[ind, tr]==-1){
          log_lik[ind] += bernoulli_lpmf(1 | probs_p[2]);
        }
        else{
          log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[2]);  
        }
        y_pred[ind,tr] = bernoulli_rng(probs_p[2]);
        if(y_pred[ind,tr]==1){
          pred_error_p = block_r[tr,2]  - est_values_p[2];
          est_values_p[2] += alphasID_t[ind]*pred_error_p;
          // update action values
        }
        else{
          pred_error_p =  block_r[tr,1]- est_values_p[1];
          est_values_p[1] += alphasID_t[ind]*pred_error_p;
          // update action values
        }
      }
      else {
        if(y[ind, tr]==-1){
          log_lik[ind] += bernoulli_lpmf(1 | probs_p[1]);
        }
        else{
          log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[1]);
        }
        y_pred[ind,tr] = bernoulli_rng(probs_p[1]);
        if(y_pred[ind,tr]==1){
          pred_error_p =  block_r[tr,1] - est_values_p[1];
          est_values_p[1] += alphasID_t[ind]*pred_error_p;
          // update action values
        }
        else{
          pred_error_p = block_r[tr,2] - est_values_p[2];
          est_values_p[2] += alphasID_t[ind]*pred_error_p;
          // update action values
        }
      }
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block
