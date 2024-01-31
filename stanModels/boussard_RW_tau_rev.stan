// Rescola-wagner model to fit the Boussard et al data set
// the model assumes the speed of learning (alpha) is given by a 
// deterministic effect (brain size treatment)
// and a random effect. The temperature (tau) is the same for all individuals
data {
  int<lower=0> N; // number of individuals 
  int<lower=0> B; // number of treatments
  int<lower=0> Tr; // number of trials per reversal
  int<lower=0> Rev; // number of reversal blocks
  int<lower=0> TotTr; // number of total trials
  array[TotTr,2] int <lower=0,upper=1> block_r;    // Reward on each trial
  array[N] int treat_ID;    // treatment ID for the second option
  array[N,TotTr] int <lower=0, upper=1> y;    // choice
}

transformed data {
  // Default value for (re-)initializing parameter vectors
  vector[2] initV;
  initV = rep_vector(0.0, 2);
}

parameters {
  // Individual parameters
  vector[N] tausID;
  vector[Rev] tausRev;
  real alpha;
  // treatment parameres
  vector[B] tausT;
  // hyper-parameters
  real mu_tau;
  real <lower=0> sigma_t;
}

transformed parameters {
  // Transform subject-level raw parameters
  array[N,Rev] real <lower=0>  tausID_t;
  real<lower=0,upper=1> alpha_t;
  for (i in 1:N) {
    for(j in 1:Rev){
      tausID_t[i,j]  = exp(mu_tau + tausT[treat_ID[i]+1] 
                                      + sigma_t*tausID[i]+tausRev[j]); 
    }
  }
  alpha_t = inv_logit(alpha);
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
  
  
  // Hyperparameters
  mu_tau  ~ normal(0, 10);
  //print("target = ", target());
  sigma_t ~ normal(0, 10);
  //print("target = ", target());
  tausRev ~ normal(0,10);
  
  // Individual parameters
  // alphasT  ~ normal(mu_alpha, sigma_a);
  tausT  ~ normal(0, 10);
  tausID  ~ normal(0, 10);
    // group parameters
  alpha ~ normal(0, 10);

  for (ind in 1:N){
    // Initialize with 0 estimated values
    est_values = initV;
    for (revBlock in 1:Rev){
      
      for (tr in 1:Tr){
        // id of the total trial
        totTrial = (revBlock-1)*Tr+tr;
        // Calculate probabilities based on values
        probs[2] = 
          inv_logit((est_values[2]-est_values[1])*tausID_t[ind,revBlock]);
        probs[1] = 1-probs[2];
        
        // choice[ind,tr]-1 ~ bernoulli(probs[2]); 
        // pred_error = reward[ind,tr]  - est_values[choice[ind,tr]];
        // est_values[choice[ind,tr]] += alphasID_t[ind]*pred_error;
        
        if(block_r[totTrial,2]==1) { // if the second option yields reward
          // use probability of second option to calculate loglikelihood of success
          y[ind,totTrial] ~ bernoulli(probs[2]); 
          
          if(y[ind,totTrial]==1){ // if choice was succesfull
            pred_error = block_r[totTrial,2]  - est_values[2];
            est_values[2] += alpha_t*pred_error;
            // update action values for option 2
          }
          else{ // if choice was unsuccesfull
            pred_error =  block_r[totTrial,1]- est_values[1];
            est_values[1] += alpha_t*pred_error;
            // update action values for option 1
          }
        }
        else { // if the first option yields reward
          // use probability of first option to calculate loglikelihood of success
          y[ind,totTrial] ~ bernoulli(probs[1]);
          
          if(y[ind,totTrial]==1){ // if choice was successful
            pred_error =  block_r[totTrial,1] - est_values[1];
            est_values[1] += alpha_t*pred_error;
            // update action values for option 1
          }
          else{ // if choice was unsuccesfull
            pred_error = block_r[totTrial,2] - est_values[2];
            est_values[2] += alpha_t*pred_error;
            // update action values for option 2
          }
        }
      } // end of trials loop
    } // end of the reversals loop
  } // end of inddividuals loop
} // end of model block

generated quantities {
  int totTrial_p;
  real pred_error_p;
  vector[2] probs_p;
  vector[2] est_values_p;
  array[N] real log_lik;
  array[N,TotTr] int y_pred;


  // Set all posterior predictions to 0 (avoids NULL values)
  for (ind in 1:N) {
    for (tr in 1:TotTr) {
      y_pred[ind, tr] = -1;
    }
  }

  for (ind in 1:N){
    est_values_p = initV;
    log_lik[ind] = 0;
    for (revBlock in 1:Rev){
      for (tr in 1:Tr){
        totTrial_p = (revBlock-1)*Tr+tr;
        probs_p[2] = 
          inv_logit((est_values_p[2]-est_values_p[1])*tausID_t[ind,revBlock]);
        probs_p[1] = 1-probs_p[2];
        if(block_r[totTrial_p,2]==1) {
          log_lik[ind] += bernoulli_lpmf(y[ind, totTrial_p] | probs_p[2]);
          y_pred[ind,totTrial_p] = bernoulli_rng(probs_p[2]);
          if(y_pred[ind,totTrial_p]==1){
            pred_error_p = block_r[totTrial_p,2]  - est_values_p[2];
            est_values_p[2] += alpha_t*pred_error_p;
            // update action values
          }
          else{
            pred_error_p =  block_r[totTrial_p,1]- est_values_p[1];
            est_values_p[1] += alpha_t*pred_error_p;
            // update action values
          }
        }
        else {
          log_lik[ind] += bernoulli_lpmf(y[ind, totTrial_p] | probs_p[1]);
          y_pred[ind,totTrial_p] = bernoulli_rng(probs_p[1]);
          if(y_pred[ind,totTrial_p]==1){
            pred_error_p =  block_r[totTrial_p,1] - est_values_p[1];
            est_values_p[1] += alpha_t*pred_error_p;
            // update action values
          }
          else{
            pred_error_p = block_r[totTrial_p,2] - est_values_p[2];
            est_values_p[2] += alpha_t*pred_error_p;
            // update action values
          }
        }
      } //end of reversal loop
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block
