data {
  int<lower=0> N; // number of individuals 
  int<lower=0> B; // number of treatments
  int<lower=0> Tr; // number of trials in total
  array[Tr,2] int <lower=0,upper=1> block_r;    // Reward on each trial
  array[N] int treat_ID;    // treatment ID for the second option
  array[N,Tr] int <lower=0, upper=1> y;    // choice
}

transformed data {
  // Default value for (re-)initializing parameter vectors
  vector[2] initV;
  initV = rep_vector(0.0, 2);
}

parameters {
  // Individual parameters
  vector[N] alphasID;
  vector [N] tauID;
  // treatment parameres
  vector[B] alphasT;
  // hyper-parameters
  vector [2] mu;
  vector <lower=0>[2] sigma;
}
// transformed parameters {
//   // Transform subject-level raw parameters
//   vector <lower=0, upper=1>[N]  alphasID_t;
//   real tau_t;
//   for (i in 1:N) {
//     alphasID_t[i]  = Phi_approx(alphasT[treat_ID[i]+1] 
//                                       + alphasID[i]);
//   }
//   tau_t = Phi_approx(tau)*20;
// }
transformed parameters {
  // Transform subject-level raw parameters
  vector <lower=0, upper=1>[N]  alphasID_t;
  vector <lower=0>[N]  tauID_t;
  for (i in 1:N) {
    alphasID_t[i]  = inv_logit(mu[1] + alphasT[treat_ID[i]+1] 
                                      + sigma[1]*alphasID[i]);
    tauID_t[i] = exp(mu[2]+ sigma[2]*tauID[i]);
  }
}
model {
  // Prediction error
  real pred_error;
  // Probabilities for both options
  vector[2] probs;
  // Estimated values
  vector[2] est_values; 
  
  // Hyperparameters
  mu  ~ normal(0, 1.0);
  //print("target = ", target());
  sigma ~ normal(0, 0.2);
  //print("target = ", target());
  
  // Individual parameters
  // alphasT  ~ normal(mu_alpha, sigma_a);
  alphasT  ~ normal(0, 1);
  alphasID  ~ normal(0, 1);
  tauID ~ normal(0,1);

  for (ind in 1:N){
    // Initialize with 0 estimated values
    est_values = initV;
    
    for (tr in 1:Tr){
      // Calculate probabilities based on values
      probs[2] = inv_logit((est_values[2]-est_values[1])*tauID_t[ind]);
      probs[1] = 1-probs[2];
      
      // Aquí vamos arreglando problemitas
      // choice[ind,tr]-1 ~ bernoulli(probs[2]); 
      // pred_error = reward[ind,tr]  - est_values[choice[ind,tr]];
      // est_values[choice[ind,tr]] += alphasID_t[ind]*pred_error;
      
      if(block_r[tr,2]==1) { // if the second option yields reward
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
      else { // if the first option yields reward
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
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block

generated quantities {

  real pred_error_p;
  vector[2] probs_p;
  vector[2] est_values_p;
  array[N] real log_lik;
  array[N,Tr] int y_pred;


  // Set all posterior predictions to 0 (avoids NULL values)
  for (ind in 1:N) {
    for (tr in 1:Tr) {
      y_pred[ind, tr] = -1;
    }
  }

  for (ind in 1:N){
    est_values_p = initV;
    log_lik[ind] = 0;
    for (tr in 1:Tr){
      probs_p[2] = inv_logit((est_values_p[2]-est_values_p[1])*tauID_t[ind]);
      probs_p[1] = 1-probs_p[2];
      if(block_r[tr,2]==1) {
        log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[2]);
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
        log_lik[ind] += bernoulli_lpmf(y[ind, tr] | probs_p[1]);
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
