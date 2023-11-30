data {
  int<lower=0> N; // number of individuals 
  int<lower=0> Tr; // number of trials in total
  array[N,Tr] int <lower=-1,upper=1> reward;    // Reward on each trial
  array[N,Tr] int <lower=0, upper=2> choice;    // choice
}

transformed data {
  // Default value for (re-)initializing parameter vectors
  vector[2] initV;
  initV = rep_vector(0.0, 2);
}

parameters {
  // Individual parameters
  vector[N] alphasID;
  vector[N] tauID;
  // hyper-parameters
  vector[2] mu;
  vector<lower=0>[2] sigma;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector <lower=0, upper=1>[N]  alphasID_t;
  vector <lower=0>[N]  tauID_t;
    for (i in 1:N) {
    alphasID_t[i]  = Phi_approx(mu[1] + sigma[1]*alphasID[i]);
    tauID_t[i] = Phi_approx(mu[2] + sigma[2] * tauID[i]) * 10;
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
  sigma ~ normal(0, 0.2);
  
  // Individual parameters
  alphasID  ~ normal(0, 1);
  tauID  ~ normal(0, 1);
  
  
  for (ind in 1:N){
    // Initialize with 0 estimated values
    est_values = initV;
    
    for (tr in 1:Tr){
      // Calculate probabilities based on values
      probs[2] = inv_logit(-(est_values[2]-est_values[1])*tauID_t[ind]);
      probs[1] = 1-probs[2];
      
      choice[ind,tr]-1 ~ bernoulli(probs[2]); 
      pred_error = reward[ind,tr]  - est_values[choice[ind,tr]];
      est_values[choice[ind,tr]] += alphasID_t[ind]*pred_error;
      
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
      probs_p[2] = inv_logit(-(est_values_p[2]-est_values_p[1])*tauID_t[ind]);
      probs_p[1] = 1-probs_p[2];
      
      log_lik[ind] += bernoulli_lpmf(choice[ind, tr]-1 | probs_p[2]);
      y_pred[ind,tr] = bernoulli_rng(probs_p[2])+1;
      
      pred_error_p = reward[ind,tr]  - est_values_p[y_pred[ind,tr]];
      est_values_p[y_pred[ind,tr]] += alphasID_t[ind]*pred_error_p;
      
      
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block
