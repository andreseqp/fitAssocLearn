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
  real tau;
  // treatment parameres
  vector[B] alphasT;
  // hyper-parameters
  real mu_alpha;
  real <lower=0> sigma_a;
}
transformed parameters {
  // Transform subject-level raw parameters
  vector <lower=0, upper=1>[N]  alphasID_t;
  real tau_t;
  for (i in 1:N) {
    alphasID_t[i]  = Phi_approx(alphasT[treat_ID[i]+1] 
                                      + alphasID[i]);
  }
  tau_t = Phi_approx(tau)*20;
}

model {
  
  real pred_error;
  vector[2] probs;
  vector[2] est_values; 
  // Hyperparameters
  mu_alpha  ~ normal(0, 1.0);
  //print("target = ", target());
  sigma_a ~ cauchy(0, 0.5);
  //print("target = ", target());
  
  // Individual parameters
  alphasT  ~ normal(mu_alpha, sigma_a);
  
  // group parameters
  tau ~ normal(0, 0.2);

  for (ind in 1:N){
    est_values = initV;
    for (tr in 1:Tr){
    probs[2] = inv_logit(-(est_values[2]-est_values[1])*tau_t);
      probs[1] = 1-probs[2];
      if(block_r[tr,2]==1) {
        y[ind,tr] ~ bernoulli(probs[2]);
        if(y[ind,tr]==1){
          pred_error = block_r[tr,2]  - est_values[2];
          est_values[2] += alphasID_t[ind]*pred_error;
          // update action values
        }
        else{
          pred_error =  block_r[tr,1]- est_values[1];
          est_values[1] += alphasID_t[ind]*pred_error;
          // update action values
        }
      }
      else {
        y[ind,tr] ~ bernoulli(probs[1]);
        if(y[ind,tr]==1){
          pred_error =  block_r[tr,1] - est_values[1];
                  //est_values[1];
          est_values[1] += alphasID_t[ind]*pred_error;
          // update action values
        }
        else{
          pred_error = block_r[tr,2] - est_values[2];
          est_values[2] += alphasID_t[ind]*pred_error;
          // update action values
        }
      }
    } // end of trials loop
  } // end of inddividuals loop
} // end of model block
