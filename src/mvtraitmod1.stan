data {
  //integers and indices
  int<lower=0> N; //total sample size
  int<lower=0> T; //total traits
  int<lower=0> Td; //total defense traits
  int<lower=0> Tf; //total foraging traits
  int<lower=0> D; //total drainages
  int<lower=0> L; //total lakes
  int<lower=0> Fd; //number of foraging fixed effects
  int<lower=0> Ff; //number of defense fixed effects
  array[N] int charr_pres; //indicator of charr presence
  array[N] int drain_id; //index of drainage per observation
  array[N] int lake_id; //index of lake per observation
  
  //predictors
  matrix[N, Fd] Xd; //design matrix for defensive traits
  matrix[N, Ff] Xf; //design matrix for foraging traits
  
  //traits
  array[N] vector[T] traits;
}

parameters {
  //fixed effects
  matrix[Fd,Td] beta_defense; //fixed effects for defensive traits
  matrix[Ff,Tf] beta_foraging; //fixed effects for foraging traits
  
  //random effects
  vector<lower=0>[T] sd_drain; //drainage scale
  vector<lower=0>[T] sd_lake; //lake scale
  vector<lower=0>[T] sd_obs; //observation-level (residual) scale
  cholesky_factor_corr[T] LR_obs; //observation-level (residual) corr
  matrix[D, T] z_drain; //unscaled drainage random effects
  matrix[L, T] z_lake; //unscaled lake random effects
}

transformed parameters{
  //random effects
  matrix[D, T] drain = z_drain * diag_matrix(sd_drain);
  matrix[L, T] lake  = z_lake  * diag_matrix(sd_lake);
  matrix[T,T] LP_obs = diag_pre_multiply(sd_obs, LR_obs);
  
  //linear predictor cbind(defense,foraging)
  matrix[N, T] mu;
  mu[, 1:Td] = Xd * beta_defense;
  mu[, (Td+1):T] = Xf * beta_foraging;
  mu += drain[drain_id] + lake[lake_id];
  
}

model {
  //likelihood function
  for(n in 1:N){
    traits[n] ~ multi_normal_cholesky(mu[n], LP_obs); 
  }

  //priors
  to_vector(beta_defense) ~ std_normal();
  to_vector(beta_foraging) ~ std_normal();
  sd_obs ~ exponential(2);
  sd_drain ~ exponential(2);
  sd_lake ~ exponential(2);
  LR_obs ~ lkj_corr_cholesky(2);
  to_vector(z_drain) ~ std_normal();
  to_vector(z_lake) ~ std_normal();
}

generated quantities {
  //trait matrices
  matrix[T, T] R = LR_obs * LR_obs';
  matrix[T, T] P = LP_obs * LP_obs';
  vector[N] log_lik;
  
  //for computing LOO
  for (n in 1:N) {
      log_lik[n] = multi_normal_cholesky_lpdf(traits[n] | 
                      mu[n], LP_obs);
    }
}



