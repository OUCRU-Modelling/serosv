data{
  int<lower=1> Nage;           // number of age groups
  real<lower=0.5> init_se;     // best guess for sensitivity
  real<lower=0.5> init_sp;     // best guess for specificity
  int<lower=0> study_size;     // study size for sensitivity, specificity estimate
  int<lower=0> posi[Nage];     // Number of positive cases
  int<lower=0> ni[Nage];       // Number of trials
}

parameters{
  real<lower=0,upper=1> est_se; // estimated sensitivity
  real<lower=0,upper=1> est_sp; // estimated specificity
  real<lower=0,upper=1> theta[Nage]; //corrected seroprevalence, denote as theta for consistency
}

model{
  vector[Nage] apparent_theta; // denote seroprevalence as theta for consistency
  // prior
  est_se ~ beta(init_se*study_size, (1-init_se)*study_size);
  est_sp ~ beta(init_sp*study_size, (1-init_sp)*study_size);


  // posterior
  for (i in 1:Nage){
      theta[i] ~ beta(1, 1); // non informative seroprevalence
      // compute apparent theta = theta*sensitivity + (1-theta)*specificity
      apparent_theta[i] = theta[i]*est_se + (1-theta[i])*(1-est_sp);

      // likelihood
      posi[i] ~ binomial(ni[i], apparent_theta[i]);

      // TODO: include condition(s) mentioned in the paper for relations between sp and se
  }
}
