data{
  int<lower=1> Nage;           // number of age groups
  real<lower=0.5> init_se;     // best guess for sensitivity
  real<lower=0.5> init_sp;     // best guess for specificity
  int<lower=0> study_size_se;  // study size for sensitivity estimate
  int<lower=0> study_size_sp;  // study size for specificity estimate
  array[Nage] int<lower=0> posi;     // Number of positive cases
  array[Nage] int<lower=0> ni;       // Number of trials
}

parameters{
  real<lower=0,upper=1> est_se; // estimated sensitivity
  real<lower=0,upper=1> est_sp; // estimated specificity
  array[Nage] real<lower=0,upper=1> theta; //corrected seroprevalence, denote as theta for consistency
}

model{
  vector[Nage] apparent_theta; // denote seroprevalence as theta for consistency
  // prior
  est_se ~ beta(init_se*study_size_se, (1-init_se)*study_size_se);
  est_sp ~ beta(init_sp*study_size_sp, (1-init_sp)*study_size_sp);


  for (i in 1:Nage){
      theta[i] ~ beta(1, 1); // non informative seroprevalence
      // compute apparent theta = theta*sensitivity + (1-theta)*specificity
      apparent_theta[i] = theta[i]*est_se + (1-theta[i])*(1-est_sp);

      // likelihood
      posi[i] ~ binomial(ni[i], apparent_theta[i]);

      // TODO: include condition(s) mentioned in the paper for relations between sp and se
  }
}
