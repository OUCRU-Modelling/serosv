data {
  int<lower=1> Nage;           // Number of age groups
  vector[Nage] age;            // Age vector
  array[Nage] int<lower=0> posi;     // Number of positive cases
  array[Nage] int<lower=0> ni;       // Number of trials
}

parameters {
  real<lower=0.00001> alpha1;                 // Parameter alpha1
  real alpha2;                 // Parameter alpha2
  real<lower=0> tau_alpha1;    // Precision (inverse variance) of alpha1
  real<lower=0> tau_alpha2;    // Precision (inverse variance) of alpha2
  real mu_alpha1;              // Mean of normal prior for alpha1
  real mu_alpha2;              // Mean of normal prior for alpha2
}

transformed parameters {
  real<lower=0> sigma_alpha1;  // Standard deviation of alpha1
  real<lower=0> sigma_alpha2;  // Standard deviation of alpha2
 
  sigma_alpha1 = sqrt(1 / tau_alpha1);
  sigma_alpha2 = sqrt(1 / tau_alpha2);
}

model {
  vector[Nage] theta; 
  for (i in 1:Nage) {
    theta[i] = inv_logit(alpha2 + alpha1 * log(age[i]));
  }
  // Likelihood
  for (i in 1:Nage) {
    posi[i] ~ binomial(ni[i], theta[i]);
  }
 
  // Priors
  alpha1 ~ normal(mu_alpha1, sigma_alpha1);  // Adjust these priors as needed
  alpha2 ~ normal(mu_alpha2, sigma_alpha2);  // Adjust these priors as needed
  tau_alpha1 ~ gamma(0.01,0.01);
  tau_alpha2 ~ gamma(0.01,0.01);
  mu_alpha1 ~ normal(0, sqrt(10000));
  mu_alpha2 ~ normal(0, sqrt(10000));
}
