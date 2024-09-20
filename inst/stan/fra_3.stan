data {
  int<lower=1> Nage;           // Number of age groups
  vector[Nage] age;            // Age vector
  int<lower=0> posi[Nage];     // Number of positive cases
  int<lower=0> ni[Nage];       // Number of trials
}

parameters {
  real<lower=0.00001> alpha1;  // Parameter alpha1
  real<lower=0.00001> alpha2;  // Parameter alpha2
  real<lower=0.00001> alpha3;  // Parameter alpha3
  real<lower=0> tau_alpha1;    // Precision (inverse variance) of alpha1
  real<lower=0> tau_alpha2;    // Precision (inverse variance) of alpha2
  real<lower=0> tau_alpha3;    // Precision (inverse variance) of alpha3
  real mu_alpha1;              // Mean of normal prior for alpha1
  real mu_alpha2;              // Mean of normal prior for alpha2
  real mu_alpha3;              // Mean of normal prior for alpha3
}

transformed parameters {
  real<lower=0> sigma_alpha1;  // Standard deviation of alpha1
  real<lower=0> sigma_alpha2;  // Standard deviation of alpha2
  real<lower=0> sigma_alpha3;  // Standard deviation of alpha3
 
  sigma_alpha1 = sqrt(1 / tau_alpha1);
  sigma_alpha2 = sqrt(1 / tau_alpha2);
  sigma_alpha3 = sqrt(1 / tau_alpha3);
}

model {
  vector[Nage] theta;

  // Priors
  alpha1 ~ normal(mu_alpha1, sigma_alpha1);
  alpha2 ~ normal(mu_alpha2, sigma_alpha2);
  alpha3 ~ normal(mu_alpha3, sigma_alpha3);
  tau_alpha1 ~ gamma(0.01, 0.01);
  tau_alpha2 ~ gamma(0.01, 0.01);
  tau_alpha3 ~ gamma(0.01, 0.01);
  mu_alpha1 ~ normal(0, sqrt(1 / 0.0001));
  mu_alpha2 ~ normal(0, sqrt(1 / 0.0001));
  mu_alpha3 ~ normal(0, sqrt(1 / 0.0001));

  // Likelihood
  for (i in 1:Nage) {
    real eta;  // Declare eta as a real variable
    eta = (alpha1 / alpha2) * age[i] * exp(-alpha2 * age[i]) +
          (1 / alpha2) * ((alpha1 / alpha2) - alpha3) * (exp(-alpha2 * age[i]) - 1) -
          alpha3 * age[i];
    theta[i] = 1 - exp(eta);
    posi[i] ~ binomial(ni[i], theta[i]);

  }

}
