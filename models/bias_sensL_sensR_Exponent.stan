data {
    int<lower=0> N; // number of observations
    vector[N] contrast_left;
    vector[N] contrast_right;
    int<lower=0,upper=1> choiceR[N];
}
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=1> n_exp; //n exponent
}
model {
  //likelihood
  for (n in 1:N)
    choiceR[n] ~ bernoulli_logit( bias + sens_left*(contrast_left[n]^n_exp) + sens_right*(contrast_right[n]^n_exp) );
}
