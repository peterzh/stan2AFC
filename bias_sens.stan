data {
    int<lower=0> N; // number of observations
    vector[N] contrast_left;
    vector[N] contrast_right;
    int<lower=0,upper=1> choiceR[N];
}
parameters {
    real bias; // bias parameter
    real sens; // sensitivity
}
model {
  //likelihood
  choiceR ~ bernoulli_logit( bias + sens*(contrast_right - contrast_left)) );
}
