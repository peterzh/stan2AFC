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
}
model {
  //likelihood
  choiceR ~ bernoulli_logit( bias + sens_left*contrast_left + sens_right*contrast_right );
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens_left.*contrast_left + p.sens_right.*contrast_right
