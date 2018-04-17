#include standard_data_blocks_2AFC.stan
parameters {
    real bias; // bias parameter
    real sens; // sensitivity
}
transformed parameters {
  vector[numTrials] z;
  z = bias + sens*(contrastRight - contrastLeft);
}
model {
  choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  vector[numTrials] log_lik;

  zTest = bias + sens*(testContrastRight - testContrastLeft);
  pRTest = exp(zTest)./(1+exp(zTest));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n] );
  }
}
