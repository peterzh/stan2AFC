#include standard_data_blocks_2AFC.stan
parameters {
    real bias; // bias parameter
    real sens; // sensitivity
}
model {
  //likelihood
  vector[numTrials] z;
  z = bias + sens*(contrastRight - contrastLeft);
  choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  zTest = bias + sens*(testContrastRight - testContrastLeft);
  pRTest = exp(zTest)./(1+exp(zTest));
}
