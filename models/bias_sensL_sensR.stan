#include standard_data_blocks_2AFC.stan
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
}
model {
  //likelihood
  choiceR ~ bernoulli_logit( bias + sens_left*contrastLeft + sens_right*contrastRight );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  zTest = bias + sens_left*testContrastLeft + sens_right*testContrastRight;
  pRTest = exp(zTest)./(1+exp(zTest));
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens_left.*contrastLeft + p.sens_right.*contrastRight
