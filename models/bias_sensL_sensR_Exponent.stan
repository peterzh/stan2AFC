#include standard_data_blocks_2AFC.stan
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp; //n exponent
}
transformed parameters {
  vector[numTrials] z;
  for (n in 1:numTrials)
  {
    z[n] = bias + sens_left*(contrastLeft[n]^n_exp) + sens_right*(contrastRight[n]^n_exp) ;
  }
}
model {
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  vector[numTrials] log_lik;

  for (c in 1:numTestContrasts)
	{
		zTest[c] = bias + sens_left*(testContrastLeft[c]^n_exp) + sens_right*(testContrastRight[c]^n_exp) ;
	}
  pRTest = exp(zTest)./(1+exp(zTest));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n] );
  }
}
