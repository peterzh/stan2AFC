#include standard_data_blocks_2AFC.stan
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp_left; //n exponent
	  real<lower=0,upper=2> n_exp_right; //n exponent
}
model {
	vector[numTrials] z;

	for (n in 1:numTrials)
	{
		z[n] = bias + sens_left*(contrastLeft[n]^n_exp_left) + sens_right*(contrastRight[n]^n_exp_right) ;
	}

	choiceR ~ bernoulli_logit( z );

}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  for (c in 1:numTestContrasts)
	{
		zTest[c] = bias + sens_left*(testContrastLeft[c]^n_exp_left) + sens_right*(testContrastRight[c]^n_exp_right) ;
	}
  pRTest = exp(zTest)./(1+exp(zTest));
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens_left.*(contrastLeft.^p.n_exp_left) + p.sens_right.*(contrastRight.^p.n_exp_right)
