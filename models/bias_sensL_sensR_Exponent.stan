#include standard_data_blocks_2AFC.stan
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp; //n exponent
}
model {
	vector[numTrials] z;
	for (n in 1:numTrials)
	{
		z[n] = bias + sens_left*(contrastLeft[n]^n_exp) + sens_right*(contrastRight[n]^n_exp) ;
	}

	choiceR ~ bernoulli_logit( z );
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens_left.*(contrastLeft.^p.n_exp) + p.sens_right.*(contrastRight.^p.n_exp)
