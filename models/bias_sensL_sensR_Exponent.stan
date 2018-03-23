data {
    int<lower=0> numTrials;
    vector[numTrials] contrast_left;
    vector[numTrials] contrast_right;
    int<lower=0,upper=1> choiceR[numTrials];
}
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
		z[n] = bias + sens_left*(contrast_left[n]^n_exp) + sens_right*(contrast_right[n]^n_exp) ;
	}
	
	choiceR ~ bernoulli_logit( z );
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens_left.*(contrast_left.^p.n_exp) + p.sens_right.*(contrast_right.^p.n_exp)
