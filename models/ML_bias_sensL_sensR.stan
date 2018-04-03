#include standard_data_blocks_2AUC.stan
parameters {
    real biasL; // bias parameter
	real biasR; // bias parameter
    real<lower=0> sensL; // sensitivity to left contrasts
	real<lower=0> sensR; // sensitivity to right contrasts
}
model {

	vector[3] z;
	vector[3] p_hat;

	z[3] = 0;

	//priors
	biasL ~ normal(0, 5);
	biasR ~ normal(0, 5);
	sensL ~ normal(0, 20);
	sensR ~ normal(0, 20);

	for (n in 1:numTrials)
	{
		z[1] = biasL + sensL * contrastLeft[n];
		z[2] = biasR + sensR * contrastRight[n];
		p_hat = softmax(z);
		choice[n] ~ categorical( p_hat );
	}


}
