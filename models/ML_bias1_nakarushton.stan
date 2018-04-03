#include standard_data_blocks_2AUC.stan
parameters {
    real biasL; // bias parameter
	real biasR; // bias parameter
    real<lower=0> sensL; // sensitivity to left contrasts
	real<lower=0> sensR; // sensitivity to right contrasts
	real<lower=0> n_exp; // n_exp
	real<lower=0,upper=1> c_50; // c50

	vector[numSessions] biasL_delta_perSession; // bias adjustment parameter for each session
	vector[numSessions] biasR_delta_perSession; // bias adjustment parameter for each session
	//real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions
}
model {
	real CL;
    real CR;
	vector[3] z;
	vector[3] p_hat;
	z[3] = 0;

	//priors
	biasL ~ normal(0, 5);
	biasR ~ normal(0, 5);
	sensL ~ normal(0, 20);
	sensR ~ normal(0, 20);
	n_exp ~ normal(0, 20);
    biasL_delta_perSession ~ normal(0, 2);
    biasR_delta_perSession ~ normal(0, 2);

	for (n in 1:numTrials)
	{
		CL = (contrastLeft[n]^n_exp) / ( (contrastLeft[n]^n_exp) + (c_50^n_exp) );
		CR = (contrastRight[n]^n_exp) / ( (contrastRight[n]^n_exp) + (c_50^n_exp) );

		z[1] = biasL + biasL_delta_perSession[sessionID[n]] + sensL * CL;
		z[2] = biasR + biasR_delta_perSession[sessionID[n]] + sensR * CR;
		p_hat = softmax(z);
		choice[n] ~ categorical( p_hat );
	}


}
