data {
	int<lower=1> numTrials;
	vector[numTrials] x;
	
	real eta;
	real rho;
	real sigma;
}
transformed data {
	//given fixed parameters, get the covariance 
	vector[numTrials] mu;
	matrix[numTrials,numTrials] L;
	cov_matrix[numTrials] Sigma;
	
	for (i in 1:numTrials)
	{
		for (j in 1:numTrials)
		{
			Sigma[i,j] = eta^2 * exp(-(rho^2) * (x[i]-x[j])^2 );
			if (i == j) 
			{
				Sigma[i,j] += sigma^2;
			}
		}
	}

	L = cholesky_decompose(Sigma);
	
	for (i in 1:numTrials) 
	{
		mu[i] = 0;
	}
	
}
parameters {
	vector[numTrials] z;
}
model {
	z ~ normal(0,1);
}
generated quantities {
	vector[numTrials] y;
	y = mu + L * z;
}