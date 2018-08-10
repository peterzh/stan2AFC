data {
	int<lower=1> numTrials;
	vector[numTrials] x;
	vector[numTrials] y;
}
transformed data {
	vector[numTrials] mu;
	for (i in 1:numTrials) {
		mu[i] = 0;
	}
}
parameters {
	//vector[numTrials] z;
	
	real<lower=0> eta_sq;
	real<lower=0> inv_rho_sq;
	real<lower=0> sigma_sq;
}
transformed parameters {
	//given fixed parameters, get the covariance 
	matrix[numTrials,numTrials] L; //cholesky factor
	cov_matrix[numTrials] Sigma;
	real<lower=0> rho_sq;
	rho_sq = inv(inv_rho_sq);
	
	for (i in 1:numTrials)
	{
		for (j in 1:numTrials)
		{
			Sigma[i,j] = eta_sq * exp(-rho_sq * (x[i]-x[j])^2 );
			if (i == j) 
			{
				Sigma[i,j] += sigma_sq;
			}
		}
	}

	L = cholesky_decompose(Sigma);
}
model {
	//hyperpriors
	eta_sq ~ cauchy(0, 5);
	inv_rho_sq ~ cauchy(0, 5);
	sigma_sq ~ cauchy(0, 5);
	y ~ multi_normal_cholesky(mu,L);
}
generated quantities {
  vector[numTrials] f = multi_normal_cholesky_rng(mu, L);
  vector[numTrials] y_hat;
  for (n in 1:numTrials)
    y_hat[n] = normal_rng(f[n], sqrt(sigma_sq));
}