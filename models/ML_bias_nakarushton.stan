data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1> numSubjects;
	
	int<lower=1, upper=numSessions> sessionID[numTrials];
	int<lower=1, upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
    vector<lower=0,upper=1>[numTrials] contrastRight;
    int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
}
parameters {
    real biasL; // bias parameter
	real biasR; // bias parameter
    real<lower=0> sensL; // sensitivity to left contrasts
	real<lower=0> sensR; // sensitivity to right contrasts
	real<lower=0> n_exp; // n_exp
	real<lower=0,upper=1> c_50; // c50
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
	
	for (n in 1:numTrials) 
	{
		CL = (contrastLeft[n]^n_exp) / ( (contrastLeft[n]^n_exp) + (c_50^n_exp) );
		CR = (contrastRight[n]^n_exp) / ( (contrastRight[n]^n_exp) + (c_50^n_exp) );
		
		z[1] = biasL + sensL * CL;
		z[2] = biasR + sensR * CR;
		p_hat = softmax(z);
		choice[n] ~ categorical( p_hat );
	}
	
	
}
