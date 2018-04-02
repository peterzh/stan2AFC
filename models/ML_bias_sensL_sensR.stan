data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
    vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo
}
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
