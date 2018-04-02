data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
    vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=0,upper=1> choiceR[numTrials]; // 0=Left, 1=Right
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
		z[n] = bias + sens_left*(contrastLeft[n]^n_exp) + sens_right*(contrastRight[n]^n_exp) ;
	}
	
	choiceR ~ bernoulli_logit( z );
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens_left.*(contrastLeft.^p.n_exp) + p.sens_right.*(contrastRight.^p.n_exp)
