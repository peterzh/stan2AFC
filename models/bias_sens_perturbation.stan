data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=2> choice[numTrials]; // 1=Left, 2=Right
	int<lower=0> numTestContrasts; //Number of query contrast points
	vector<lower=0,upper=1>[numTestContrasts] testContrastLeft;
	vector<lower=0,upper=1>[numTestContrasts] testContrastRight;
	int<lower=0,upper=1> perturbation[numTrials];
}
transformed data {
	int<lower=0,upper=1> choiceR[numTrials]; // 0=Left, 1=Right
	for (n in 1:numTrials){
		choiceR[n] = choice[n] - 1; 
	}
}
parameters {
	real bias; // bias parameter
	real sens; // sensitivity
	real bias_perturbation;
	real sens_perturbation;
}
transformed parameters {
	vector[numTrials] z;
	for (n in 1:numTrials)
	{
		if (perturbation[n] == 1) {
			z[n] = (bias + bias_perturbation) + (sens + sens_perturbation)*(contrastRight[n] - contrastLeft[n]);
		} else {
			z[n] = bias + sens*(contrastRight[n] - contrastLeft[n]);
		}
	}
}
model {
	bias_perturbation ~ normal(0, 10); //prior on the perturbations
	sens_perturbation ~ normal(0, 10);
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
	real zTest[numTestContrasts];
	real pTest[numTestContrasts];
	real zTest_perturbation[numTestContrasts];
	real pTest_perturbation[numTestContrasts];
	real log_lik[numTrials];
	
	for (c in 1:numTestContrasts) {
		zTest[c] = bias + sens*(testContrastRight[c] - testContrastLeft[c]);
		pTest[c] = exp(zTest[c])/(1 + exp(zTest[c]));
		
		zTest_perturbation[c] = (bias + bias_perturbation) + (sens + sens_perturbation)*(testContrastRight[c] - testContrastLeft[c]);
		pTest_perturbation[c] = exp(zTest_perturbation[c])/(1 + exp(zTest_perturbation[c]));
	}

	for (n in 1:numTrials){
		log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n] );
	}
}