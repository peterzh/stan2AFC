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
}
transformed data {
	int<lower=0,upper=1> choiceR[numTrials]; // 0=Left, 1=Right
	for (n in 1:numTrials){
		choiceR[n] = choice[n] - 1; 
	}
}
parameters {
	real bias; // grand bias parameter (avg over all sessions)
	real sens; // grand sensitivity (avg over all sessions)

	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions

}
transformed parameters {
	vector[numTrials] z;
	z = bias + bias_delta_perSession[sessionID] + sens*(contrastRight - contrastLeft);
}
model {
	bias_delta_perSession ~ normal(0, bias_delta_perSession_sd); //define hyperprior on bias deltas
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
	real zTest[numTestContrasts,numSessions];
	real pTest[numTestContrasts,numSessions];
	real zTestGrandAverage[numTestContrasts];
	real pTestGrandAverage[numTestContrasts];
	real log_lik[numTrials];

	for (c in 1:numTestContrasts) {
		for (sess in 1:numSessions)
		{
			zTest[c,sess] = bias + bias_delta_perSession[sess] + sens*(testContrastRight[c] - testContrastLeft[c]);
			pTest[c,sess] = exp(zTest[c,sess])./(1+exp(zTest[c,sess]));
		}
		zTestGrandAverage[c] = bias + sens*(testContrastRight[c] - testContrastLeft[c]);
		pTestGrandAverage[c] = exp(zTestGrandAverage[c])./(1+exp(zTestGrandAverage[c]));
	}
	
	for (n in 1:numTrials){
		log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n]);
	}
}
