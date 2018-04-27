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
	int<lower=1,upper=numSubjects> subjectID_perSession[numSessions];
	for (n in 1:numTrials){
		choiceR[n] = choice[n] - 1; 
	}
	
	//identify the subject for each session. Poorly coded here! I shouldn't need to go through all trials for this..
	for (n in 1:numTrials){
		subjectID_perSession[sessionID[n]] = subjectID[n];
	}
}
parameters {
	real bias; // grand bias parameter (avg over all subjects and all sessions)
	real sens; // grand sens parameter

	vector[numSubjects] bias_delta_perSubject; // bias adjustment parameter for each subject
	real<lower=0> bias_delta_perSubject_sd; // SD of bias adjustment values subjects

	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions

	vector[numSubjects] sens_delta_perSubject; // sens adjustment parameter for each subject
	real<lower=0> sens_delta_perSubject_sd; // SD of sens adjustment values subjects
}
transformed parameters {
	vector[numTrials] z;
	z = bias + bias_delta_perSession[sessionID] + bias_delta_perSubject[subjectID] + rows_dot_product(sens + sens_delta_perSubject[subjectID], contrastRight - contrastLeft);
}
model {
	bias_delta_perSubject ~ normal(0, bias_delta_perSubject_sd);
	bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	sens_delta_perSubject ~ normal(0, sens_delta_perSubject_sd);
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
	real zTest[numTestContrasts,numSessions];
	real pTest[numTestContrasts,numSessions];
	real zTestSubjectAverage[numTestContrasts,numSubjects];
	real pTestSubjectAverage[numTestContrasts,numSubjects];
	real zTestGrandAverage[numTestContrasts];
	real pTestGrandAverage[numTestContrasts];
	real log_lik[numTrials];

	for (c in 1:numTestContrasts) {
		for (sess in 1:numSessions)
		{
			zTest[c,sess] = bias + bias_delta_perSession[sess] + bias_delta_perSubject[subjectID_perSession[sess]] + (sens + sens_delta_perSubject[subjectID_perSession[sess]])*(testContrastRight[c] - testContrastLeft[c]);
			pTest[c,sess] = exp(zTest[c,sess])./(1+exp(zTest[c,sess]));
		}
		
		for (subj in 1:numSubjects) {
			zTestSubjectAverage[c,subj] = bias + bias_delta_perSubject[subj] + (sens + sens_delta_perSubject[subj])*(testContrastRight[c] - testContrastLeft[c]);
			pTestSubjectAverage[c,subj] = exp(zTestSubjectAverage[c,subj])./(1+exp(zTestSubjectAverage[c,subj]));
		}
		zTestGrandAverage[c] = bias + sens*(testContrastRight[c] - testContrastLeft[c]);
		pTestGrandAverage[c] = exp(zTestGrandAverage[c])./(1+exp(zTestGrandAverage[c]));
	}

	for (n in 1:numTrials){
		log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n]);
	}
}
