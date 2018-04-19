data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=0,upper=1> choice[numTrials]; // 0=Left, 1=Right
	int<lower=0> numTestContrasts; //Number of query contrast points
	vector<lower=0,upper=1>[numTestContrasts] testContrastLeft;
	vector<lower=0,upper=1>[numTestContrasts] testContrastRight;
}
parameters {
  real bias; // grand bias parameter (avg over all subjects and all sessions)
	real sens; // grand sens parameter

	vector[numSubjects] bias_delta_perSubject; // bias adjustment parameter for each subject
	real<lower=0> bias_delta_perSubject_sd; // SD of bias adjustment values subjects

	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions

	vector[numSubjects] sens_delta_perSubject; // bias adjustment parameter for each subject
	real<lower=0> sens_delta_perSubject_sd; // SD of bias adjustment values subjects

	vector[numSessions] sens_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> sens_delta_perSession_sd; // SD of bias adjustment values across sessions
}
model {

	vector[numTrials] B;
	vector[numTrials] S;
	vector[numTrials] C;

	bias_delta_perSubject ~ normal(0, bias_delta_perSubject_sd);
  bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	sens_delta_perSubject ~ normal(0, sens_delta_perSubject_sd);
  sens_delta_perSession ~ normal(0, sens_delta_perSession_sd);

	B = bias + bias_delta_perSubject[subjectID] + bias_delta_perSession[sessionID];
	S = sens + sens_delta_perSubject[subjectID] + sens_delta_perSession[sessionID];
	C = contrastRight - contrastLeft; // contrast difference on each trial
	choice ~ bernoulli_logit( B + rows_dot_product(S,C) );

}
generated quantities {
  //vector[numTestContrasts] zTest[numSessions];
  //vector[numTestContrasts] pRTest[numSessions];

  vector[numTestContrasts] zTestSubjectAverage[numSubjects];
  vector[numTestContrasts] pRTestSubjectAverage[numSubjects];

  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;

  /*
  for (sess in 1:numSessions)
  {
    subjIDthisSession =
    zTest[sess] = bias + bias_delta_perSubject[subjIDthisSession] + bias_delta_perSession[sess] + (sens + sens_delta_perSubject[subjIDthisSession] + sens_delta_perSession[sess])*(testContrastRight - testContrastLeft);
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }
  */

  for (subj in 1:numSubjects)
  {
    zTestSubjectAverage[subj] = bias + bias_delta_perSubject[subj] + (sens + sens_delta_perSubject[subj])*(testContrastRight - testContrastLeft);
    pRTestSubjectAverage[subj] = exp(zTestSubjectAverage[subj])./(1+exp(zTestSubjectAverage[subj]));
  }

  //grand average
  zTestGrandAverage = bias + sens*(testContrastRight - testContrastLeft);
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));
}
