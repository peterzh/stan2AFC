data {
  int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
  vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=0,upper=1> choiceR[numTrials]; // 0=Left, 1=Right

  int<lower=0> numTestContrasts; //Number of query contrast points
  vector<lower=0,upper=1>[numTestContrasts] testContrastLeft;
  vector<lower=0,upper=1>[numTestContrasts] testContrastRight;
}
parameters {
  real bias; // grand bias parameter (avg over all sessions)
  real sens; // grand sensitivity (avg over all sessions)

	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions

	vector[numSessions] sens_delta_perSession; // sens adjustment parameter for each session
	real<lower=0> sens_delta_perSession_sd; // SD of sens adjustment values across sessions
}
transformed parameters {
  vector[numTrials] z;
  z = bias + bias_delta_perSession[sessionID] + rows_dot_product(sens + sens_delta_perSession[sessionID], contrastRight - contrastLeft);
}
model {
	bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
  sens_delta_perSession ~ normal(0, sens_delta_perSession_sd);
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest[numSessions];
  vector[numTestContrasts] pRTest[numSessions];
  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;
  vector[numTrials] log_lik;

  for (sess in 1:numSessions)
  {
    zTest[sess] = bias + bias_delta_perSession[sess] + (sens + sens_delta_perSession[sess])*(testContrastRight - testContrastLeft);
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }

  zTestGrandAverage = bias + sens*(testContrastRight - testContrastLeft);
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n]);
  }
}
