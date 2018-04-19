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
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp; //n exponent
    
	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
  	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions
	
	vector[numSessions] sens_left_delta_perSession; // bias adjustment parameter for each session
  	real<lower=0> sens_left_delta_perSession_sd; // SD of bias adjustment values across sessions
	
	vector[numSessions] sens_right_delta_perSession; // bias adjustment parameter for each session
  	real<lower=0> sens_right_delta_perSession_sd; // SD of bias adjustment values across sessions
	
	vector[numSessions] n_delta_perSession; // bias adjustment parameter for each session
  	real<lower=0> n_delta_perSession_sd; // SD of bias adjustment values across sessions
}
transformed parameters {
  vector[numTrials] z;
  for (n in 1:numTrials)
	{
    z[n] = bias + bias_delta_perSession[sessionID[n]] + (sens_left + sens_left_delta_perSession[sessionID[n]])*contrastLeft[n]^(n_exp + n_delta_perSession[sessionID[n]]) + (sens_right + sens_right_delta_perSession[sessionID[n]])*contrastRight[n]^(n_exp + n_delta_perSession[sessionID[n]]) ;
	}
}
model {
  bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
  sens_left_delta_perSession ~ normal(0, sens_left_delta_perSession_sd);
  sens_right_delta_perSession ~ normal(0, sens_right_delta_perSession_sd);
  n_delta_perSession ~ normal(0, n_delta_perSession_sd);
  choice ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] testContrastLeftTransformed;
  vector[numTestContrasts] testContrastRightTransformed;
  vector[numTestContrasts] zTest[numSessions];
  vector[numTestContrasts] pRTest[numSessions];
  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;
  vector[numTrials] log_lik;

  for (sess in 1:numSessions)
  {
    for (c in 1:numTestContrasts)
	{
		testContrastLeftTransformed[c] = testContrastLeft[c]^(n_exp + n_delta_perSession[sess]);
		testContrastRightTransformed[c] = testContrastRight[c]^(n_exp + n_delta_perSession[sess]);
	}
    zTest[sess] = bias + bias_delta_perSession[sess] + (sens_left + sens_left_delta_perSession[sess])*testContrastLeftTransformed + (sens_right + sens_right_delta_perSession[sess])*testContrastRightTransformed;
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }

  for (c in 1:numTestContrasts)
	{
		testContrastLeftTransformed[c] = testContrastLeft[c]^n_exp;
		testContrastRightTransformed[c] = testContrastRight[c]^n_exp;
	}
	
  zTestGrandAverage = bias + sens_left*testContrastLeftTransformed + sens_right*testContrastRightTransformed;
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choice[n] | z[n]);
  }
}
