data {
	int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=1,upper=3> choice[numTrials]; // 1=Left, 2=Right, 3=NoGo

	int<lower=0> numTestContrasts; //Number of query contrast points
	vector<lower=0,upper=1>[numTestContrasts] testContrastLeft;
	vector<lower=0,upper=1>[numTestContrasts] testContrastRight;
}
parameters {
    real bias_left; // bias parameter
	real bias_right; // bias parameter
    real<lower=0> sens_left; // sensitivity to left contrasts
	real<lower=0> sens_right; // sensitivity to right contrasts
	real<lower=0,upper=2> n_exp; //n exponent
}
transformed parameters {
  vector[3] z[numTrials];  // log odds LvNG
  
  for (n in 1:numTrials)
  {
    z[n][1] = bias_left + sens_left*(contrastLeft[n]^n_exp); 
	z[n][2] = bias_right + sens_right*(contrastRight[n]^n_exp); 
	z[n][3] = 0;
  }
}
model {
	//priors, help with MNL models
	bias_left ~ normal(0, 5);
	bias_right ~ normal(0, 5);
	sens_left ~ normal(0, 20);
	sens_right ~ normal(0, 20);
	
	for (n in 1:numTrials) {
		choice[n] ~ categorical_logit( z[n] );
	}
}
generated quantities {
  vector[3] zTest[numTestContrasts];
  vector[3] pTest[numTestContrasts];
  vector[numTrials] log_lik;

  for (c in 1:numTestContrasts)
	{
		zTest[c][1] = bias_left + sens_left*(testContrastLeft[c]^n_exp); 
		zTest[c][2] = bias_right + sens_right*(testContrastRight[c]^n_exp); 
		zTest[c][3] = 0;
		pTest[c] = softmax( zTest[c] );
	}

  for (n in 1:numTrials){
    log_lik[n] = categorical_logit_lpmf(choice[n] | z[n] );
  }
}

