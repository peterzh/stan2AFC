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
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp; //n exponent
}
transformed parameters {
  vector[numTrials] z;
  for (n in 1:numTrials)
  {
    z[n] = bias + sens_left*(contrastLeft[n]^n_exp) + sens_right*(contrastRight[n]^n_exp) ;
  }
}
model {
	choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  vector[numTrials] log_lik;

  for (c in 1:numTestContrasts)
	{
		zTest[c] = bias + sens_left*(testContrastLeft[c]^n_exp) + sens_right*(testContrastRight[c]^n_exp) ;
	}
  pRTest = exp(zTest)./(1+exp(zTest));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n] );
  }
}
