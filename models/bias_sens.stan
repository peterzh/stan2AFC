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
    real bias; // bias parameter
    real sens; // sensitivity
}
transformed parameters {
  vector[numTrials] z;
  z = bias + sens*(contrastRight - contrastLeft);
}
model {
  choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] zTest;
  vector[numTestContrasts] pRTest;
  vector[numTrials] log_lik;

  zTest = bias + sens*(testContrastRight - testContrastLeft);
  pRTest = exp(zTest)./(1+exp(zTest));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n] );
  }
}
