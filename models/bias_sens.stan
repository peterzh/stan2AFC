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
    real bias; // bias parameter
    real sens; // sensitivity
}
model {
  //likelihood
  vector[numTrials] z;
  z = bias + sens*(contrastRight - contrastLeft);
  choiceR ~ bernoulli_logit( z );
}
generated quantities {
  //generate predicted probabilities for many contrast values
  vector[1000] cDiffQuery;
  vector[1000] zQuery;
  vector[1000] pRQuery;
  cDiffQuery = (cumulative_sum(rep_vector(1,1000))-500)/500;
  zQuery = bias + sens*cDiffQuery;
  pRQuery = exp(zQuery)./(1+exp(zQuery));
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)
