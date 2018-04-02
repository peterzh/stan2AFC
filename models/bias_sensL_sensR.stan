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
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
}
model {
  //likelihood
  choiceR ~ bernoulli_logit( bias + sens_left*contrastLeft + sens_right*contrastRight );
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens_left.*contrastLeft + p.sens_right.*contrastRight
