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
  choiceR ~ bernoulli_logit( bias + sens*(contrastRight - contrastLeft) );
}

//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)
