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
    real bias; // grand bias parameter (avg over all sessions)
    real sens; // grand sensitivity (avg over all sessions)
	
	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions

}
model {
	vector[numTrials] B;
	vector[numTrials] C;
	
    bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	
	B = bias + bias_delta_perSession[sessionID]; // grand mean bias + per-session bias deviations
	C = contrastRight - contrastLeft; // contrast difference on each trial
	choiceR ~ bernoulli_logit( B + sens*C );
	
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens.*(contrastRight - contrastLeft)