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
	
	vector[numSessions] sens_delta_perSession; // sens adjustment parameter for each session
	real<lower=0> sens_delta_perSession_sd; // SD of sens adjustment values across sessions

}
model {
	vector[numTrials] B;
	vector[numTrials] S;
	vector[numTrials] C;

	bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
    sens_delta_perSession ~ normal(0, sens_delta_perSession_sd);
	
	B = bias + bias_delta_perSession[sessionID]; // grand mean bias + per-session bias deviations
	S = sens + sens_delta_perSession[sessionID]; // grand mean sens + per-session sens deviations
	C = contrastRight - contrastLeft; // contrast difference on each trial
	choiceR ~ bernoulli_logit( B + rows_dot_product(S,C) );
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)