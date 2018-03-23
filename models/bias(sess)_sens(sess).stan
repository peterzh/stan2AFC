data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1, upper=numSessions> sessionID[numTrials];
	
	vector[numTrials] contrast_left;
    vector[numTrials] contrast_right;
    int<lower=0,upper=1> choiceR[numTrials];
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
	C = contrast_right - contrast_left; // contrast difference on each trial
	choiceR ~ bernoulli_logit( B + rows_dot_product(S,C) );
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens.*(contrast_right - contrast_left)