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

}
model {

    bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	
	choiceR ~ bernoulli_logit( bias + bias_delta_perSession[sessionID] + sens*(contrast_right - contrast_left) );
	
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens.*(contrast_right - contrast_left)