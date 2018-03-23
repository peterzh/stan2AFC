data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1> numSubjects;
	int<lower=1, upper=numSessions> sessionID[numTrials];
	int<lower=1, upper=numSubjects> subjID[numTrials];
	
	vector[numTrials] contrast_left;
    vector[numTrials] contrast_right;
    int<lower=0,upper=1> choiceR[numTrials];
}
parameters {
    real bias; // grand bias parameter (avg over all subjects and all sessions)
	real sens; // grand sens parameter
	
	vector[numSubjects] bias_delta_perSubject; // bias adjustment parameter for each subject
	real<lower=0> bias_delta_perSubject_sd; // SD of bias adjustment values subjects
	
	vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions
	
	vector[numSubjects] sens_delta_perSubject; // bias adjustment parameter for each subject
	real<lower=0> sens_delta_perSubject_sd; // SD of bias adjustment values subjects
	
	vector[numSessions] sens_delta_perSession; // bias adjustment parameter for each session
	real<lower=0> sens_delta_perSession_sd; // SD of bias adjustment values across sessions
}
model {
	real b;
	real s;
	real z;
	
    bias_delta_perSubject ~ normal(0, bias_delta_perSubject_sd);
    bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	sens_delta_perSubject ~ normal(0, sens_delta_perSubject_sd);
    sens_delta_perSession ~ normal(0, sens_delta_perSession_sd);
	
    for (n in 1:numTrials)
    {
		b = bias + bias_delta_perSubject[subjID[n]] + bias_delta_perSession[sessionID[n]];
		s = sens + sens_delta_perSubject[subjID[n]] + sens_delta_perSession[sessionID[n]];
        z = b + s*(contrast_right[n] - contrast_left[n]);
        choiceR[n] ~ bernoulli_logit(z);
    }
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens.*(contrast_right - contrast_left)