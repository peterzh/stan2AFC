#include standard_data_blocks_2AFC.stan
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

	vector[numTrials] B;
	vector[numTrials] S;
	vector[numTrials] C;

	bias_delta_perSubject ~ normal(0, bias_delta_perSubject_sd);
    bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
	sens_delta_perSubject ~ normal(0, sens_delta_perSubject_sd);
    sens_delta_perSession ~ normal(0, sens_delta_perSession_sd);

	B = bias + bias_delta_perSubject[subjectID] + bias_delta_perSession[sessionID];
	S = sens + sens_delta_perSubject[subjectID] + sens_delta_perSession[sessionID];
	C = contrastRight - contrastLeft; // contrast difference on each trial
	choiceR ~ bernoulli_logit( B + rows_dot_product(S,C) );

}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)
