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
generated quantities {
  //vector[numTestContrasts] zTest[numSessions];
  //vector[numTestContrasts] pRTest[numSessions];

  vector[numTestContrasts] zTestSubjectAverage[numSubjects];
  vector[numTestContrasts] pRTestSubjectAverage[numSubjects];

  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;

  /*
  for (sess in 1:numSessions)
  {
    subjIDthisSession =
    zTest[sess] = bias + bias_delta_perSubject[subjIDthisSession] + bias_delta_perSession[sess] + (sens + sens_delta_perSubject[subjIDthisSession] + sens_delta_perSession[sess])*(testContrastRight - testContrastLeft);
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }
  */

  for (subj in 1:numSubjects)
  {
    zTestSubjectAverage[subj] = bias + bias_delta_perSubject[subj] + (sens + sens_delta_perSubject[subj])*(testContrastRight - testContrastLeft);
    pRTestSubjectAverage[subj] = exp(zTestSubjectAverage[subj])./(1+exp(zTestSubjectAverage[subj]));
  }

  //grand average
  zTestGrandAverage = bias + sens*(testContrastRight - testContrastLeft);
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)
