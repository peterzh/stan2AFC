#include standard_data_blocks_2AFC.stan
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
generated quantities {
  vector[numTestContrasts] zTest[numSessions];
  vector[numTestContrasts] pRTest[numSessions];
  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;

  for (sess in 1:numSessions)
  {
    zTest[sess] = bias + bias_delta_perSession[sess] + (sens + sens_delta_perSession[sess])*(testContrastRight - testContrastLeft);
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }

  zTestGrandAverage = bias + sens*(testContrastRight - testContrastLeft);
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));
}
//Z function(s):
//@(p,contrastLeft,contrastRight) p.bias + p.sens.*(contrastRight - contrastLeft)
