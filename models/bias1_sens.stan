#include standard_data_blocks_2AFC.stan
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
generated quantities { //todo: modify to generate per-session fits
  vector[numTestContrasts] zTest[numSessions];
  vector[numTestContrasts] pRTest[numSessions];
  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;

  for (sess in 1:numSessions)
  {
    zTest[sess] = bias + bias_delta_perSession[sess] + sens*(testContrastRight - testContrastLeft);;
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }

  zTestGrandAverage = bias + sens*(testContrastRight - testContrastLeft);
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));
}
//Z function(s):
//@(p,contrast_left,contrast_right) p.bias + p.sens.*(contrastRight - contrastLeft)
