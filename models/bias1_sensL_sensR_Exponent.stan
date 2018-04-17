#include standard_data_blocks_2AFC.stan
parameters {
    real bias;
    real<upper=0> sens_left;
    real<lower=0> sens_right;
    real<lower=0,upper=2> n_exp; //n exponent
    vector[numSessions] bias_delta_perSession; // bias adjustment parameter for each session
  	real<lower=0> bias_delta_perSession_sd; // SD of bias adjustment values across sessions
}
transformed parameters {
  vector[numTrials] z;
  for (n in 1:numTrials)
	{
    z[n] = bias + bias_delta_perSession[sessionID[n]] + sens_left*contrastLeft[n]^n_exp + sens_right*contrastRight[n]^n_exp ;
	}
}
model {
  bias_delta_perSession ~ normal(0, bias_delta_perSession_sd);
  choiceR ~ bernoulli_logit( z );
}
generated quantities {
  vector[numTestContrasts] testContrastLeftTransformed;
  vector[numTestContrasts] testContrastRightTransformed;
  vector[numTestContrasts] zTest[numSessions];
  vector[numTestContrasts] pRTest[numSessions];
  vector[numTestContrasts] zTestGrandAverage;
  vector[numTestContrasts] pRTestGrandAverage;
  vector[numTrials] log_lik;

  for (c in 1:numTestContrasts)
  {
    testContrastLeftTransformed[c] = testContrastLeft[c]^n_exp;
    testContrastRightTransformed[c] = testContrastRight[c]^n_exp;
  }

  for (sess in 1:numSessions)
  {
    zTest[sess] = bias + bias_delta_perSession[sess] + sens_left*testContrastLeftTransformed + sens_right*testContrastRightTransformed;
    pRTest[sess] = exp(zTest[sess])./(1+exp(zTest[sess]));
  }

  zTestGrandAverage = bias + sens_left*testContrastLeftTransformed + sens_right*testContrastRightTransformed;
  pRTestGrandAverage = exp(zTestGrandAverage)./(1+exp(zTestGrandAverage));

  for (n in 1:numTrials){
    log_lik[n] = bernoulli_logit_lpmf(choiceR[n] | z[n]);
  }
}
