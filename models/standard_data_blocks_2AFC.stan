data {
  int<lower=1> numTrials;
	int<lower=1> numSessions;
	int<lower=1> numSubjects;
	int<lower=1,upper=numSessions> sessionID[numTrials];
	int<lower=1,upper=numSubjects> subjectID[numTrials];
	vector<lower=0,upper=1>[numTrials] contrastLeft;
	vector<lower=0,upper=1>[numTrials] contrastRight;
	int<lower=0,upper=1> choice[numTrials]; // 0=Left, 1=Right

  int<lower=0> numTestContrasts; //Number of query contrast points
  vector<lower=0,upper=1>[numTestContrasts] testContrastLeft;
  vector<lower=0,upper=1>[numTestContrasts] testContrastRight;
}
