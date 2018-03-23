// Standard stan model code for 2AFC datasets
data {
    int<lower=1> numTrials;
	int<lower=1> numSessions; 
	int<lower=1> numSubjects;
	
	int<lower=1, upper=numSessions> sessionID[numTrials];
	int<lower=1, upper=numSubjects> subjID[numTrials];
	vector[numTrials] contrast_left;
    vector[numTrials] contrast_right;
    int<lower=0,upper=1> choiceR[numTrials]; // 0 for left choice, 1 for right choice
}