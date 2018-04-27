
%% get 2AFC data
sessions = d.Session & 'choice_type="2AFC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('contrastLeft', contrast_left,...
              'contrastRight', contrast_right,...
              'choice', response,...
              'sessionID', sessionID,...
              'subjectID', subjID);
          
data2 = getrow(data, data.subjectID==1);
[~,~,data2.sessionID]=unique(data2.sessionID);
[~,~,data2.subjectID]=unique(data2.subjectID);
          
          
%% get 2AUC data
sessions = d.Session & 'choice_type="2AUC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('contrastLeft', contrast_left,...
              'contrastRight', contrast_right,...
              'choice', response,...
              'sessionID', sessionID,...
              'subjectID', subjID);
data2 = getrow(data, data.subjectID==1);
[~,~,data2.sessionID]=unique(data2.sessionID);
[~,~,data2.subjectID]=unique(data2.subjectID);


%% Simulate data
N=10000;
bias=2;
sens=2.5;
cDiff = repmat([-1;-0.25;0;0.25;1],N/5,1);
z = bias+sens*cDiff;
pR = exp(z)./(1+exp(z));


data3 = struct('contrastLeft', (cDiff<0).*(-cDiff),...
              'contrastRight', (cDiff>0).*cDiff,...
              'choiceR', binornd(1,pR),...
              'sessionID', ones(N,1),...
              'subjectID', ones(N,1));
          
%% Load model
b = behavModel(data2,'bias1_sensL_sensR_Exponent_perturbation');
b.getPosterior;
b.plotPosteriorMarginals;
b.plotPsych;