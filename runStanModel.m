
%% get 2AFC data
sessions = d.Session & 'choice_type="2AFC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
response = response-1;
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('numTrials',length(response),...
              'numSubjects', max(subjID),...
              'numSessions', max(sessionID),...
              'contrastLeft', contrast_left,...
              'contrastRight', contrast_right,...
              'choiceR', response,...
              'sessionID', sessionID,...
              'subjectID', subjID);
          
          
%% get 2AUC data
sessions = d.Session & 'choice_type="2AUC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('numTrials',length(response),...
              'numSubjects', max(subjID),...
              'numSessions', max(sessionID),...
              'contrastLeft', contrast_left,...
              'contrastRight', contrast_right,...
              'choice', response,...
              'sessionID', sessionID,...
              'subjectID', subjID);
  
%% Add test contrasts
data.numTestContrasts = 1000;
data.testContrastLeft = [linspace(1,0,data.numTestContrasts/2)'; zeros(data.numTestContrasts/2,1)];
data.testContrastRight = [zeros(data.numTestContrasts/2,1); linspace(0,1,data.numTestContrasts/2)'];


%% Load model
b = behavModel(data,'bias_sens');
b.getPosterior;
b.plotPosteriorMarginals;