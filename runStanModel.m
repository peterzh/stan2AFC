
%% get 2AFC data
sessions = d.Session & 'choice_type="2AFC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
response = response-1;
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('contrastLeft', contrast_left,...
              'contrastRight', contrast_right,...
              'choiceR', response,...
              'sessionID', sessionID,...
              'subjectID', subjID);
          
data2 = getrow(data, data.subjectID==12);
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

%% Load model
b = behavModel(data,'bias_sens');
b.getPosterior;
b.plotPosteriorMarginals;