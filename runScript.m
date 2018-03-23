
%Get 2AFC detection sessions from one mouse
sessions = d.Session & 'stimulus_type="Detection"' & 'choice_type="2AFC"' & 'performance>0.7';
sessions = proj(sessions,'CONCAT(session_date,"_",session_num,"_",mouse_name)->expRef');
trials = (d.Trial * d.TrialStim * d.TrialResponse) * sessions;
[mouse_name,expRef, contrast_left, contrast_right, response] = fetchn( trials, 'mouse_name', 'expRef', 'contrast_left', 'contrast_right', 'response');
response(response==1)=0;
response(response==2)=1;
    
[~,~,sessionID] = unique(expRef);
[~,~,subjID] = unique(mouse_name);

data = struct('numTrials',length(response),...
              'numSubjects', max(subjID),...
              'numSessions', max(sessionID),...
              'contrast_left', contrast_left,...
              'contrast_right', contrast_right,...
              'choiceR', response,...
              'sessionID', sessionID,...
              'subjID', subjID);

b = behavModel(data,'hier_bias_sens');
b.getPosterior;
b.plotPosteriorMarginals;
% b.plotPosterior;
% b.plotPsych;
drawnow;