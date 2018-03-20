
%Get 2AFC detection behaviour
sessions = fetch(d.Session & 'project_id LIKE "%2AFC%"' & 'num_trials>200');

for s = 1:length(sessions)
    trials = d.Trial & sessions(s);
    [contrast_left, contrast_right] = fetchn( trials * d.TrialStim, 'contrast_left', 'contrast_right');
    [response] = fetchn( trials * d.TrialResponse, 'response');
    
    %Change response coding 0=Left, 1=Right
    response(response==1)=0;
    response(response==2)=1;
    
    data = struct('N',length(response),'contrast_left',contrast_left,...
                    'contrast_right',contrast_right,'choiceR',response);
    
    b = behavModel(data,'bias_sensL_sensR_exponentL_exponentR');
    b.plotPosterior;
    b.plotPsych;
    drawnow;
end

