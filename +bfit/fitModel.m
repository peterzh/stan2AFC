function fit = fitModel(modelName,data)
%modelname: string
%data: struct

%Choose where to save the fit
root_dir = fullfile( fileparts(mfilename('fullpath')),'..' , 'fits');
[file,pathname] = uiputfile(fullfile(root_dir,'*.mat'));
file = fullfile(pathname,file);

%Get model object
stanModelObj = bfit.getCompiledModel(modelName);

%Prepare the dataset to be compatible with stan
%Convert data structure to stan data format
data.numTrials = length(data.contrastLeft);

if isfield(data,'sessionID')
    data.numSessions = max(data.sessionID);
    assert(min(data.sessionID)==1,'SessionID does not start at 1');
end

if isfield(data,'subjectID')
    data.numSubjects = max(data.subjectID);
    assert(min(data.subjectID)==1,'subjectID does not start at 1');
    ss = unique([data.subjectID data.sessionID],'rows');
    data.subjID_session = ss(:,1);
end

if isfield(data,'perturbation')
    data.numPerturbations = max(data.perturbation);
    assert(min(data.perturbation)==0,'Perturbation variable does not contain zeros');
end

if isfield(data,'perturbationTime')
    assert( all(diff(data.perturbationTime(data.perturbation==1))>=0), 'Trials must be ordered by perturbation time'); 
end

%Get posterior
%Fit model on data
fitObj = stan('fit',stanModelObj,'method','sample','data',data,'iter',1000,'chains',4,'verbose',true);
fitObj.block;

%Get all parameter values
p = fitObj.extract('permuted',false);

%Convert to scalar struct
fields = fieldnames(p);
p1 = p(1);
for i = 1:length(fields)
    p1.(fields{i}) = cat(1,p.(fields{i}));
end

p1 = rmfield(p1, fields(contains(fields,{'__'})));
p1 = rmfield(p1, fields(contains(fields,{'logOdds'})));
posterior = p1;

warning('TODO: Assess convergence');
% waic=mstan.waic(p1.log_lik);


%Save fit
fit = struct;
fit.modelName = modelName;
fit.posterior = posterior;
fit.data = data;

save(file,'-struct','fit','-v7.3');
end