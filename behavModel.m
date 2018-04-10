classdef behavModel < handle
    properties
        modelName;
        data;
    end
    
    properties(Access=public)
        stanModelObj;
        z;
        MAP;
        Posterior;
        PosteriorType;
        data_stan;
    end
    
    methods
        function obj = behavModel(data,modelName)
            obj.modelName = modelName;
            obj.data = data;
            
            %Define directory containing the compiled model
            models_dir = fullfile( fileparts(mfilename('fullpath')) , 'models');
            
            %If compiled model .mat exists, load it
            modelFile = fullfile(models_dir, [modelName '.mat']);
            
            if exist( modelFile, 'file')
                load(modelFile,'sm','z');
                fprintf('Loaded model %s\n', modelName);
                
            else %If it doesn't exist, compile and save it
                fprintf('Compiled model %s not found, compiling now...\n', modelName);
                
                %Get stan model file
                stanFile = strrep(modelFile, '.mat', '.stan');
                
                %Create symbolic link to ensure compiler finds the #include
                %headers
                [~,~]=system(['ln -s ' models_dir '/standard_data_blocks_2AFC.stan'  ' ' mstan.stan_home '/standard_data_blocks_2AFC.stan']);
                [~,~]=system(['ln -s ' models_dir '/standard_data_blocks_2AUC.stan'  ' ' mstan.stan_home '/standard_data_blocks_2AUC.stan']);
                
                %Compile model
                sm = StanModel('file',stanFile,'working_dir',tempdir);
                %                 sm = StanModel('file',stanFile,'working_dir',fileparts(stanFile));
                %                 sm.compile('model','STANCFLAGS =include_paths C:\Users\Peter\Documents\MATLAB\stan2AFC\models');
                sm.compile();
                
                %Read Z function from the stan file as well
                c = textread(stanFile,'%c')';
                zIdx = strfind(c,'//@');
                
                z = {};
                for i = 1:length(zIdx)
                    z{i} = eval( c( (2+zIdx(i)):end ) );
                end
                
                %Save copy of compiled model .mat object
                save(modelFile,'sm','z');
            end
            
            %Store stanModel object in local behavModel object
            obj.stanModelObj = sm;
            obj.z = z;
            
            %Convert data structure to stan data format
            obj.data_stan = obj.data;
            obj.data_stan.numTrials = length(obj.data.choiceR);
            obj.data_stan.numSubjects = max(obj.data.subjectID);
            obj.data_stan.numSessions = max(obj.data.sessionID);
            obj.data_stan.numTestContrasts = 1000;
            obj.data_stan.testContrastLeft = [linspace(1,0,obj.data_stan.numTestContrasts/2)'; zeros(obj.data_stan.numTestContrasts/2,1)];
            obj.data_stan.testContrastRight = [zeros(obj.data_stan.numTestContrasts/2,1); linspace(0,1,obj.data_stan.numTestContrasts/2)'];
            
        end
        
        function simulateAndFit
            %Simulate data and fit model
            error('todo');
        end
        
        function p = getMAP(obj)
            %Get maximum a-posteriori estimate of parameters
            fitObj = stan('fit',obj.stanModelObj,'method','optimize','data',obj.data_stan,'verbose',true);
            fitObj.block;
            
            %Output MAP parameter values
            p = fitObj.extract('permuted',false);
            p = rmfield(p, 'lp__');
            obj.MAP = p;
            
        end
        
        function p = getPosterior(obj)
            
            %Fit model on data
            fitObj = stan('fit',obj.stanModelObj,'method','sample','data',obj.data_stan,'iter',1000,'chains',4,'verbose',true);
            fitObj.block;
            
            %Get all parameter values
            p = fitObj.extract;
            fields = fieldnames(p);
            p = rmfield(p, fields(contains(fields,'__')));
            obj.Posterior = p;
            obj.PosteriorType = 'HMC';
            
            warning('TODO: Assess convergence');
            
%             %Print parameter summaries
%             f = fieldnames(p);
%             params = f(~contains(f,'Test'));
%             for prm = 1:length(params)
%                 dat = p.(params{prm});
%                 fprintf('%s :: mean %.3f :: 2.5%% %.3f :: 97.5%% %.3f\n',params{prm}, mean(dat), quantile(dat,0.025), quantile(dat,0.975))
%             end
%             keyboard;
        end
        
        function p = getPosterior_Variational(obj)
            %Fit model on data
            fitObj = stan('fit',obj.stanModelObj,'method','variational','data',obj.data_stan,'iter',1000,'chains',4,'verbose',true);
            fitObj.block;
            
            %Get all parameter values
            p = fitObj.extract;
            fields = fieldnames(p);
            p = rmfield(p, fields(contains(fields,'__')));
            obj.Posterior = p;
            obj.PosteriorType = 'Variational';
        end
        
        function plotFullPosterior(obj)
            
            %             if ~isempty(obj.Posterior)
            %                 error('posterior not estimated');
            p = obj.Posterior;
            p = struct2array(p);
            pNames = fieldnames(obj.Posterior);
            numVars = structfun(@(f) size(f,2), obj.Posterior);
            %Expand the param labels if any variable in the posterior
            %struct is more than 1D
            pNames_expanded = pNames;
            for param = 1:length(numVars)
                pNames_expanded{param} = repmat(pNames(param),numVars(param),1);
            end
            pNames_expanded = cat(1,pNames_expanded{:});
            
            %Plot posterior distribution
            figure('color','w');
            [S,AX,BigAx,H,HAx] = plotmatrix(p,'k.');
            set(AX,'box','off');
            set(HAx,'box','off');
            set(H,'DisplayStyle','stairs');
            delete(AX(find(triu(ones(length(pNames_expanded)),1))));
            
            for param = 1:length(pNames_expanded)
                xlabel(AX(end,param),pNames_expanded{param},'interpreter','none');
                ylabel(AX(param,1),pNames_expanded{param},'interpreter','none');
            end
        end
        
        function plotPosteriorMarginals(obj)
            
            %             if ~isempty(obj.Posterior)
            %                 error('posterior not estimated');
            p = obj.Posterior;
            pNames = fieldnames(obj.Posterior);
            pNames(end-1:end)=[]; %Remove zTest and pRTest variables
            
            figure('color','w');
            for param = 1:length(pNames)
                ha = subplot(length(pNames),1, param);
                
                temp = p.(pNames{param});
                
                hold(ha,'on');
                for var = 1:size(temp,2)
                    histogram( ha, temp(:,var),...
                        'DisplayStyle','stairs');
                end
                title(ha, pNames{param}, 'interpreter', 'none');
            end
            
        end
        
        function plotPsych(obj)
            %Plot psychometric data/curves. Plotting depends on the model
            %hierarchy type
            
            posterior_credible_intervals = [0.025 0.975];
            
            if isfield(obj.data_stan,'choice')
                error('Plotting not coded for 2AUC data');
            end
            
            p = obj.Posterior;
            m = obj.MAP;
            
            if ~isempty(p)
                fprintf('Posterior estimation: %s\n',obj.PosteriorType);
            end
            
            %Different plots depending on the hierarchy type
            %Non-hierarchical model
            if contains(obj.modelName,'2')
                hier = 2;
            elseif contains(obj.modelName,'1')
                hier = 1;
            else
                hier = 0;
            end
            
            switch(hier)
                case 0 %No hierarchy, just fit all data in a single model
                    ps = []; ms = [];
                    if ~isempty(p)
                        ps = quantile(p.pRTest,posterior_credible_intervals,1);
                        fprintf('Posterior estimation: %s\n',obj.PosteriorType);
                    end
                    
                    if ~isempty(m)
                        ms = m.pRTest;
                    end
                    
                    fig = figure('color','w');
                    ha=axes;
                    obj.util_plotSingle(ha, obj.data_stan, ps, ms);
                    
                case 1 %Per-session deviations
                    %Make a plot for each session, if there aren't too many
                    %sessions, and then make a plot for the grand average
                                        
                    fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495]);
                    numEdgePlots = ceil(sqrt(obj.data_stan.numSessions));
                    for session = 1:obj.data_stan.numSessions
                        ha = subplot(2*numEdgePlots,numEdgePlots,session);

                        %Get subset data
                        sessIdx = obj.data_stan.sessionID==session;
                        data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
                            'testContrastRight',obj.data_stan.testContrastRight,...
                            'contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
                            'contrastRight',obj.data_stan.contrastRight(sessIdx),...
                            'choiceR',obj.data_stan.choiceR(sessIdx));
                        
                        ps = []; ms = [];
                        if ~isempty(p)
                            ps = quantile(squeeze(p.pRTest(:,session,:)),posterior_credible_intervals,1);
                        end
                        if ~isempty(m)
                            ms = m.pRTest(session,:);
                        end
                        obj.util_plotSingle(ha, data_subset, ps, ms);
                        
                        xlabel(ha,''); ylabel(ha,'');
%                         set(ha,'xcolor','none','ycolor','none');
                        set(ha,'xtick','','ytick','');
                        title(ha,sprintf('session %d',session));
                    end
                    
                    %Create grand average plot
                    ha = subplot(2,1,2);
                    ps = []; ms = [];
                    if ~isempty(p)
                        ps = quantile(p.pRTestGrandAverage,posterior_credible_intervals,1);
                    end
                    if ~isempty(m)
                        ms = m.pRTestGrandAverage;
                    end
                    obj.util_plotSingle(ha, obj.data_stan, ps, ms);
                    title(ha,'grand average (pooled data across sessions)');
                    
                case 2 %Per session & per subject deviations
                    
                    fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495]);
                    numEdgePlots = ceil(sqrt(obj.data_stan.numSubjects));
                    
                    for subj = 1:obj.data_stan.numSubjects
                        ha = subplot(2*numEdgePlots,numEdgePlots,subj);

                        %Get subset data
                        subjIdx = obj.data_stan.subjectID==subj;
                        data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
                            'testContrastRight',obj.data_stan.testContrastRight,...
                            'contrastLeft',obj.data_stan.contrastLeft(subjIdx),...
                            'contrastRight',obj.data_stan.contrastRight(subjIdx),...
                            'choiceR',obj.data_stan.choiceR(subjIdx));
                        
                        ps = []; ms = [];
                        if ~isempty(p)
                            keyboard;
                            ps = quantile(p.pRTestSubjectAverage(subj,:),posterior_credible_intervals,1);
                        end
                        
                        if ~isempty(m)
                            ms = m.pRTestSubjectAverage(subj,:);
                        end
                        obj.util_plotSingle(ha, data_subset, ps, ms);
                        xlabel(ha,''); ylabel(ha,'');
                        set(ha,'xcolor','none','ycolor','none');
                        title(ha,sprintf('subject %d',subj));
                    end
                    
                    %Create grand average plot
                    ha = subplot(2,1,2);
                    ps = []; ms = [];
                    if ~isempty(p)
                        ps = quantile(p.pRTestGrandAverage,posterior_credible_intervals,1);
                    end
                    if ~isempty(m)
                        ms = m.pRTestGrandAverage;
                    end
                    obj.util_plotSingle(ha, obj.data_stan, ps, ms);
                    title(ha,'grand average');
                    
            end

            
        end
        
        
        function util_plotSingle(~,axis_handle, dataStruct, PosteriorPrediction, MAPPrediction)
            hold(axis_handle,'on');
            
            cDiff = dataStruct.testContrastRight - dataStruct.testContrastLeft;
            %Plot posterior prediction
            if ~isempty(PosteriorPrediction)
                fx = fill(axis_handle,[cDiff; flipud(cDiff)], [PosteriorPrediction(2,:) fliplr( PosteriorPrediction(1,:) ) ], 'k');
                fx.EdgeAlpha=0;
                fx.FaceColor = [1 1 1]*0.8;
            end
            
            %Plot MAP prediction
            if ~isempty(MAPPrediction)
                mx=plot(axis_handle,cDiff, MAPPrediction, 'k--');
            end
            
            %Plot data
            cDiffData = dataStruct.contrastRight - dataStruct.contrastLeft;
            cVals = unique(cDiffData);
            ph=[];
            for c = 1:length(cVals)
                r = dataStruct.choiceR(cDiffData == cVals(c));
                [ph(c),pci] = binofit(sum(r==1,1),length(r));
                l(1)=line(axis_handle,[1 1]*cVals(c),pci);
                set(l,'Color',[1 1 1]*0,'Linewidth',0.5);
                plot(axis_handle,cVals(c),ph(c),'.','markersize',20,'color',[1 1 1]*0);
            end
            
            set(axis_handle,'xlim',[-1 1],'ylim',[0 1]);
            ylabel(axis_handle,'pR');
            xlabel(axis_handle,'CR - CL');
            %
            %             if ~isempty(PosteriorPrediction) && ~isempty(MAPPrediction)
            %                 lx = legend([fx,mx],{obj.PosteriorType,'MAP'},'Location','SouthEast');
            %                 lx.Box='off';
            %             end
            %
            hold(axis_handle,'off');
            set(axis_handle,'dataaspectratio',[1 1 1]);
        end
        
        
        
    end
    
    
    
end