classdef behavModel < handle
    properties(Access=public)
        modelName;
        data;        
        MAP;
        Posterior;
        PosteriorType;
    end
    
    properties(Access=public)
        stanModelObj;
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
                load(modelFile,'sm');
                fprintf('Loaded model %s\n', modelName);
                
            else %If it doesn't exist, compile and save it
                fprintf('Compiled model %s not found, compiling now...\n', modelName);
                
                %Get stan model file
                stanFile = strrep(modelFile, '.mat', '.stan');
                
                %Compile model
                sm = StanModel('file',stanFile,'working_dir',tempdir);
                %                 sm = StanModel('file',stanFile,'working_dir',fileparts(stanFile));
                %                 sm.compile('model','STANCFLAGS =include_paths C:\Users\Peter\Documents\MATLAB\stan2AFC\models');
                sm.compile();
                
                
                %Save copy of compiled model .mat object
                save(modelFile,'sm');
            end
            
            %Store stanModel object in local behavModel object
            obj.stanModelObj = sm;
            
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
        
        function [p1,waic] = getPosterior(obj)
            
            %Fit model on data
            fitObj = stan('fit',obj.stanModelObj,'method','sample','data',obj.data_stan,'iter',1000,'chains',4,'verbose',true);
            fitObj.block;
                        
            %Get all parameter values
            p = fitObj.extract('permuted',false);
            
            %Convert to scalar struct
            fields = fieldnames(p);
            p1 = p(1);
            for i = 1:length(fields)
                p1.(fields{i}) = cat(1,p.(fields{i}));
            end
            
            p1 = rmfield(p1, fields(contains(fields,'__')));
            obj.Posterior = p1;
            obj.PosteriorType = 'HMC';
            
            warning('TODO: Assess convergence');
%             [str,tab] = fitObj.print();

            waic=mstan.waic(p1.log_lik);

        end
        
        function p = getPosterior_Variational(obj)
            %Fit model on data
            fitObj = stan('fit',obj.stanModelObj,'method','variational','data',obj.data_stan,'iter',1000,'chains',4,'verbose',true);
            fitObj.block;
            
            %Get all parameter values
            p = fitObj.extract('permuted',false);
            fields = fieldnames(p);
            p = rmfield(p, fields(contains(fields,'__')));
            obj.Posterior = p;
            obj.PosteriorType = 'Variational';
        end
        
        function plotPosteriorFull(obj)
            
            if isempty(obj.Posterior)
                error('posterior not estimated');
            end
                        
            p = obj.Posterior;
            pNames = fieldnames(p);
            p = rmfield(p, pNames(find(strcmp(pNames,'z')):end));
            pNames = fieldnames(p);
            
            
            %Expand the param labels if any variable in the posterior
            %struct is more than 1D
            pNames_expanded = pNames;
            for param = 1:length(pNames)
                numVars = size(p.(pNames{param}),2);
                pNames_expanded{param} = repmat(pNames(param),numVars,1);
            end
            pNames_expanded = cat(1,pNames_expanded{:});
            
            %Plot posterior distribution
            figure('color','w');
            [S,AX,BigAx,H,HAx] = plotmatrix(struct2array(p),'k.');
            set(AX,'box','off');
            set(HAx,'box','off');
            set(H,'DisplayStyle','bar', 'edgealpha', 0);
            delete(AX(find(triu(ones(length(pNames_expanded)),1))));
            
            for param = 1:length(pNames_expanded)
                xlabel(AX(end,param),pNames_expanded{param},'interpreter','none', 'Rotation', 45);
                ylabel(AX(param,1),pNames_expanded{param},'interpreter','none', 'Rotation', 45);
            end
        end
        
        function plotPosteriorMarginals(obj)
            
            %             if ~isempty(obj.Posterior)
            %                 error('posterior not estimated');
            p = obj.Posterior;
            pNames = fieldnames(obj.Posterior);
            pNames(find(strcmp(pNames,'z')):end)=[];
            
            figure('color','w');
            ha = tight_subplot( ceil(length(pNames)/3), 3, 0.05, 0.05, 0.05);
            for param = 1:length(pNames)
%                 ha = subplot(length(pNames),1, param);
                
                temp = p.(pNames{param});
                
                hold(ha(param),'on');
                for var = 1:size(temp,2)
                    histogram( ha(param), temp(:,var),...
                        'DisplayStyle','bar', 'edgealpha', 0);
                end
                title(ha(param), pNames{param}, 'interpreter', 'none');
            end
            set(ha,'XTickLabelMode','auto');
            delete(ha(param+1:end));
            
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
                    
                    fig = figure('color','w','name',obj.modelName);
                    ha=axes;
                    obj.util_plotSingle(ha, obj.data_stan, ps, ms);
                    
                case 1 %Per-session deviations
                    %Make a plot for each session, if there aren't too many
                    %sessions, and then make a plot for the grand average
                                        
                    fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495],'name',obj.modelName);
                    numEdgePlots = ceil(sqrt(obj.data_stan.numSessions));
                    ha = tight_subplot( numEdgePlots, numEdgePlots, 0.02, [0.4 0.03], 0.01);
                    for session = 1:obj.data_stan.numSessions
%                         ha = subplot(2*numEdgePlots,numEdgePlots,session);

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
                        obj.util_plotSingle(ha(session), data_subset, ps, ms);
                        
                        xlabel(ha(session),''); ylabel(ha(session),'');
%                         set(ha,'xcolor','none','ycolor','none');
                        set(ha(session),'xtick','','ytick','');
                        title(ha(session),sprintf('session %d',session));
                    end
                    delete(ha(session+1:end))
                    
                    %Create grand average plot
                    ha = tight_subplot( 1, 1, 0.01, [0.05 0.6], 0.1);
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
                    
                    fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495],'name',obj.modelName);
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
                            ps = quantile(squeeze(p.pRTestSubjectAverage(:,subj,:)),posterior_credible_intervals,1);
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
%             set(axis_handle,'dataaspectratio',[1 1 1]);
        end
        
        
        
    end
    
    
    
end