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
        data_type;
        Z;
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
                load(modelFile, 'sm','Z');
                fprintf('Loaded model %s\n', modelName);
                
            else %If it doesn't exist, compile and save it
                fprintf('Compiled model %s not found, compiling now...\n', modelName);
                
                %Get stan model file
                stanFile = strrep(modelFile, '.mat', '.stan');
                
                %Get Z function
                fid = fopen(stanFile);
                text = textscan(fid, '%s'); text = strjoin(text{1},' ');
                fclose(fid);
                
                %Try loading Z1 and Z2
                Z = {};
                s=strsplit(text,'<Z1>');
                Z{1} = eval(s{2});
                s=strsplit(text,'<Z2>');
                try
                    Z{2} = eval(s{2});
                catch
                end

                %Compile model
                sm = StanModel('file',stanFile,'working_dir',tempdir);
                sm.compile();
               
                %Save copy of compiled model .mat object
                save(modelFile,'sm','Z');
            end
            
            %Store stanModel object in local behavModel object
            obj.stanModelObj = sm;
            obj.Z = Z;
            
            %IDENTIFY FEATURES
            %1) Hierarchical?
            if contains(obj.modelName,'2')
                obj.data_type.hierarchy = 'perSession+perSubject';
            elseif contains(obj.modelName,'1')
                obj.data_type.hierarchy = 'perSession';
            else
                obj.data_type.hierarchy = 'none';
            end
            
            %2) Discrimination or Detection?
            if any(min(obj.data.contrastLeft,obj.data.contrastRight) > 0)
                obj.data_type.stimulusType = 'Discrimination';
            else
                obj.data_type.stimulusType = 'Detection';
            end
            
            %3) Forced or unforced choice?
            if max(obj.data.choice) == 2
                obj.data_type.choiceType = '2AFC';
            elseif max(obj.data.choice) == 3
                obj.data_type.choiceType = '2AUC';
            end
            
            %Convert data structure to stan data format
            obj.data_stan = obj.data;
            obj.data_stan.numTrials = length(obj.data.contrastLeft);
            obj.data_stan.numSubjects = max(obj.data.subjectID);
            obj.data_stan.numSessions = max(obj.data.sessionID);
            
            %Add test contrasts, depending on whether the original data was
            %detection or discrimination
            
            if contains(obj.data_type.stimulusType,'Detection')
                numTestContrasts = 40;
                obj.data_stan.numTestContrasts = numTestContrasts;
                obj.data_stan.testContrastLeft = [linspace(1,0,numTestContrasts/2)'; zeros(numTestContrasts/2,1)];
                obj.data_stan.testContrastRight = [zeros(numTestContrasts/2,1); linspace(0,1,numTestContrasts/2)'];
                
            elseif contains(obj.data_type.stimulusType,'Discrimination')
                
                %Create different pedestals for test contrasts, based on
                %the number of pedestals in the dataset
                cont = [obj.data.contrastLeft obj.data.contrastRight];
                cVals = unique(min(cont,[],2));
                
                if length(cVals)>4
                    %Take top 4 instead
                    tab=tabulate(min(cont,[],2));
                    tab=sortrows(tab,3,'descend');
                    cVals = tab(1:4,1);
                    cVals = sort(cVals);
                end
                
                CL = []; CR = [];
                for ped = 1:length(cVals)
                    c = (cVals(ped):0.01:1)';
                    CL = [CL; c; ones(size(c))*cVals(ped) ];
                    CR = [CR; ones(size(c))*cVals(ped); c];
                end
  
                obj.data_stan.numTestContrasts = numel(CL);
                obj.data_stan.testContrastLeft = CL;
                obj.data_stan.testContrastRight = CR;
            end
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

            p = obj.Posterior;
            pNames = fieldnames(obj.Posterior);
            pNames(find(strcmp(pNames,'z')):end)=[];
            
            figure('color','w');
            ha = tight_subplot( ceil(length(pNames)/3), 3, 0.1, 0.05, 0.05);
            set(ha,'ycolor','none');
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
        
%         function ph = calculatePhat(~,Z,Posterior,contrasts)
%              numTestContrasts = 40;
%             testContrastLeft = [linspace(1,0,numTestContrasts/2)'; zeros(numTestContrasts/2,1)];
%              testContrastRight = [zeros(numTestContrasts/2,1); linspace(0,1,numTestContrasts/2)'];
%              
%              
%              for iter = 1:length(Posterior.bias)
%                  z1 = Z{1}(Posterior.bias(iter), Posterior.sens(iter), testContrastRight,testContrastLeft);
%                  pR = exp(z1)./(1+exp(z1));
%                  
%              end
%             
%             
%             keyboard;
%         end
        
        function plotPsych(obj)
            %Plot psychometric data/curves. Plotting depends on the model
            %hierarchy type
            
            credible_interval = [0.025 0.975];
            
            p = obj.Posterior;
            
            ps = [];
            if ~isempty(p)
                fprintf('Posterior estimation: %s\n',obj.PosteriorType);
%                 ph = obj.calculatePhat(obj.Z,p,'detection');
                
                %Query contrast points
                ps = struct;
                ps.interval = quantile(p.pTest,credible_interval,1);
                ps.mean = mean(p.pTest,1);
            end
            
            %Different plots depending on the hierarchy type
            switch(obj.data_type.hierarchy)
                case 'none' %No hierarchy, just fit all data in a single model
                     
                    fig = figure('color','w','name',obj.modelName);
                    ha=axes;
                    obj.util_plotDetection(ha, obj.data_stan, ps);
%                     obj.util_plotDiscrimination(ha, obj.data_stan, ps);

                case 'perSession' %Per-session deviations
                    %Make a plot for each session, if there aren't too many
                    %sessions, and then make a plot for the grand average
                                        
                    fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495],'name',obj.modelName);
                    numEdgePlots = ceil(sqrt(obj.data_stan.numSessions));
                    ha = tight_subplot( numEdgePlots, numEdgePlots, 0.02, [0.4 0.03], 0.01);
                    for session = 1:obj.data_stan.numSessions

                        %Get subset data
                        sessIdx = obj.data_stan.sessionID==session;
                        data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
                            'testContrastRight',obj.data_stan.testContrastRight,...
                            'contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
                            'contrastRight',obj.data_stan.contrastRight(sessIdx),...
                            'choice',obj.data_stan.choice(sessIdx));
                        
                        ps_thisSess = [];
                        if ~isempty(ps)
                            ps_thisSess = struct('interval', permute(ps.interval(:,:,session,:), [1 2 4 3] ),...
                                                 'mean', permute(ps.mean(:,:,session,:), [1 2 4 3]) );
                        end
                        obj.util_plotDetection(ha(session), data_subset, ps_thisSess );
%                         obj.util_plotDiscrimination(ha(session), data_subset, ps_thisSess);
                        
                        xlabel(ha(session),''); ylabel(ha(session),'');
                        set(ha(session),'xtick','','ytick','');
                        title(ha(session),sprintf('session %d',session));
                    end
                    delete(ha(session+1:end))
                    
                    %Create grand average plot
                    ha = tight_subplot( 1, 1, 0.01, [0.05 0.6], 0.1);

                    ps_grandAv = [];
                    if ~isempty(ps)
                        ps_grandAv = struct('interval', quantile(p.pTestGrandAverage,credible_interval,1),...
                            'mean', mean(p.pTestGrandAverage,1) );
                    end
%                     obj.util_plotDetection(ha, obj.data_stan, ps_grandAv);
                    obj.util_plotDiscrimination(ha, obj.data_stan, ps_grandAv);
                    title(ha,'grand average (pooled data across sessions)');
                    
                otherwise %Per session & per subject deviations
                    error('Not coded. Need prediction from model'); 
            end

            
        end
        
        function util_plotDiscrimination(~,axis_handle, dataStruct, Posterior)
            colours = [ 0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.4940    0.1840    0.5560];
            
            %Get original axis position
            set(axis_handle,'units','normalized');
            bounds = get(axis_handle,'position');
            
            %Delete original axis
%             axis_handle.delete;
            set(axis_handle,'xcolor','none','ycolor','none');
            
            %replace with many subplots within the bounds of the original
            %axis
            ha = tight_subplot(4,3,0, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
            for i = 1:length(ha)
                hold(ha(i),'on');
            end
            
            %Define the four pedestal values in the data
            cont = [dataStruct.testContrastLeft dataStruct.testContrastRight];
            cVals = unique(min(cont,[],2));
            
            %Plot posterior prediction
            if ~isempty(Posterior)
                
                %For each pedestal value
                for ped = 1:length(cVals)
                    %Get subset of test Contrasts which are compatible. If
                    %not exact then find nearest
                    
                    pedIdx = min(dataStruct.testContrastLeft,dataStruct.testContrastRight) == cVals(ped);
                    cDiff_ped = dataStruct.testContrastRight(pedIdx) - dataStruct.testContrastLeft(pedIdx);
                    [cDiff_ped, sortIdx]=sort(cDiff_ped);
                    
                    interval_ped = Posterior.interval(:,pedIdx,:);
                    interval_ped = interval_ped(:,sortIdx,:);
                    mean_ped = Posterior.mean(:,pedIdx,:);
                    mean_ped = mean_ped(:,sortIdx,:);
                    
                    ii=1;
                    for r = [1 3 2]
                        fx = fill(ha( 3*(ped-1) + ii ),[cDiff_ped; flipud(cDiff_ped)], [interval_ped(2,:,r) fliplr( interval_ped(1,:,r) ) ], 'k');
                        fx.FaceAlpha=0.3;
                        fx.EdgeAlpha=0;
                        fx.FaceColor = colours(r,:);
                        plot(ha( 3*(ped-1) + ii ),cDiff_ped, mean_ped(:,:,r), '-', 'color', colours(r,:));
                        ii = ii + 1;
                    end
                    
                end
                
            end
            
            
            for i = 1:length(ha)
                hold(ha(i),'off');
                set(ha(i),'xlim',[-1 1],'ylim',[0 1]);
            end
            
            %             warning('todo');
        end
        
        function util_plotDetection(~,axis_handle, dataStruct, Posterior)
            hold(axis_handle,'on');
            colours = [ 0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.4940    0.1840    0.5560];
            
            %Get subset of data which is compatible with detection
            detectIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight)==0;
            dataStruct.contrastLeft = dataStruct.contrastLeft(detectIdx);
            dataStruct.contrastRight = dataStruct.contrastRight(detectIdx);
            dataStruct.choice = dataStruct.choice(detectIdx);
            if any(~detectIdx)
                warning('DISCRIMINATION contrast conditions in dataset. Only plotting DETECTION subset of the data');
            end
            
            detectIdx = min(dataStruct.testContrastRight,dataStruct.testContrastLeft)==0;
            dataStruct.testContrastLeft = dataStruct.testContrastLeft(detectIdx);
            dataStruct.testContrastRight = dataStruct.testContrastRight(detectIdx);
            dataStruct.numTestContrasts = sum(detectIdx);
            
            %Plot posterior prediction
            if ~isempty(Posterior)
                Posterior.interval = Posterior.interval(:,detectIdx,:);
                Posterior.mean = Posterior.mean(:,detectIdx,:);
                
                cDiff = dataStruct.testContrastRight - dataStruct.testContrastLeft;
                
                %Resort
                [cDiff, sortIdx]=sort(cDiff);
                Posterior.interval = Posterior.interval(:,sortIdx,:);
                Posterior.mean = Posterior.mean(:,sortIdx,:);

                if size(Posterior.interval,3)==1
                    fx = fill(axis_handle,[cDiff; flipud(cDiff)], [Posterior.interval(2,:) fliplr( Posterior.interval(1,:) ) ], 'k');
                    fx.EdgeAlpha=0;
                    fx.FaceColor = [1 1 1]*0.8;
                    
                    plot(axis_handle,cDiff, Posterior.mean, '-');

                elseif size(Posterior.interval,3)==3

                    for r = 1:3
                        fx = fill(axis_handle,[cDiff; flipud(cDiff)], [Posterior.interval(2,:,r) fliplr( Posterior.interval(1,:,r) ) ], 'k');
                        fx.FaceAlpha=0.3;
                        fx.EdgeAlpha=0;
                        fx.FaceColor = colours(r,:);
                        plot(axis_handle,cDiff, Posterior.mean(:,:,r), '-', 'color', colours(r,:));
                    end
                end
                
            end
            
            %Plot data
            cDiffData = dataStruct.contrastRight - dataStruct.contrastLeft;
            cVals = unique(cDiffData);
            for c = 1:length(cVals)
                resp = dataStruct.choice(cDiffData == cVals(c));
                [ph,pci] = binofit(sum([resp==1 resp==2 resp==3],1),length(resp));
                
                if max(dataStruct.choice)==2 %Forced choice
                    L=line(axis_handle,[1 1]*cVals(c),pci(2,:));
                    set(L,'Color',[1 1 1]*0,'Linewidth',0.5);
                    plot(axis_handle,cVals(c),ph(2),'.','markersize',20,'color',[1 1 1]*0);
                elseif max(dataStruct.choice)==3 %Unforced choice
                    
                    for r = 1:3
                        L=line(axis_handle,[1 1]*cVals(c),pci(r,:));
                        set(L,'Color',colours(r,:),'Linewidth',0.5);
                        plot(axis_handle,cVals(c),ph(r),'.','markersize',20,'color',colours(r,:));
                    end
                    
                end
                
                
            end
            
            set(axis_handle,'xlim',[-1 1],'ylim',[0 1]);
            ylabel(axis_handle,'pR');
            xlabel(axis_handle,'CR - CL');
            hold(axis_handle,'off');
%             set(axis_handle,'dataaspectratio',[1 1 1]);
        end
        
    end
    
    
    
end