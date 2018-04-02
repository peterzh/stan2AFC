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
            
%             obj.getMAP;
%             obj.getPosterior;
        end
        
        function simulateAndFit
            %Simulate data and fit model
            error('todo');
        end
        
        function p = getMAP(obj)
            %Get maximum a-posteriori estimate of parameters
            fitObj = stan('fit',obj.stanModelObj,'method','optimize','data',obj.data,'verbose',true);
            fitObj.block;
            
            %Output MAP parameter values
            p = fitObj.extract('permuted',false);
            p = rmfield(p, 'lp__');
            obj.MAP = p;
            
        end
        
        function p = getPosterior(obj)
            
            %Fit model on data
            fitObj = stan('fit',obj.stanModelObj,'method','sample','data',obj.data,'iter',1000,'chains',4,'verbose',true);
            fitObj.block;
            
            %Get all parameter values
            p = fitObj.extract;
            fields = fieldnames(p);
            p = rmfield(p, fields(contains(fields,'__')));
            obj.Posterior = p;
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
            %Plot per-session curves (+data) and then the grand average
            %curves (no data)
            p = obj.Posterior;
            
            fig = figure('color','w');
            
            numPlots = obj.data.numSessions+1;
            
            for sess = 1:obj.data.numSessions
                subplot( ceil(sqrt(numPlots)) , ceil(sqrt(numPlots)), sess);
                hold on; set(gca,'box','off','xlim',[-1 1],'ylim',[0 1]);
                
                sessIdx = obj.data.sessionID == sess;
                %Plot data
                cont = [obj.data.contrast_left(sessIdx) obj.data.contrast_right(sessIdx)];
                contDiff = cont(:,2) - cont(:,1);
                resp = obj.data.choice(sessIdx);
                cVals = unique(contDiff);
                ph=[];
                for c = 1:length(cVals)
                    r = resp(contDiff == cVals(c));
                    [ph(c,:),pci] = binofit(sum(r==2,1),length(r));
                    l(1)=line(gca,[1 1]*cVals(c),pci);
                    set(l,'Color',[1 1 1]*0.5,'Linewidth',0.5);
                end
                plot(cVals,ph,'.','markersize',15,'color',[1 1 1]*0.5);
                
                
                %Plot model prediction
                cEval = [linspace(1,0,1000)' zeros(1000,1); zeros(1000,1) linspace(0,1,1000)'];
                numIter = length(p.bias);
                pR = nan(numIter, 2000);
                for iter = 1:numIter
                    pSet = structfun(@(f) f(iter,:), p, 'uni', 0);
                    
                    pSetSess = struct;
                    pSetSess.bias = pSet.bias + pSet.bias_delta_perSession(sess);
                    pSetSess.sens = pSet.sens + pSet.sens_delta_perSession(sess);
                    zSet = obj.z{1}(pSetSess, cEval(:,1), cEval(:,2));
                    pR(iter,:) = 1./(1+exp(-zSet));
                end
                bounds = quantile(pR,[0.025 0.975]);
                cDiff = cEval(:,2) - cEval(:,1);
                %
                fx = fill([cDiff; flipud(cDiff)], [bounds(2,:) fliplr( bounds(1,:) ) ], 'k');
                fx.FaceAlpha = 0.2;
                
%                 zSet = obj.z{1}(obj.MAP, cEval(:,1), cEval(:,2));
%                 pR = 1./(1+exp(-zSet));
%                 plot(ha_ped(1),cDiff, pR, 'k--');
            end
            
            
            %Now plot grand avg
            subplot( ceil(sqrt(numPlots)) , ceil(sqrt(numPlots)), numPlots);
            
            pR = nan(numIter, 2000);
            for iter = 1:numIter
                pSet = structfun(@(f) f(iter,:), p, 'uni', 0);
                
                pSetSess = struct;
                pSetSess.bias = pSet.bias ;
                pSetSess.sens = pSet.sens ;
                zSet = obj.z{1}(pSetSess, cEval(:,1), cEval(:,2));
                pR(iter,:) = 1./(1+exp(-zSet));
            end
            bounds = quantile(pR,[0.025 0.975]);
            cDiff = cEval(:,2) - cEval(:,1);
            %
            fx = fill([cDiff; flipud(cDiff)], [bounds(2,:) fliplr( bounds(1,:) ) ], 'k');
            fx.FaceAlpha = 0.2;
            
            
        end
        
        
        
        
    end
    
    
    
end