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
            
            obj.getMAP;
            obj.getPosterior;
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
            p = fitObj.extract;
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
        
        function plotPosterior(obj)
            
            %             if ~isempty(obj.Posterior)
            %                 error('posterior not estimated');
            p = obj.Posterior;
            pNames = fieldnames(obj.Posterior);
            
            %Plot posterior distribution
            figure('color','w');
            [S,AX,BigAx,H,HAx] = plotmatrix(struct2array(p),'k.');
            set(AX,'box','off');
            set(HAx,'box','off');
            set(H,'DisplayStyle','stairs');
            delete(AX(find(triu(ones(length(pNames)),1))));
            
            for param = 1:length(pNames)
                xlabel(AX(end,param),pNames{param},'interpreter','none');
                ylabel(AX(param,1),pNames{param},'interpreter','none');
            end
        end
        
        function plotPsych(obj)
            p = obj.Posterior;
            cDiff = obj.data.contrast_right - obj.data.contrast_left;
            %Plot data
            figure('color','w');
            axes; hold on;
            cVal = unique(cDiff);
            for c = 1:length(cVal)
                r = obj.data.choiceR( cDiff == cVal(c));
                plot(cVal(c),mean(r),'k.','markersize',20);
            end
            xlabel('CR - CL');
            ylabel('pR');
            title(sprintf('N=%d',obj.data.N));
            
            
            %Plot psychometric fit
            cEval = [linspace(1,0,1000)' zeros(1000,1); zeros(1000,1)  linspace(0,1,1000)'];
            numIter = length(p.bias);
            pR = nan(numIter, 2000);
            for iter = 1:numIter
                pSet = structfun(@(f) f(iter), p, 'uni', 0);
                zSet = obj.z{1}(pSet, cEval(:,1), cEval(:,2));
                pR(iter,:) = 1./(1+exp(-zSet));
            end
            bounds = quantile(pR,[0.025 0.975]);
            cDiff = cEval(:,2) - cEval(:,1);
            
            fx = fill( [cDiff; flipud(cDiff)], [bounds(2,:) fliplr( bounds(1,:) ) ], 'k');
            fx.FaceAlpha = 0.2;
            
            
            %Now draw MAP estimate
            zSet = obj.z{1}(obj.MAP, cEval(:,1), cEval(:,2));
            pR = 1./(1+exp(-zSet));
            plot(cDiff, pR, 'k--');
        end
        
        
        
        
    end
    
    
    
end