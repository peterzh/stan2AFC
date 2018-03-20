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
            
            fig = figure('color','w');
            ha = tight_subplot(4,2,0.05,[0.05 0.05],[0.05 0.01]);
            set(ha([1 3 5 7]),'dataaspectratio',[1 1 1],'xlim',[0 1],'ylim',[0 1]);
            set(ha,'xcolor','none','ycolor','none');
            
            cont = [obj.data.contrast_left obj.data.contrast_right];
            resp = obj.data.choiceR;
            %                 cVals = unique(cont(:));
            cVals = unique(min(cont,[],2));
            
            if length(cVals)>4
                %Take top 4 instead
                tab=tabulate(min(cont,[],2));
                tab=sortrows(tab,3,'descend');
                cVals = tab(1:4,1);
                cVals = sort(cVals);
                %                     keyboard;
            end
            
            numPedestals = length(cVals);
            cols = [0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.4940    0.1840    0.5560];
            
            set(ha,'xlim',[-1 1]*(max(cont(:))+0.1),'ylim',[0 1]);
            
            for ped = 1:numPedestals
                
                %Add cartoon on the left of the contrast conditions
                %plotted
                ha_cartoon = ha( 2*(ped-1) + 1);
                ped_idx = min(cont,[],2)==cVals(ped);
                ped_c = unique(cont(ped_idx,:),'rows');
                
                plot(ha_cartoon, ped_c(:,2), ped_c(:,1), 'k.','markersize',15);
                set(ha_cartoon,'dataaspectratio',[1 1 1],'xlim',[0 max(cont(:))+0.1],'ylim',[0 max(cont(:))+0.1],'box','off');
                set(ha_cartoon,'ytick',unique(ped_c(:)),'xtick',unique(ped_c(:)));
                
                ha_ped = [ ha(2*(ped-1) + 2) ];
                set(ha_ped(1),'yticklabelmode','auto','ycolor','k');
                if ped == numPedestals
                    set(ha_ped,'xticklabelmode','auto','xcolor','k');
                end
                
                hold(ha_ped(1),'on');
                set(ha_ped,'colororder',cols);
                
                %Plot actual datapoints
                ped_idx = min(cont,[],2)==cVals(ped);
                ped_c_diff = diff(cont(ped_idx,:),[],2);
                ped_r = resp(ped_idx);
                uC = unique(ped_c_diff);
                ph=[];
                for c = 1:length(uC)
                    r = ped_r(ped_c_diff==uC(c));
                    [ph(c,:),pci] = binofit(sum(r==1,1),length(r));
                    l(1)=line(ha_ped(1),[1 1]*uC(c),pci);
                    set(l,'Color',cols(2,:),'Linewidth',0.5);

                end
%                 set(ha(ped),'ColorOrderIndex',1);
                    plot(ha_ped(1),uC,ph,'.','markersize',15,'color',cols(2,:));

                
                %Plot predictions
                cEval = [linspace(max(abs(uC))+0.1,0,1000)' zeros(1000,1); zeros(1000,1) linspace(0,max(abs(uC))+0.1,1000)'] + cVals(ped);
                numIter = length(p.bias);
                pR = nan(numIter, 2000);
                for iter = 1:numIter
                    pSet = structfun(@(f) f(iter), p, 'uni', 0);
                    zSet = obj.z{1}(pSet, cEval(:,1), cEval(:,2));
                    pR(iter,:) = 1./(1+exp(-zSet));
                end
                bounds = quantile(pR,[0.025 0.975]);
                cDiff = cEval(:,2) - cEval(:,1);
                %
                fx = fill(ha_ped(1), [cDiff; flipud(cDiff)], [bounds(2,:) fliplr( bounds(1,:) ) ], 'k');
                fx.FaceAlpha = 0.2;
                
                zSet = obj.z{1}(obj.MAP, cEval(:,1), cEval(:,2));
                pR = 1./(1+exp(-zSet));
                plot(ha_ped(1),cDiff, pR, 'k--');
                
                
                
%                 p_hat = obj.calculatePhat(obj.parameterFits,testCont);
%                 set(gca,'ColorOrderIndex',1);
%                 for ch=1:3
%                     plot(ha_ped(ch),diff(testCont,[],2),p_hat(:,ch),'linewidth',1,'color',cols(ch,:));
%                 end
%                 

                
                if ped==1
                    title(ha_ped(1),'pRight');
                end
                
                
                
            end
            
%             set(get(fig,'children'),'fontsize',6);
            
            
            
        end
        
        
        
        
    end
    
    
    
end