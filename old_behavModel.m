      

%         function plotPsych(obj)
%             %Plot psychometric data/curves.
%             
%             credible_interval = [0.025 0.975];
%             
%             p = obj.Posterior;
%             
%             ps = [];
%             if ~isempty(p)
%                 fprintf('Posterior estimation: %s\n',obj.PosteriorType);
%                 %                 ph = obj.calculatePhat(obj.Z,p,'detection');
%                 
%                 %Query contrast points
%                 ps = struct;
%                 ps.interval = quantile(p.pTest,credible_interval,1);
%                 ps.mean = mean(p.pTest,1);
%             end
%             
%             %One figure per subject
%             prop=[]; ii=1;
%             for subj = 1:obj.data_stan.numSubjects
%                 fig = figure('color','w','units','normalized','position',[0.2131 0.0524 0.4530 0.8495],'name',[obj.modelName ' subject ' num2str(subj)]);
%                 
%                 sessionIDs = find(obj.data_stan.subjID_session==subj);
%                 numEdgePlots = ceil(sqrt(length(sessionIDs)));
%                 ha = tight_subplot( numEdgePlots, numEdgePlots, 0.02, 0.01, 0.01);
%                 for session = 1:length(sessionIDs)
%                     
%                     %Get subset data
%                     sessIdx = obj.data_stan.sessionID==sessionIDs(session);
%                     
%                     %If perturbation data exists, then only include p=0
%                     %data
%                     if isfield(obj.data_stan,'perturbation')
%                         sessIdx = sessIdx & obj.data_stan.perturbation==0;
%                     end
%                     
%                     data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
%                         'testContrastRight',obj.data_stan.testContrastRight,...
%                         'contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
%                         'contrastRight',obj.data_stan.contrastRight(sessIdx),...
%                         'choice',obj.data_stan.choice(sessIdx));
%                     
%                     ps_thisSess = [];
%                     if ~isempty(ps)
%                         if ndims(ps.interval) > 3
%                         ps_thisSess = struct('interval', permute(ps.interval(:,:,sessionIDs(session),:), [1 2 4 3] ),...
%                             'mean', permute(ps.mean(:,:,sessionIDs(session),:), [1 2 4 3]) );
%                         else
%                             ps_thisSess = ps;
%                         end
%                     end
%                     %                         obj.util_plotDetection(ha(session), data_subset, ps_thisSess );
%                     obj.util_plotDiscrimination(ha(session), data_subset, ps_thisSess);
%                     title(ha(session),sprintf('Mouse %d (sess %d)',subj,sessionIDs(session)));
%                 end
%                 delete(ha(session+1:end))
%                 set(findobj(fig,'Type','Line'),'markersize',5);
%                 
%             end
%             
%             %If subject average exists, plot that
%             if isfield(p,'pTest_subjectAverage')
%                 colours = [ 0    0.4470    0.7410;
%                     0.8500    0.3250    0.0980;
%                     0.4940    0.1840    0.5560];
%                 
%                     ps = struct;
%                     ps.interval = quantile(p.pTest_subjectAverage,credible_interval,1);
%                     ps.mean = mean(p.pTest_subjectAverage,1);
%                 
%                 fig = figure('color','w','name',[obj.modelName ' subject average']);
%                 haSubj = tight_subplot(2,3,0.05,0.01,0.01);
%                 for subj = 1:obj.data_stan.numSubjects
%                     
%                     data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
%                         'testContrastRight',obj.data_stan.testContrastRight);
%                     
% %                     fig = figure('color','w','name',[obj.modelName ' subject average ' num2str(subj)]);
%                     
%                     psSubj=struct('interval', permute(ps.interval(:,:,subj,:),[1 2 4 3]),...
%                                    'mean',permute(ps.mean(:,:,subj,:),[1 2 4 3]) );
%                     haGrand=obj.util_plotDiscrimination(haSubj(subj), data_subset, psSubj);
%                     
%                     %Compute the grand average probabilities
%                     cont = [obj.data_stan.testContrastLeft obj.data_stan.testContrastRight];
%                     cVals = unique(min(cont,[],2));
%                     
%                     %Plot the averages across sessions for this subject
%                     sessThisSubject = unique(obj.data_stan.sessionID(obj.data_stan.subjectID==subj));
%                     for sess = 1:length(sessThisSubject)
%                         
%                         %Get subset data
%                         sessIdx = obj.data_stan.sessionID==sessThisSubject(sess);
%                         
%                         %If perturbation data exists, then only include p=0
%                         %data
%                         if isfield(obj.data_stan,'perturbation')
%                             sessIdx = sessIdx & obj.data_stan.perturbation==0;
%                         end
%                         dataStruct = struct('contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
%                             'contrastRight',obj.data_stan.contrastRight(sessIdx),...
%                             'choice',obj.data_stan.choice(sessIdx));
%                         
%                         for ped = 1:length(cVals)
%                             
%                             pedIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight) == cVals(ped);
%                             cDiff_ped = dataStruct.contrastRight(pedIdx) - dataStruct.contrastLeft(pedIdx);
%                             resp = dataStruct.choice(pedIdx);
%                             
%                             cDiff_ped_unique = unique(cDiff_ped);
%                             for c = 1:length(cDiff_ped_unique)
%                                 R = resp(cDiff_ped == cDiff_ped_unique(c));
%                                 [ph,pci] = binofit(sum([R==1 R==2 R==3],1),length(R));
%                                 
%                                 ii=1;
%                                 for r = [1 3 2]
%                                     hold(haGrand( 3*(ped-1) + ii ),'on');
%                                     px=plot(haGrand( 3*(ped-1) + ii ), cDiff_ped_unique(c), ph(r), '.','markersize',10,'Color',colours(r,:));
%                                     ii = ii + 1;
%                                 end
%                                 
%                             end
%                             
%                         end
%                     end
%                     
%                     %Now add the mean of means
%                     for i = 1:length(haGrand)
%                         dots = findobj(haGrand(i),'Marker','.');
%                         
%                         xx = cell2mat(get(dots,'XData'));
%                         yy = cell2mat(get(dots,'YData'));
%                         
%                         cD = unique(xx);
%                         for c = 1:length(cD)
%                             plot(haGrand(i), cD(c), mean(yy(xx==cD(c))), '.','markersize',20, 'Color', [1 1 1]*0 );
%                             lx=line(haGrand(i), [1 1]*cD(c), quantile(yy(xx==cD(c)),[0.25 0.75]), 'Color', [1 1 1]*0 );
%                         end
%                     end
%                     
%                     %Now add axis labels
%                     set(haGrand(10:12),'XTickLabelMode','auto');
%                     set(haGrand([1 4 7 10]),'YTickLabelMode','auto');
%                 end
%             end
%             
%             %If grand average exists, plot that
%             if isfield(p,'pTest_grandAverage')
%                 colours = [ 0    0.4470    0.7410;
%                     0.8500    0.3250    0.0980;
%                     0.4940    0.1840    0.5560];
%                 
%                 ps = struct;
%                 ps.interval = quantile(p.pTest_grandAverage,credible_interval,1);
%                 ps.mean = mean(p.pTest_grandAverage,1);
%                 
%                 data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
%                     'testContrastRight',obj.data_stan.testContrastRight);
%                 
%                 fig = figure('color','w','name',[obj.modelName ' grand average ']);
%                 ax=axes;
%                 haGrand=obj.util_plotDiscrimination(ax, data_subset, ps);
%                 
%                 %Compute the grand average probabilities
%                 cont = [obj.data_stan.testContrastLeft obj.data_stan.testContrastRight];
%                 cVals = unique(min(cont,[],2));
%                 
%                 %Plot the averages across sessions and subjects
%                 for sess = 1:obj.data_stan.numSessions
%                     
%                     %Get subset data
%                     sessIdx = obj.data_stan.sessionID==sess;
%                     
%                     %If perturbation data exists, then only include p=0
%                     %data
%                     if isfield(obj.data_stan,'perturbation')
%                         sessIdx = sessIdx & obj.data_stan.perturbation==0;
%                     end
%                     dataStruct = struct('contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
%                         'contrastRight',obj.data_stan.contrastRight(sessIdx),...
%                         'choice',obj.data_stan.choice(sessIdx));
%                     
%                     for ped = 1:length(cVals)
%                         
%                         pedIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight) == cVals(ped);
%                         cDiff_ped = dataStruct.contrastRight(pedIdx) - dataStruct.contrastLeft(pedIdx);
%                         resp = dataStruct.choice(pedIdx);
%                         
%                         cDiff_ped_unique = unique(cDiff_ped);
%                         for c = 1:length(cDiff_ped_unique)
%                             R = resp(cDiff_ped == cDiff_ped_unique(c));
%                             [ph,pci] = binofit(sum([R==1 R==2 R==3],1),length(R));
%                             
%                             ii=1;
%                             for r = [1 3 2]
%                                 hold(haGrand( 3*(ped-1) + ii ),'on');
%                                 px=plot(haGrand( 3*(ped-1) + ii ), cDiff_ped_unique(c), ph(r), '.','markersize',10,'Color',colours(r,:));
%                                 ii = ii + 1;
%                             end
%                             
%                         end
%                         
%                     end
%                 end
%                 
%                 %Now add the mean of means
%                 for i = 1:length(haGrand)
%                     dots = findobj(haGrand(i),'Marker','.');
%                     
%                     xx = cell2mat(get(dots,'XData'));
%                     yy = cell2mat(get(dots,'YData'));
%                     
%                     cD = unique(xx);
%                     for c = 1:length(cD)
%                         plot(haGrand(i), cD(c), mean(yy(xx==cD(c))), '.','markersize',20, 'Color', [1 1 1]*0 );
%                         lx=line(haGrand(i), [1 1]*cD(c), quantile(yy(xx==cD(c)),[0.25 0.75]), 'Color', [1 1 1]*0 );
%                     end
%                 end
%                 
%                 %Now add axis labels
%                 set(haGrand(10:12),'XTickLabelMode','auto');
%                 set(haGrand([1 4 7 10]),'YTickLabelMode','auto');
%             end
%             
%             if isfield(p,'pTest_grandAverage_perturbation')
%                 numRegions = size(p.pTest_grandAverage_perturbation,3);
%                 ps = struct;
%                 ps.interval = quantile(p.pTest_grandAverage,credible_interval,1);
%                 ps.mean = mean(p.pTest_grandAverage,1);
%                 data_subset = struct('testContrastLeft',obj.data_stan.testContrastLeft,...
%                     'testContrastRight',obj.data_stan.testContrastRight);
%                 
%                 fig = figure('color','w','name',[obj.modelName ' grand average perturbation']);
%                 ha=tight_subplot(2,ceil(numRegions/2),0.05,0.01,0.01);
%                
%                 for i = 1:numRegions
%                     haGrand=obj.util_plotDiscrimination(ha(i), data_subset, ps);
%                     set(findobj(haGrand,'type','line'),'color',[0 0 0]);
%                     set(findobj(haGrand,'type','patch'),'FaceColor',[0 0 0]);
%                     
%                     %add perturbation fit
%                     pt = struct;
%                     pt.interval = quantile( permute(p.pTest_grandAverage_perturbation(:,:,i,:),[1 2 4 3]) ,credible_interval,1);
%                     pt.mean = mean(permute(p.pTest_grandAverage_perturbation(:,:,i,:),[1 2 4 3]),1);
%                 
%                     fTemp=figure;
%                     ax=axes;
%                     haPerturb=obj.util_plotDiscrimination(ax, data_subset, pt);
%                     set(findobj(haPerturb,'type','line'),'color',[1 0 0]);
%                     set(findobj(haPerturb,'type','patch'),'FaceColor',[1 0 0]);
%                     
%                     %Now add data
%                     %Overlay data, average over all sessions and subjects:
%                     cont = [obj.data_stan.testContrastLeft obj.data_stan.testContrastRight];
%                     cVals = unique(min(cont,[],2));
%                     for session = 1:obj.data_stan.numSessions
%                         
%                         %Get subset data
%                         sessIdx = obj.data_stan.sessionID==session & obj.data_stan.perturbation==i;
%                         
%                         dataStruct = struct('contrastLeft',obj.data_stan.contrastLeft(sessIdx),...
%                             'contrastRight',obj.data_stan.contrastRight(sessIdx),...
%                             'choice',obj.data_stan.choice(sessIdx));
%                         
%                         for ped = 1:length(cVals)
%                             
%                             pedIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight) == cVals(ped);
%                             cDiff_ped = dataStruct.contrastRight(pedIdx) - dataStruct.contrastLeft(pedIdx);
%                             resp = dataStruct.choice(pedIdx);
%                             
%                             cDiff_ped_unique = unique(cDiff_ped);
%                             for c = 1:length(cDiff_ped_unique)
%                                 R = resp(cDiff_ped == cDiff_ped_unique(c));
%                                 
%                                 if ~isempty(R)
%                                     [ph,pci] = binofit(sum([R==1 R==2 R==3],1),length(R));
%                                     
%                                     ii=1;
%                                     for r = [1 3 2]
%                                         hold(haPerturb( 3*(ped-1) + ii ),'on');
%                                         px=plot(haPerturb( 3*(ped-1) + ii ), cDiff_ped_unique(c), ph(r), '.','markersize',10,'Color',[1 0 0]);
%                                         ii = ii + 1;
%                                     end
%                                 end
%                             end
%                             
%                         end
%                     end
%                     
%                     for iii = 1:length(haPerturb)
%                         dots = findobj(haPerturb(iii),'Marker','.');
%                         
%                         xx = cell2mat(get(dots,'XData'));
%                         yy = cell2mat(get(dots,'YData'));
%                         
%                         delete(dots);
%                         
%                         cD = unique(xx);
%                         for c = 1:length(cD)
%                             plot(haPerturb(iii), cD(c), mean(yy(xx==cD(c))), '.','markersize',10, 'Color', [1 0 0] );
% %                             lx=line(haPerturb(iii), [1 1]*cD(c), quantile(yy(xx==cD(c)),[0.025 0.975]), 'Color', [1 1 1]*0 );
%                         end
%                     end
%                     
%                     
%                     
%                     for ii = 1:length(haPerturb)
%                         Lstuff = get(haPerturb(ii),'children');                        
%                         copyobj(Lstuff, haGrand(ii))
%                     end
%                     
%                     close(fTemp)
%                     
%                  
%                 end
%             end
%             
%             
%             
%         end
%         
%         function ha=util_plotDiscrimination(~,axis_handle, dataStruct, Posterior)
%             colours = [ 0    0.4470    0.7410;
%                 0.8500    0.3250    0.0980;
%                 0.4940    0.1840    0.5560];
%             
%             %Get original axis position
%             set(axis_handle,'units','normalized');
%             bounds = get(axis_handle,'position');
%             
%             set(axis_handle,'xcolor','none','ycolor','none');
%             
%             %replace with many subplots within the bounds of the original
%             %axis
%             ha = tight_subplot(4,3,0.01, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
%             
%             set(ha(1:9),'xcolor','none');
%             set(ha([2 3 5 6 8 9 11 12]),'ycolor','none');
%             
%             for i = 1:length(ha)
%                 hold(ha(i),'on');
%             end
%             
%             %Define the four pedestal values in the data
%             cont = [dataStruct.testContrastLeft dataStruct.testContrastRight];
%             cVals = unique(min(cont,[],2));
%             
%             %Plot posterior prediction
%             if ~isempty(Posterior)
%                 
%                 %For each pedestal value
%                 for ped = 1:length(cVals)
%                     %Get subset of test Contrasts which are compatible. If
%                     %not exact then find nearest
%                     
%                     pedIdx = min(dataStruct.testContrastLeft,dataStruct.testContrastRight) == cVals(ped);
%                     cDiff_ped = dataStruct.testContrastRight(pedIdx) - dataStruct.testContrastLeft(pedIdx);
%                     [cDiff_ped, sortIdx]=sort(cDiff_ped);
%                     
%                     interval_ped = Posterior.interval(:,pedIdx,:);
%                     interval_ped = interval_ped(:,sortIdx,:);
%                     mean_ped = Posterior.mean(:,pedIdx,:);
%                     mean_ped = mean_ped(:,sortIdx,:);
%                     
%                     ii=1;
%                     for r = [1 3 2]
%                         fx = fill(ha( 3*(ped-1) + ii ),[cDiff_ped; flipud(cDiff_ped)], [interval_ped(2,:,r) fliplr( interval_ped(1,:,r) ) ], 'k');
%                         fx.FaceAlpha=0.3;
%                         fx.EdgeAlpha=0;
%                         fx.FaceColor = colours(r,:);
%                         plot(ha( 3*(ped-1) + ii ),cDiff_ped, mean_ped(:,:,r), '-', 'color', colours(r,:));
%                         ii = ii + 1;
%                     end
%                     
%                 end
%                 
%             end
%             
%             if isfield(dataStruct,'contrastLeft')
%                 %Plot original data
%                 for ped = 1:length(cVals)
%                     
%                     pedIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight) == cVals(ped);
%                     cDiff_ped = dataStruct.contrastRight(pedIdx) - dataStruct.contrastLeft(pedIdx);
%                     resp = dataStruct.choice(pedIdx);
%                     
%                     cDiff_ped_unique = unique(cDiff_ped);
%                     for c = 1:length(cDiff_ped_unique)
%                         R = resp(cDiff_ped == cDiff_ped_unique(c));
%                         [ph,pci] = binofit(sum([R==1 R==2 R==3],1),length(R));
%                         
%                         ii=1;
%                         for r = [1 3 2]
%                             
%                             px=plot(ha( 3*(ped-1) + ii ), cDiff_ped_unique(c), ph(r), 'k.', 'markersize',20);
%                             
%                             ii = ii + 1;
%                         end
%                         
%                     end
%                     
%                 end
%             end
%             
%             for i = 1:length(ha)
%                 hold(ha(i),'off');
%                 set(ha(i),'xlim',[-1 1],'ylim',[0 1]);
%             end
%             
%         end
%         
%         function util_plotDetection(~,axis_handle, dataStruct, Posterior)
%             hold(axis_handle,'on');
%             colours = [ 0    0.4470    0.7410;
%                 0.8500    0.3250    0.0980;
%                 0.4940    0.1840    0.5560];
%             
%             %remove any perturbation data if it exists
%             if isfield(dataStruct, 'perturbation')
%                 dataStruct.contrastLeft = dataStruct.contrastLeft(dataStruct.perturbation==0);
%                 dataStruct.contrastRight = dataStruct.contrastRight(dataStruct.perturbation==0);
%                 dataStruct.choice = dataStruct.choice(dataStruct.perturbation==0);
%                 dataStruct = rmfield(dataStruct,'perturbation');
%             end
%             
%             %Get subset of data which is compatible with detection
%             detectIdx = min(dataStruct.contrastLeft,dataStruct.contrastRight)==0;
%             dataStruct.contrastLeft = dataStruct.contrastLeft(detectIdx);
%             dataStruct.contrastRight = dataStruct.contrastRight(detectIdx);
%             dataStruct.choice = dataStruct.choice(detectIdx);
%             dataStruct.sessionID = dataStruct.sessionID(detectIdx);
%             dataStruct.subjectID = dataStruct.subjectID(detectIdx);
%             if any(~detectIdx)
%                 warning('DISCRIMINATION contrast conditions in dataset. Only plotting DETECTION subset of the data');
%             end
%             
%             detectIdx = min(dataStruct.testContrastRight,dataStruct.testContrastLeft)==0;
%             dataStruct.testContrastLeft = dataStruct.testContrastLeft(detectIdx);
%             dataStruct.testContrastRight = dataStruct.testContrastRight(detectIdx);
%             dataStruct.numTestContrasts = sum(detectIdx);
%             
%             %Plot posterior prediction
%             if ~isempty(Posterior)
%                 Posterior.interval = Posterior.interval(:,detectIdx,:);
%                 Posterior.mean = Posterior.mean(:,detectIdx,:);
%                 
%                 cDiff = dataStruct.testContrastRight - dataStruct.testContrastLeft;
%                 
%                 %Resort
%                 [cDiff, sortIdx]=sort(cDiff);
%                 Posterior.interval = Posterior.interval(:,sortIdx,:);
%                 Posterior.mean = Posterior.mean(:,sortIdx,:);
%                 
%                 if size(Posterior.interval,3)==1
%                     fx = fill(axis_handle,[cDiff; flipud(cDiff)], [Posterior.interval(2,:) fliplr( Posterior.interval(1,:) ) ], 'k');
%                     fx.EdgeAlpha=0;
%                     fx.FaceColor = [1 1 1]*0.8;
%                     
%                     plot(axis_handle,cDiff, Posterior.mean, '-');
%                     
%                 elseif size(Posterior.interval,3)==3
%                     
%                     for r = 1:3
%                         fx = fill(axis_handle,[cDiff; flipud(cDiff)], [Posterior.interval(2,:,r) fliplr( Posterior.interval(1,:,r) ) ], 'k');
%                         fx.FaceAlpha=0.3;
%                         fx.EdgeAlpha=0;
%                         fx.FaceColor = colours(r,:);
%                         plot(axis_handle,cDiff, Posterior.mean(:,:,r), '-', 'color', colours(r,:));
%                     end
%                 end
%                 
%             end
%             
%             %Plot data
%             cDiffData = dataStruct.contrastRight - dataStruct.contrastLeft;
%             cVals = unique(cDiffData);
%             for c = 1:length(cVals)
%                 resp = dataStruct.choice(cDiffData == cVals(c));
%                 [ph,pci] = binofit(sum([resp==1 resp==2 resp==3],1),length(resp));
%                 
%                 if max(dataStruct.choice)==2 %Forced choice
%                     %                     L=line(axis_handle,[1 1]*cVals(c),pci(2,:));
%                     %                     set(L,'Color',[1 1 1]*0,'Linewidth',0.5);
%                     plot(axis_handle,cVals(c),ph(2),'.','markersize',20,'color',[1 1 1]*0);
%                 elseif max(dataStruct.choice)==3 %Unforced choice
%                     
%                     for r = 1:3
%                         %                         L=line(axis_handle,[1 1]*cVals(c),pci(r,:));
%                         %                         set(L,'Color',colours(r,:),'Linewidth',0.5);
%                         plot(axis_handle,cVals(c),ph(r),'.','markersize',20,'color',colours(r,:));
%                     end
%                     
%                 end
%                 
%                 
%             end
%             
%             set(axis_handle,'xlim',[-1 1],'ylim',[0 1]);
%             ylabel(axis_handle,'pR');
%             xlabel(axis_handle,'CR - CL');
%             hold(axis_handle,'off');
%             %             set(axis_handle,'dataaspectratio',[1 1 1]);
%         end
