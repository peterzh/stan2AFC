function plot2AUCDiscrimination(posterior,data)
colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];

%Compute data probabilities
if ~isfield(data,'perturbation')
    data.perturbation = zeros(size(data.response));
end

%Transform to pedestal/diff representation
data.pedestal = min(data.contrastLeft,data.contrastRight);
data.cDiff = data.contrastRight - data.contrastLeft;

[counts,~,~,labels] = crosstab(data.cDiff,...
    data.pedestal,...
    data.choice,...
    data.perturbation,...
    data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob_ave = nanmean(prob,5);
cdiffs = cellfun(@str2num,labels(1:size(prob,1),1));
peds = cellfun(@str2num,labels(1:size(prob,2),2));

for subj = 1:data.numSubjects
    
    sessionsThisSubj = find(data.subjID_session == subj);
    
    f=figure('color','w');
    ha_background = tight_subplot( ceil(sqrt( length(sessionsThisSubj) )), ceil(sqrt( length(sessionsThisSubj) )) , [0.04 0.01], [0.01 0.04], 0.01 );
    
    for sid = 1:length(sessionsThisSubj)        
        sessNum = sessionsThisSubj(sid);
        
        %Parameters for this session
        BL = posterior.bias(:,1) + posterior.b_sess(:,1,sessNum) + posterior.b_subj(:,1,subj);
        BR = posterior.bias(:,2) + posterior.b_sess(:,2,sessNum) + posterior.b_subj(:,2,subj);
        SL = posterior.sens(:,1) + posterior.b_sess(:,3,sessNum) + posterior.b_subj(:,3,subj);
        SR = posterior.sens(:,2) + posterior.b_sess(:,4,sessNum) + posterior.b_subj(:,4,subj);
        N = posterior.n_exp      + posterior.b_sess(:,5,sessNum) + posterior.b_subj(:,5,subj);
        
        %For each session, add new subplots
        set(ha_background(sid),'units','normalized');
        bounds = get(ha_background(sid),'position');
        set(ha_background(sid),'xcolor','none','ycolor','none');
        
        title(ha_background(sid),sprintf('M%d S%d',subj,sessNum));
        
        %replace with many subplots within the bounds of the original
        ha = tight_subplot(4,3,0.001, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
        set(ha(1:9),'xcolor','none');
        set(ha([2 3 5 6 8 9 11 12]),'ycolor','none');
        for i = 1:length(ha)
            hold(ha(i),'on');
        end
        
        for ped = 1:length(peds)
            CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
            CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
            ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
            ph_interval = quantile(ph,[0.025 0.975]);
            ph_mean = mean(ph,1);
            
            for r = 1:3
                %plot datapoints
                plot(ha( 3*(ped-1) + r ),cdiffs,prob(:,ped,r,1,sessNum),'.','markersize',10,'color',colours(r,:));
                
                fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
                fx.FaceAlpha=0.3;
                fx.EdgeAlpha=0;
                fx.FaceColor = colours(r,:);
                plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
            end
            
        end
        
    end
    
    set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);

end
% 
% 
% f=figure('color','w');
% ha_background = tight_subplot( ceil(sqrt(max(data.sessionID))), ceil(sqrt(max(data.sessionID))), [0.04 0.01], [0.01 0.04], 0.01 );
% Plot curves for each session
% for session = 1:max(data.sessionID)
%     Parameters for this session
%     BL = posterior.bias(:,1) + posterior.b_sess(:,1,session) + posterior.b_subj(:,1,data.subjID_session(session));
%     BR = posterior.bias(:,2) + posterior.b_sess(:,2,session) + posterior.b_subj(:,2,data.subjID_session(session));
%     SL = posterior.sens(:,1) + posterior.b_sess(:,3,session) + posterior.b_subj(:,3,data.subjID_session(session));
%     SR = posterior.sens(:,2) + posterior.b_sess(:,4,session) + posterior.b_subj(:,4,data.subjID_session(session));
%     N = posterior.n_exp      + posterior.b_sess(:,5,session) + posterior.b_subj(:,5,data.subjID_session(session));
%     
%     For each session, add new subplots
%     set(ha_background(session),'units','normalized');
%     bounds = get(ha_background(session),'position');
%     set(ha_background(session),'xcolor','none','ycolor','none');
%     
%     title(ha_background(session),sprintf('M%d S%d',data.subjID_session(session),session));
%     
%     replace with many subplots within the bounds of the original
%     ha = tight_subplot(4,3,0.001, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
%     set(ha(1:9),'xcolor','none');
%     set(ha([2 3 5 6 8 9 11 12]),'ycolor','none');
%     for i = 1:length(ha)
%         hold(ha(i),'on');
%     end
%     
%     for ped = 1:length(peds)
%         CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
%         CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
%         ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
%         ph_interval = quantile(ph,[0.025 0.975]);
%         ph_mean = mean(ph,1);
%         
%         for r = 1:3
%             plot datapoints
%             plot(ha( 3*(ped-1) + r ),cdiffs,prob(:,ped,r,1,session),'.','markersize',10,'color',colours(r,:));
%             
%             fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
%             fx.FaceAlpha=0.3;
%             fx.EdgeAlpha=0;
%             fx.FaceColor = colours(r,:);
%             plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
%         end
%         
%     end
%     
% end

% set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);


%Plot grand average

%Parameters of the 'average session'
BL = posterior.bias(:,1) ;
BR = posterior.bias(:,2) ;
SL = posterior.sens(:,1) ;
SR = posterior.sens(:,2) ;
N = posterior.n_exp ;

%Parameters drawn from the variation across subjects and sessions
% BL = posterior.bias(:,1) + posterior.sd_subj(:,1).*randn(2000,1) + posterior.sd_sess(:,1).*randn(2000,1);
% BR = posterior.bias(:,2) + posterior.sd_subj(:,2).*randn(2000,1) + posterior.sd_sess(:,2).*randn(2000,1);
% SL = posterior.sens(:,1) + posterior.sd_subj(:,3).*randn(2000,1) + posterior.sd_sess(:,3).*randn(2000,1);
% SR = posterior.sens(:,2) + posterior.sd_subj(:,4).*randn(2000,1) + posterior.sd_sess(:,4).*randn(2000,1);
% N = posterior.n_exp      + posterior.sd_subj(:,5).*randn(2000,1) + posterior.sd_sess(:,5).*randn(2000,1);

f=figure;
ha = tight_subplot(4,3,0.001, 0.01,0.01);
for i = 1:length(ha)
    hold(ha(i),'on');
end

numSess = size(prob,5);
sessJitter = 0.01*randn(numSess,1);
for ped = 1:length(peds)
    CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
    CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
    ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
    ph_interval = quantile(ph,[0.025 0.975]);
    ph_mean = nanmean(ph,1);
    
    for r = 1:3
        %plot datapoints
        plot(ha( 3*(ped-1) + r ),cdiffs,prob_ave(:,ped,r,1),'.','markersize',10,'color',colours(r,:));
        plot(ha( 3*(ped-1) + r ),cdiffs + sessJitter',permute( prob(:,ped,r,1,:), [1 5 2 3 4]) ,'.','markersize',1,'color',colours(r,:));
        
        
        
        fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3;
        fx.EdgeAlpha=0;
        fx.FaceColor=colours(r,:);
        plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
    end
    
end
set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);



%% Try out the other kind of curves
[counts,~,~,labels] = crosstab(data.contrastLeft,...
    data.contrastRight,...
    data.choice,...
    data.perturbation,...
    data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob = prob(:,:,:,1,:); %only laser off conditon
prob_ave = nanmean(prob,5);


cls = cellfun(@str2num,labels(1:size(prob,1),1));
crs = cellfun(@str2num,labels(1:size(prob,2),2));

colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];


figure;
ha = tight_subplot(2,2,0.05,0.05,0.05);
for i = 1:length(ha)
    hold(ha(i),'on');
end
set(ha([1 3]),'xdir','reverse','yaxislocation','right');

%Varying CL while holding CR = 0
CL = linspace(0,1,100);
CR = zeros(1,100);
ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
ph_interval = quantile(ph,[0.025 0.975]);
ph_mean = mean(ph,1);
for r=1:3
    plot(ha(1),CL, ph_mean(1,:,r), '-','Color',colours(r,:));
    fx = fill(ha(1),[CL fliplr(CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
    fx.FaceAlpha=0.3;
    fx.EdgeAlpha=0;
    fx.FaceColor = colours(r,:);
    
    plot(ha(1),cls,prob_ave(:,1,r),'.','markersize',15,'color',colours(r,:));
    plot(ha(1),cls,squeeze(prob(:,1,r,1,:)),'.','markersize',1,'color',colours(r,:));
end

%Varying CR while holding CL = 0
CR = linspace(0,1,100);
CL = zeros(1,100);
ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
ph_interval = quantile(ph,[0.025 0.975]);
ph_mean = mean(ph,1);
for r=1:3
    plot(ha(2),CR, ph_mean(1,:,r), '-','Color',colours(r,:));
    fx = fill(ha(2),[CR fliplr(CR)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
    fx.FaceAlpha=0.3;
    fx.EdgeAlpha=0;
    fx.FaceColor = colours(r,:);
    
    plot(ha(2),cls,prob_ave(1,:,r),'.','markersize',15,'color',colours(r,:));
    plot(ha(2),cls,squeeze(prob(1,:,r,1,:)),'.','markersize',1,'color',colours(r,:));
end

%Varying CL while holding CR at different values, only plotting one side
for cr = 1:length(crs)
    CL = linspace(0,1,100);
    CR = ones(1,100)*crs(cr);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_mean = mean(ph,1);
    
    r=1;
    plot(ha(3),CL, ph_mean(1,:,r), '-','Color',[colours(r,:) cr/length(crs)],'linewidth',1);
    
    px=scatter(ha(3), cls,prob_ave(:,cr,r),'filled','SizeData',15);
    px.MarkerFaceColor = colours(r,:);
    alpha(px,cr/length(crs))
end

%Varying CR while holding CL at different values, only plotting one side
for cl = 1:length(cls)
    CR = linspace(0,1,100);
    CL = ones(1,100)*cls(cl);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_mean = mean(ph,1);
    
    r=2;
    plot(ha(4),CR, ph_mean(1,:,r), '-','Color',[colours(r,:) cl/length(cls)],'linewidth',1);
        
    px=scatter(ha(4), crs,prob_ave(cl,:,r),'filled','SizeData',15);
    px.MarkerFaceColor = colours(r,:);
    alpha(px,cl/length(cls))
end
set(ha,'xlim',[0 1],'ylim',[0 1],'dataaspectratio',[1 1 1],'xticklabelmode','auto')

xlabel(ha(3),'ContrastLeft');
xlabel(ha(4),'ContrastRight');
ylabel(ha(3),'pChoice');

title(ha(1),'CR = 0');
title(ha(2),'CL = 0');

title(ha(3),'CR = [0, 0.1, 0.24, 0.54]');
title(ha(4),'CL = [0, 0.1, 0.24, 0.54]');

set(ha,'xtick',cls);

%% Plot per-session variation in nogo (c=0) behaviour
pNG = squeeze(prob(1,1,3,1,:));
subjID = data.subjID_session;

[~,sortIdx]=sortrows([subjID pNG],[1 2]);
% [~,sortIdx] = sort(pNG);

subjCols = parula(max(subjID));

f=figure('color','w');
axes; hold on;
for i = 1:length(sortIdx)
    thisSessIdx = sortIdx(i);
    BL = posterior.bias(:,1) + posterior.b_sess(:,1,thisSessIdx) + posterior.b_subj(:,1,data.subjID_session(thisSessIdx));
    BR = posterior.bias(:,2) + posterior.b_sess(:,2,thisSessIdx) + posterior.b_subj(:,2,data.subjID_session(thisSessIdx));
    SL = posterior.sens(:,1) + posterior.b_sess(:,3,thisSessIdx) + posterior.b_subj(:,3,data.subjID_session(thisSessIdx));
    SR = posterior.sens(:,2) + posterior.b_sess(:,4,thisSessIdx) + posterior.b_subj(:,4,data.subjID_session(thisSessIdx));
    N = posterior.n_exp      + posterior.b_sess(:,5,thisSessIdx) + posterior.b_subj(:,5,data.subjID_session(thisSessIdx));
    
    plot(i, prob(1,1,3,1,thisSessIdx), '.', 'markersize', 20, 'color', subjCols(subjID(thisSessIdx),:));
    
    %model prediction
    ph=bplot.CN(BL,SL,BR,SR,N, 0 , 0 );
    pM = mean(ph(:,1,3),1);
    pQ = quantile(ph(:,1,3),[0.025 0.975],1);
    
    fx=fill(i + [-1 1 1 -1]*0.45, pM + [-1 -1 1 1]*0.001, 'k' ); fx.EdgeAlpha=0;
    fx.FaceColor = subjCols(subjID(thisSessIdx),:);
    
    
    fx=fill(i + [-1 1 1 -1]*0.45, [pQ(1) pQ(1) pQ(2) pQ(2)] , 'k' ); alpha(fx,0.1); fx.EdgeAlpha=0.1;
    fx.FaceColor = subjCols(subjID(thisSessIdx),:);
end

ylabel('pNG on zero contrast');
xlabel('Session');
ylim([0 1]);

%% Overlay per-session curves for NoGos at different detection contrasts
f=figure;
ha = tight_subplot(1,1,0.001, 0.01,0.01);
for i = 1:length(ha)
    hold(ha(i),'on');
end

CL = [linspace(1,0,100) zeros(1,100)];
CR = [zeros(1,100) linspace(0,1,100)];

for sess = 1:max(data.sessionID)
    BL = posterior.bias(:,1) + posterior.b_sess(:,1,sess) + posterior.b_subj(:,1,data.subjID_session(sess));
    BR = posterior.bias(:,2) + posterior.b_sess(:,2,sess) + posterior.b_subj(:,2,data.subjID_session(sess));
    SL = posterior.sens(:,1) + posterior.b_sess(:,3,sess) + posterior.b_subj(:,3,data.subjID_session(sess));
    SR = posterior.sens(:,2) + posterior.b_sess(:,4,sess) + posterior.b_subj(:,4,data.subjID_session(sess));
    N = posterior.n_exp      + posterior.b_sess(:,5,sess) + posterior.b_subj(:,5,data.subjID_session(sess));
    
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    pNG = mean(ph(:,:,3),1);
    
    plot( CR-CL, pNG, 'k-');
end


end