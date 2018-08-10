
%% Get data from omnibusGLM
% o1=omnibusLaserGLM('galvo_unilateral',{'Ochoa','Kornberg','Medawar','Burnet','Beadle','Bovet'},'session_date>"2018-03-18"');
o1=omnibusLaserGLM('sparse_unilateral',{'Spemann','Whipple','Morgan','Murphy','Chomsky'},'');
d3=o1.data{end};
% perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1','LeftM1','RightM1'};
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1'};

d3 = getrow(d3, any(d3.laserRegion == perturbationRegion,2) | d3.laserType==0);
dat = struct('contrastLeft', d3.stimulus(:,1),...
              'contrastRight', d3.stimulus(:,2),...
              'choice', d3.response,...
              'sessionID', d3.sessionID,...
              'subjectID', d3.subjectID);
dat.perturbation = zeros(size(dat.choice));
for p = 1:length(perturbationRegion)
    dat.perturbation(d3.laserRegion == perturbationRegion{p}) = p;
end

%% Gaussian process pulse model
o=omnibusLaserGLM('galvoPulse_unilateral',{'Ochoa','Kornberg','Medawar','Burnet','Beadle','Bovet'},'session_date>"2018-03-18"');
d3=o.data{end};

perturbationRegion = {'LeftVIS'};
d3 = getrow(d3, any(d3.laserRegion == perturbationRegion,2) | d3.laserType==0);
dat = struct('contrastLeft', d3.stimulus(:,1),...
              'contrastRight', d3.stimulus(:,2),...
              'choice', d3.response,...
              'sessionID', d3.sessionID,...
              'subjectID', d3.subjectID);
dat.perturbation = zeros(size(dat.choice));
for p = 1:length(perturbationRegion)
    dat.perturbation(d3.laserRegion == perturbationRegion{p}) = p;
end
dat.perturbationTime = d3.laserOnset;

%Sort trials by perturbation time, makes plotting easier later
dat.perturbationTime(dat.perturbation==0) = nan;
[~,sortIdx]=sort(dat.perturbationTime);
dat = structfun( @(f) f(sortIdx,:), dat, 'uni', 0);

fit = bfit.fitModel('Perturbation_GP',dat);



%% Model plot functions
%Parameters must be column vectors, contrast must be row vector
fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\2018-06-21_1703_Two-Level-Perturbation.mat");

%% Plot perturbation scatter

% m = mean([b.Posterior.bias_left_perturbation(:,p) b.Posterior.bias_right_perturbation(:,p)],1);
%         S = cov([b.Posterior.bias_left_perturbation(:,p) b.Posterior.bias_right_perturbation(:,p)]);
%         G = gmdistribution(m,S);
%         F = @(x,y) pdf(G,[x y]);
%
%         fax=ezcontour(F);
%         fax.LevelList=[0.2];
%         fax.LineColor = areaCols(p,:);
%         fax.LineWidth=2;
%
%

figure('color','w');
subplot(1,2,1); %Bias perturbation 
hold on;
% plot(fit.posterior.delta(:,1,1), fit.posterior.delta(:,2,1),'k.'); %Left VIS
m = mean([fit.posterior.delta(:,1,1), fit.posterior.delta(:,2,1)],1);
S = cov([fit.posterior.delta(:,1,1), fit.posterior.delta(:,2,1)]);
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fax=fcontour(F);
fax.LineColor = areaCols(1,:);
fax.LineWidth=1;


% plot(fit.posterior.delta(:,1,3), fit.posterior.delta(:,2,3),'k.'); %Left M2
m = mean([fit.posterior.delta(:,1,3), fit.posterior.delta(:,2,3)],1);
S = cov([fit.posterior.delta(:,1,3), fit.posterior.delta(:,2,3)]);
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fax=fcontour(F);
fax.LineColor = areaCols(3,:);
fax.LineWidth=1;

line([0 0],[-4 4]);
line([-4 4],[0 0]);
set(gca,'dataaspectratio',[1 1 1]);
xlabel('BL'); ylabel('BR');

subplot(1,2,2); %Bias perturbation 
hold on;
% plot(fit.posterior.delta(:,1,1), fit.posterior.delta(:,2,1),'k.'); %Left VIS
m = mean([fit.posterior.delta(:,3,1), fit.posterior.delta(:,4,1)],1);
S = cov([fit.posterior.delta(:,3,1), fit.posterior.delta(:,4,1)]);
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fax=fcontour(F);
fax.LineColor = areaCols(1,:);
fax.LineWidth=1;


% plot(fit.posterior.delta(:,1,3), fit.posterior.delta(:,2,3),'k.'); %Left M2
m = mean([fit.posterior.delta(:,3,3), fit.posterior.delta(:,4,3)],1);
S = cov([fit.posterior.delta(:,3,3), fit.posterior.delta(:,4,3)]);
G = gmdistribution(m,S);
F = @(x,y) pdf(G,[x y]);
fax=fcontour(F);
fax.LineColor = areaCols(3,:);
fax.LineWidth=1;

line([0 0],[-4 4]);
line([-4 4],[0 0]);
set(gca,'dataaspectratio',[1 1 1]);
xlabel('SL'); ylabel('SR');




%% Detection plots for perturbation conditions
CL = [linspace(1,0.1,200), linspace(0.1,0,400), zeros(1,600)];
CR = [zeros(1,600), linspace(0,0.1,400) linspace(0.1,1,200)];

%Non-laser global fit on detection trials
BL = fit.posterior.bias(:,1) ;
BR = fit.posterior.bias(:,2);
SL = fit.posterior.sens(:,1);
SR = fit.posterior.sens(:,2);
N = fit.posterior.n_exp;
NL_ph=PHAT(BL,SL,BR,SR,N,CL,CR);
NL_phM=mean(NL_ph,1);
NL_phQ=quantile(NL_ph,[0.025 0.975],1);

%non-laser detection trial data
idx=min([fit.data_stan.contrastLeft fit.data_stan.contrastRight],[],2)==0;
cDiff = fit.data_stan.contrastRight(idx) - fit.data_stan.contrastLeft(idx);
resp = fit.data_stan.choice(idx);
perturb = fit.data_stan.perturbation(idx);
sess = fit.data_stan.sessionID(idx);
[counts,~,~,lab] = crosstab(cDiff,resp,perturb,sess);
prob = counts./sum(counts,2);%Convert to probability over choices
prob = nanmean(prob,4); %Average over sessions
cDiffUnique=cellfun(@str2num,lab(1:size(prob,1),1));

figure;
ha = tight_subplot(8,3,0.01,0.05,0.05);

figure;
ha2 = tight_subplot(8,3,0.01,0.05,0.05);

% for region = 1:6
for region = 1:8
    %Global fit
    BL = fit.posterior.bias(:,1) + fit.posterior.delta(:,1,region);
    BR = fit.posterior.bias(:,2) + fit.posterior.delta(:,2,region);
    SL = fit.posterior.sens(:,1) + fit.posterior.delta(:,3,region);
    SR = fit.posterior.sens(:,2) + fit.posterior.delta(:,4,region);
    N = fit.posterior.n_exp;
    
    ph=PHAT(BL,SL,BR,SR,N,CL,CR);
    phM=mean(ph,1);
    phQ=quantile(ph,[0.025 0.975],1);
    
    phDiff = ph-NL_ph;
    phDiffM=mean(phDiff,1);
    phDiffQ=quantile(phDiff,[0.025 0.975],1);
    
    for r = 1:3
        hold( ha(3*(region-1) + r), 'on');
        fx = fill(ha(3*(region-1) + r),[CR-CL fliplr(CR-CL)], [NL_phQ(1,:,r) fliplr( NL_phQ(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0;
        plot(ha(3*(region-1) + r), CR-CL,NL_phM(1,:,r),'k-');
        
        fx = fill(ha(3*(region-1) + r),[CR-CL fliplr(CR-CL)], [phQ(1,:,r) fliplr( phQ(2,:,r) ) ], 'r');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0;
        plot(ha(3*(region-1) + r),CR-CL,phM(1,:,r),'r-');
        
        plot(ha(3*(region-1) + r),cDiffUnique,prob(:,r,1),'k.','markersize',10);
        plot(ha(3*(region-1) + r),cDiffUnique,prob(:,r,1+region),'r.','markersize',10);

        
        hold( ha2(3*(region-1) + r), 'on');
        line(ha2(3*(region-1) + r),[-1 1],[0 0],'color',[1 1 1]*0.5,'linestyle','-');
        fx = fill(ha2(3*(region-1) + r),[CR-CL fliplr(CR-CL)], [phDiffQ(1,:,r) fliplr( phDiffQ(2,:,r) ) ], 'r');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0;
        plot(ha2(3*(region-1) + r),CR-CL,phDiffM(1,:,r),'r-');
        plot(ha2(3*(region-1) + r),cDiffUnique,prob(:,r,1+region) - prob(:,r,1),'r.','markersize',10);
    end
end
set([ha;ha2],'xcolor','none','ycolor','none');
set([ha(end-2);ha2(end-2)],'ycolor','k','xcolor','k','xticklabelmode','auto','yticklabelmode','auto')

set(ha2,'ylim',[-1 1]*0.5,'xlim',[-1 1]);
set(ha,'ylim',[0 1],'xlim',[-1 1]);
title(ha(1),'pLeft');
title(ha(2),'pRight');
title(ha(3),'pNoGo');
title(ha2(1),'\Delta pLeft');
title(ha2(2),'\Delta pRight');
title(ha2(3),'\Delta pNoGo');

%% Plot model surfaces
numC = 50;
Cs = linspace(0,1,numC);
[CR,CL] = meshgrid(Cs); 
CL = CL(:)';
CR = CR(:)';


%Non-laser global fit on detection trials
BL = fit.posterior.bias(:,1);
BR = fit.posterior.bias(:,2);
SL = fit.posterior.sens(:,1);
SR = fit.posterior.sens(:,2);
N = fit.posterior.n_exp;
NL_ph = PHAT(BL,SL,BR,SR,N,CL,CR);
NL_matrix = reshape(  permute(mean(NL_ph,1),[2 3 1]) ,numC,numC,3);


figure;
ha = tight_subplot(5,3,0.01,0.05,0.05);

for r = 1:3
    
    imagesc(ha(r),Cs,Cs,NL_matrix(:,:,r));
    caxis(ha(r),[0 1]);
end

for region = 1:4
    %Global fit
    BL = fit.posterior.bias(:,1) + fit.posterior.delta(:,1,region);
    BR = fit.posterior.bias(:,2) + fit.posterior.delta(:,2,region);
    SL = fit.posterior.sens(:,1) + fit.posterior.delta(:,3,region);
    SR = fit.posterior.sens(:,2) + fit.posterior.delta(:,4,region);
    N = fit.posterior.n_exp;
    
    ph = PHAT(BL,SL,BR,SR,N,CL,CR);
    phDiff = ph-NL_ph;
    
    phDiff_matrix = reshape(  permute(mean(phDiff,1),[2 3 1]) ,numC,numC,3);
    
    for r = 1:3
        
        imagesc(ha(3*(region) + r),Cs,Cs,phDiff_matrix(:,:,r));
        caxis(ha(3*(region) + r),[-1 1]*0.5);
    end
end
set(ha,'ydir','normal','dataaspectratio',[1 1 1]);
                                cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                                colormap(flipud(cmap));

%% Plot perturbation scatter
% param_axes = {'bias_left_perturbation', 'bias_right_perturbation';
%                'sens_left_perturbation', 'sens_right_perturbation';
%                'bias_left_perturbation','sens_left_perturbation';
%                'bias_right_perturbation','sens_right_perturbation'};
           
% m = mean([b.Posterior.bias_left_perturbation(:,p) b.Posterior.bias_right_perturbation(:,p)],1);
%         S = cov([b.Posterior.bias_left_perturbation(:,p) b.Posterior.bias_right_perturbation(:,p)]);
%         G = gmdistribution(m,S);
%         F = @(x,y) pdf(G,[x y]);
% 
%         fax=ezcontour(F);
%         fax.LevelList=[0.2];
%         fax.LineColor = areaCols(p,:);
%         fax.LineWidth=2;
%         
% 
figure('color','w');

subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;
px=gobjects;





for p = 1:length(perturbationRegion)
    subplot(1,2,1);
    
    xx=squeeze(b.Posterior.delta(:,1,p));
    yy=squeeze(b.Posterior.delta(:,2,p));
    px(p)=scatter(xx, yy, 'filled','SizeData',10);
    m = mean([xx yy],1);
    S = cov([xx yy]);
    G = gmdistribution(m,S);
    F = @(x,y) pdf(G,[x y]);
    %
    fax=ezcontour(F);
    fax.LevelList=[0.5];
    fax.LineColor = px(p).CData;
    fax.LineWidth=2;

    subplot(1,2,2);
    xx=squeeze(b.Posterior.delta(:,3,p));
    yy=squeeze(b.Posterior.delta(:,4,p));
    pxx=scatter(xx, yy, 'filled','SizeData',10);
    m = mean([xx yy],1);
    S = cov([xx yy]);
    G = gmdistribution(m,S);
    F = @(x,y) pdf(G,[x y]);
    %
    fax=ezcontour(F);
    fax.LevelList=[0.2];
    fax.LineColor = pxx.CData;
    fax.LineWidth=2;

end

subplot(1,2,1);
axis equal;
legend(px,perturbationRegion);
xlabel('BL perturbation');
ylabel('BR perturbation');
xlim auto;
ylim auto;
line([0 0],ylim(gca)); line(xlim(gca),[0 0]);

subplot(1,2,2);
axis equal;
xlabel('SL perturbation');
ylabel('SR perturbation');
xlim auto;
ylim auto;
line([0 0],ylim(gca)); line(xlim(gca),[0 0]);


%Draw nonlaser psych curves
credible_interval = [0.025 0.975];
psNL = struct;
psNL.interval = quantile(b.Posterior.pTest,credible_interval,1);
psNL.mean = mean(b.Posterior.pTest,1);

idx = b.data_stan.perturbation==0;
data_subset_NL = struct('testContrastLeft',b.data_stan.testContrastLeft,...
    'testContrastRight',b.data_stan.testContrastRight,...
    'contrastLeft',b.data.contrastLeft(idx),...
    'contrastRight',b.data.contrastRight(idx),...
    'choice',b.data.choice(idx));

for p = 1:length(perturbationRegion)
    ps = struct;
    ps.interval = quantile( permute(b.Posterior.pTest_perturbation(:,:,p,:),[1 2 4 3] ),credible_interval,1);
    ps.mean = mean(permute(b.Posterior.pTest_perturbation(:,:,p,:),[1 2 4 3] ),1);

    idx = b.data_stan.perturbation==p;
    data_subset = struct('testContrastLeft',b.data_stan.testContrastLeft,...
        'testContrastRight',b.data_stan.testContrastRight,...
        'contrastLeft',b.data.contrastLeft(idx),...
        'contrastRight',b.data.contrastRight(idx),...
        'choice',b.data.choice(idx));
    
    f2=figure;
    axL=axes;
    b.util_plotDiscrimination(axL, data_subset, ps);
    set(findobj(f2,'Type','Line'),'Color',[1 0 0]);
    set(findobj(f2,'Type','Patch'),'FaceColor',[1 0 0]);
    
    f1=figure;
    axNL=axes;
    b.util_plotDiscrimination(axNL, data_subset_NL, psNL);
    set(findobj(f1,'Type','Line'),'Color',[0 0 0]);
    set(findobj(f1,'Type','Patch'),'FaceColor',[0 0 0]);
    
    axs1 = get(f1,'children');
    axs2 = get(f2,'children');
    
    for ax = 1:length(axs1)
        nLstuff = get(axs1(ax),'children');
        Lstuff = get(axs2(ax),'children');
        copyobj(Lstuff, axs1(ax))
    end
    
    close(f2);
    set(f1,'name',perturbationRegion{p});
end

%% Plot hierarchical correlation matrices
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
       linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
   
pNames = fieldnames(b.Posterior);
corrParams = pNames(contains(fieldnames(b.Posterior),'corr')); %Get corr parameters

figure('color','w','name','correlation of parameter deviations in hierarchy');
for c = 1:length(corrParams)
    
    thisC = b.Posterior.(corrParams{c});
    
    subplot(1,length(corrParams),c);
    hx=imagesc(permute(mean(thisC,1),[2 3 1])); caxis([-1 1]); title('Correlation over sessions');

    title(corrParams{c},'interpreter','none');
    
    xlabel('Parameter ID'); ylabel('Parameter ID');
    set(gca,'xtick',1:5,'ytick',1:5,'dataaspectratio',[1 1 1]);
    colorbar;
end
colormap(flipud(cmap));

%% Scatter posterior deltas
% areaCols = [      0    0.4470    0.7410; %LeftVIS
%         0.3010    0.7450    0.9330; %RightVIS
%         0.4660    0.6740    0.1880; %LeftM2
%     0.8500    0.3250    0.0980; %RightM2
%     0.4940    0.1840    0.5560; %purple
%     0.8 0.8 0.8]; %grey

areaCols = [94 60 153;
            178 171 210;
            230 97 1;
            253 184 99]/255;%LeftVIS

figure('color','w');
subplot(1,2,1); hold on;
fax2=gobjects;
for p = 1:4

    for subj = 1:b.data_stan.numSubjects

        B_perturb = [b.Posterior.BL_perturbation(:,p) + b.Posterior.BL_perturbation_subj(:,subj,p),...
            b.Posterior.BR_perturbation(:,p) + b.Posterior.BR_perturbation_subj(:,subj,p)];
        
        m = mean(B_perturb,1);
        S = cov(B_perturb);
        G = gmdistribution(m,S);
        F = @(x,y) pdf(G,[x y]);
        fax=ezcontour(F);
        fax.LevelList=[0.2];
        fax.LineWidth=1;
        fax.LineColor=areaCols(p,:);
    end

    B_perturb = [b.Posterior.BL_perturbation(:,p),...
        b.Posterior.BR_perturbation(:,p)];
    m = mean(B_perturb,1);
    S = cov(B_perturb);
    G = gmdistribution(m,S);
    F = @(x,y) pdf(G,[x y]);
    fax2(p)=ezcontour(F);
    fax2(p).LevelList=[0.2];
    fax2(p).LineWidth=5;
    fax2(p).LineColor=areaCols(p,:);
end

axis auto;
line([0 0],ylim(gca),'Color',[0 0 0],'linestyle','--'); line(xlim(gca),[0 0],'Color',[0 0 0],'linestyle','--');
title('');
xlabel('BL perturbation'); ylabel('BR perturbation');
legend(fax2,perturbationRegion(1:4));


subplot(1,2,2); hold on; 
fax2=gobjects;
for p = 1:4

    for subj = 1:b.data_stan.numSubjects

        S_perturb = [b.Posterior.SL_perturbation(:,p) + b.Posterior.SL_perturbation_subj(:,subj,p),...
            b.Posterior.SR_perturbation(:,p) + b.Posterior.SR_perturbation_subj(:,subj,p)];
        
        m = mean(S_perturb,1);
        S = cov(S_perturb);
        G = gmdistribution(m,S);
        F = @(x,y) pdf(G,[x y]);
        fax=ezcontour(F);
        fax.LevelList=[0.2];
        fax.LineWidth=1;
        fax.LineColor=areaCols(p,:);
    end

    S_perturb = [b.Posterior.SL_perturbation(:,p),...
        b.Posterior.SR_perturbation(:,p)];
    m = mean(S_perturb,1);
    S = cov(S_perturb);
    G = gmdistribution(m,S);
    F = @(x,y) pdf(G,[x y]);
    fax2(p)=ezcontour(F);
    fax2(p).LevelList=[0.2];
    fax2(p).LineWidth=5;
    fax2(p).LineColor=areaCols(p,:);
end

axis auto;
line([0 0],ylim(gca),'Color',[0 0 0],'linestyle','--'); line(xlim(gca),[0 0],'Color',[0 0 0],'linestyle','--');
title('');
xlabel('SL perturbation'); ylabel('SR perturbation');
legend(fax2,perturbationRegion(1:4));
