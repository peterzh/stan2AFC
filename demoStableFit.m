clear all; close all;

%Demo more stable fitting
load('2AFC_testData_Session_oneSubject.mat');

%Get data from one session
sessID = 1;
data_sess = structfun(@(f) f(data2.sessionID==sessID), data2, 'uni', 0);

%Many-session data with the sessID data removed, to be replaced with
%subsets from that session
data_notsess = structfun(@(f) f(data2.sessionID~=sessID), data2, 'uni', 0);

sampleSize = round(linspace(10,150,10));
for p = 1:length(sampleSize)
    
    trialIdx = randsample( length(data_sess.choiceR), sampleSize(p));    
    %Portion of data from that session
    data_sess_part = structfun(@(f) f(trialIdx), data_sess, 'uni', 0);
    
    %Portion of data from that session + all other sessions' data
    data_sess_hier_part = table2struct([struct2table(data_notsess); struct2table(data_sess_part)],'ToScalar',true);
    
    %For plotting...

    %Fit single session on its own
    b = behavModel(data_sess_part,'bias_sens');
    b.getPosterior;
%     b.getMAP;
    ha=subplot(2,length(sampleSize),p);

    data_sess_part.testContrastLeft = b.data_stan.testContrastLeft;
    data_sess_part.testContrastRight = b.data_stan.testContrastRight;
    
%     b.util_plotSingle(ha, data_sess_part, [], b.MAP.pRTest);
    b.util_plotSingle(ha, data_sess_part, quantile(b.Posterior.pRTest,[0.025 0.975]), []);
%     b.util_plotSingle(ha, data_sess_part, quantile(b.Posterior.pRTest,[0.025 0.975]), b.MAP.pRTest);

    title(ha,sprintf('%d trials', sampleSize(p))); drawnow;
    
    %Fit session within a larger hierarchical framework
    b = behavModel(data_sess_hier_part,'bias1_sens');
    b.getPosterior;
%     b.getMAP;
    ha=subplot(2,length(sampleSize),p+length(sampleSize));
%     b.util_plotSingle(ha, data_sess_part, [], b.MAP.pRTest(sessID,:));
    b.util_plotSingle(ha, data_sess_part, quantile(squeeze(b.Posterior.pRTest(:,sessID,:)),[0.025 0.975]), []);
%     b.util_plotSingle(ha, data_sess_part, quantile(squeeze(b.Posterior.pRTest(:,sessID,:)),[0.025 0.975]), b.MAP.pRTest(sessID,:));
    drawnow;
end





