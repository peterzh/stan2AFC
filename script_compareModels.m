
models = {'bias1_sens1','bias1_sens','bias1_sensL_sensR_Exponent','bias1_sensL1_sensR1_Exponent','bias1_sensL1_sensR1_Exponent1'};

figure;
axis; hold on;
data
for m = 1:length(models)
    b = behavModel(data2,models{m});
    post = b.getPosterior;

%     b.plotPsych; drawnow;
    [total,se,pointwise] = mstan.waic(post.log_lik);
    
    lx = line( [1 1]*m, [total.waic-se.waic total.waic+se.waic]);
    plot( m, total.waic, 'ko'); drawnow;
end

set(gca,'xtick',1:length(models),'xticklabels',models);
xlim([0 length(models)+1]);