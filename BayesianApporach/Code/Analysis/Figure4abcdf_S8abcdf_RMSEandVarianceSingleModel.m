% Figure 4
species = {'Saccharomyces_cerevisiae','Yarrowia_lipolytica'};
currentpath = pwd;
cd ../../Results/model_bayesian
for k = 1:length(species)
    spec = species{k};
    cd(spec)
    nfound = length(dir('kcat_genra*.txt'));
    for i = 1:nfound
        tmp = readmatrix(['kcat_genra',num2str(i),'.txt'],'FileType','text','Delimiter',',');
        theta(i) = tmp(end,end);
    end
    x = [1:1:length(theta)];
    figure % plot RMSE
    hold on;
    fig1 = plot(x,theta,'--^','LineWidth',0.75,'MarkerEdgeColor',[33, 102, 172]/255,...
    'MarkerFaceColor',[33, 102, 172]/255,'MarkerSize',3);
  
    set(gca,'FontSize',6,'FontName','Helvetica');
    
    ylabel('RMSE in each generation','FontSize',7,'FontName','Helvetica','Color','k');
    xlabel('No. generation','FontSize',7,'FontName','Helvetica','Color','k');
    set(gcf,'position',[0 200 150 150]);
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    box off
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    saveas(fig1,['RMSE_DROP',species{k},'.pdf']);
    clear theta
    
    % plot varinence change
    color = [33, 102, 172;178, 24, 43]./255;
    num = [1,nfound];
    figure
    for j = 1:2
        hold on
        tmp = readmatrix(['kcat_genra',num2str(num(j)),'.txt'],'FileType','text','Delimiter',',');
        theta_100 = tmp(end,:); % is the rmse error
        kcat_100 = tmp(1:end-2,:);
        tot_prot_weight = tmp(end-1,1);
        % recalculate the sigma and mu
        ss = num2cell(kcat_100',1);
        [a,b] = arrayfun(@updateprior,ss);
        var_result{j} = b;
        fig2 = histogram(b,50,'FaceAlpha',0.3,'FaceColor',color(j,:),'EdgeColor',color(j,:));
    end
    set(gca,'FontSize',6,'FontName','Helvetica');
    ylim([0,1000])
    ylabel('Varaiance distribution','FontSize',7,'FontName','Helvetica','Color','k');
    xlabel('{\itk}_c_a_t variance in log10 scale','FontSize',7,'FontName','Helvetica','Color','k');
    leg = legend({'DL','Bayesian'});
    leg.ItemTokenSize = [10,2];
    set(gcf,'position',[0 200 150 150]);
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    set(leg,'box','off');
    box off
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    saveas(fig2,['VarianceDistribution',species{k},'.pdf']);
    cd ../
    
end


% identify rxns and enzymes are tuned toghther and plot
% species is the species list
% inputpath is the result data file char array
cd ../
for j = 1:length(species)
    cd ('model_bayesian/')
    cd(species{j})
    nfound = length(dir('kcat_genra*.txt'));
    tmp = readmatrix(['kcat_genra',num2str(nfound),'.txt'],'FileType','text','Delimiter',',');
    kcat_posterior = tmp(1:end-2,:);
    cd ../../model_dl
    z = load([species{j},'_dl.mat']);
    enzymedata = z.enzymedata;
    kcat_prior = arrayfun(@getrSample, enzymedata.kcat,enzymedata.kcat_var,enzymedata.enzyme_ec_kcat_range(:,1),enzymedata.enzyme_ec_kcat_range(:,2),repmat(100,length(enzymedata.kcat),1),'UniformOutput',false);
    kcat_prior = cell2mat(kcat_prior);
    [~,p]= ttest2(kcat_posterior,kcat_prior,'Dim',2,'Vartype','unequal','Alpha',0.01);
    p_adj = pval_adjust(p,'sidak');
    sig_enzyme{j} = enzymedata.rxn_list(p_adj < 0.01);
    sig_enzyme_count(j) = length(enzymedata.rxn_list(p_adj < 0.01));
    [~,p] = vartest2(kcat_posterior,kcat_prior,'Tail','left','Dim',2,'Alpha',0.01);
    p_adj = pval_adjust(p,'sidak');
    sig_sigma{j} =  enzymedata.rxn_list(p_adj < 0.01);
    sig_sigma_count(j) = length(enzymedata.rxn_list(p_adj < 0.01));
    cd ../model_bayesian
    cd(species{j})
    
    figure
    fig3 = bar([sig_enzyme_count(j),sig_sigma_count(j)],0.5,'FaceColor',[33, 102, 172]/255,'EdgeColor',[33, 102, 172]/255,'LineWidth',0.75);
    fig3.FaceAlpha = 0.3;
    set(gca,'FontSize',6,'FontName','Helvetica');
    xticklabels({'Mean value','Variance',7,'FontName','Helvetica','Color','k'})
    ylabel('Number of changed enzyme','FontSize',7,'FontName','Helvetica','Color','k');
    xtickangle(30)
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    set(gcf,'position',[0 200 150 150]);
    box off
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    saveas(fig3,['ChangedEnzyme',species{k},'.pdf']);
    cd('../../')
    
    
    % plot figure for correlation of kcat prior and kcat posterior
    figure
    hold on
    for i = 1:100
        fig4 = scatter(log10(enzymedata.kcat),log10(kcat_posterior(:,i)),20,'o','LineWidth',0.75,'MarkerFaceColor',[189,189,189]./255,'MarkerEdgeColor',[189,189,189]./255,'MarkerFaceAlpha',0.3);
    end
    ss = num2cell(kcat_posterior',1);
    [a,b] = arrayfun(@updateprior,ss);
    [RHOtmp,PVALtmp] = corr(log10(enzymedata.kcat),log10(a)','Type','Pearson');
    fig4 = scatter(log10(enzymedata.kcat),log10(a),20,'o','LineWidth',0.75,'MarkerFaceColor',[33, 102, 172]/255,'MarkerEdgeColor',[33, 102, 172]/255,'MarkerFaceAlpha',0.3);
    x = [0:0.001:9];
    fig4 = plot(x,x,'--k','LineWidth',0.75);
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('DL {\itk}_c_a_t in log10 scale','FontSize',7,'FontName','Helvetica','Color','k')
    ylabel('Bayesian {\itk}_c_a_t in log10 scale','FontSize',7,'FontName','Helvetica','Color','k');
    text(1,10.5,['p = ' num2str(round(PVALtmp,2))],'FontSize',7,'FontName','Helvetica','Color','k')
    text(1,11.5,['r = ' num2str(round(RHOtmp,2))],'FontSize',7,'FontName','Helvetica','Color','k')
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    set(gcf,'position',[0 200 150 150]);
    box off;
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    saveas(fig4,['KcatCorr',species{j},'.pdf']);
    cd ../
end
save('res_sig_prior_posterior.mat','sig_sigma_count','sig_sigma','sig_enzyme','sig_enzyme_count','species')

%% Figure for pca plot of kcat distribution
for j = 1:length(species)
    cd ('model_bayesian/')
    cd(species{j})
    load(['res_ForKcatPCA',species{j},'.mat']) % result from function PCAsampledKcatsOneSpecies in analysis
    figure
    hold on;
    [x,y] = sort(theta_all);
    h = scatter(score(:,1),score(:,2),20,'o','filled','LineWidth',0.75,'MarkerEdgeColor',[55,55,55]./255,'MarkerFaceColor',[55,55,55]./255,'MarkerFaceAlpha',0.3);legend off;
    h = scatter(score(y(end-100:end),1),score(y(end-100:end),2),20,'o','filled','LineWidth',0.75,'MarkerEdgeColor',[33, 102, 172]/255,'MarkerFaceColor',[33, 102, 172]/255,'MarkerFaceAlpha',0.3);legend off;
    h = scatter(score(y(1:100),1),score(y(1:100),2),40,'o','filled','LineWidth',0.75,'MarkerEdgeColor',[178, 24, 43]/255,'MarkerFaceColor',[178, 24, 43]/255,'MarkerFaceAlpha',0.3);legend off;
    leg = legend({'Intermediates','DL','Bayesian'});
    set(leg,'box','off');
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel(['PC1 (',num2str(round(explained(1),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
    ylabel(['PC2 (',num2str(round(explained(2),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    set(gcf,'position',[0 200 150 150]);
    box off;
    ax1 = gca;
    ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
        'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
        'XAxisLocation','top', 'XTick', [],...
        'YAxisLocation','right', 'YTick', []);
    linkaxes([ax1, ax2])
    saveas(h,['PCASampledKcatChange',species{j},'.pdf']);
    cd ../
end

