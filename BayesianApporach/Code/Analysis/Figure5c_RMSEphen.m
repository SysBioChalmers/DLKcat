% Figure 4c
load('species_withdata.mat')
currentpath = pwd;
cd ../../Results/model_build_files/
for i = 1:length(species_withdata)
    cd('model_classic')
    i
    z = load([species_withdata{i},'_classic.mat']);
    model = z.model;
    enzymedata = z.enzymedata;
    strain = z.strain;
    growthdata = z.growthdata;
    max_growth = z.max_growth;
    rxn2block = z.rxn2block;
    [~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
    tot_prot_weight = tot_prot_weight*0.5;
    
    [rmse_final2,~,~,~,~] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
    result(i,1) = rmse_final2;
    coverage(i,1) = length(enzymedata.kcat);
    cd ../model_dl
    z = load([species_withdata{i},'_dl.mat']);
    model = z.model;
    enzymedata = z.enzymedata;
    
    [rmse_final2,exp,simulated3,growthdata,max_growth] = abc_matlab_max(model,enzymedata,enzymedata.kcat,tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block);
    result(i,2) = rmse_final2;
    coverage(i,2) = length(enzymedata.kcat);
    cd ../model_bayesian
    cd(species_withdata{i})
    a = dir('kcat_genra*.txt');
    if ~isempty(a)
        a = {a.name};
        tmp = strrep(a,'kcat_genra','');
        tmp = strrep(tmp,'.txt','');
        tmp = cellfun(@str2double,tmp);
        tmp = max(tmp);
        a = ['kcat_genra',num2str(tmp),'.txt'];
        tmp = readmatrix(a,'FileType','text','Delimiter',',');
        tot_prot_weight = tmp(end-1,1);
        theta_100 = tmp(end,:); % is the rmse error
        result(i,3) = mean(theta_100);
        coverage(i,3) = coverage(i,2);
    else
        result(i,3) = result(i,2);
        coverage(i,3) = coverage(i,2);
    end
    
    cd ../../
end
save('res_RMSEPhenotype.mat','coverage','result','species_withdata')
result_final = num2cell(result(:));
result_final(:,2) = [repmat({'Classic'},length(species_withdata),1);repmat({'DL'},length(species_withdata),1);repmat({'Posterior'},length(species_withdata),1)];
writecell(result_final,'res_RMSEPhenotype.txt','Delimiter',',','QuoteStrings',false)
violin = violinplot(cell2mat(result_final(:,1)),result_final(:,2),'ShowNotches',false,'ShowMean' ,false,'ViolinAlpha',1,'EdgeColor',[0,0,0],'ShowData',false,'BoxColor',[1,1,1]);
violin(1).ViolinColor = [253,224,221]./255;
violin(2).ViolinColor = [250,159,181]./255;
violin(3).ViolinColor = [197,27,138]./255;
violin(2).MedianPlot.SizeData = 1;
violin(1).MedianPlot.SizeData = 1;
violin(3).MedianPlot.SizeData = 1;
xticklabels({'classic','DL','Posterior'})
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['RMSE for phenotype prediction'],'FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
% remove the tick from right axis and top axis
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
           'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
           'XAxisLocation','top', 'XTick', [],... 
           'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])
cd(currentpath)
