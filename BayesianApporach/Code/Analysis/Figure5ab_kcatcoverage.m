% kcatCoverage Figure 5a b
% this function alculate the kcat coverage foe the enzymes in the model
% enzyme coverage
load('StrianData.mat')
species = StrianData.strains;
current_path = pwd;
cd ../../Results/model_build_files/
for i = 1:length(species)
    i
    cd('model_classic')
    z = load([species{i},'_classic.mat']);
    model = z.model;
    enzymedata = z.enzymedata;
    coverage(i,3) = length(find(~cellfun(@isempty,model.grRules)));
    strain = z.strain;
    growthdata = z.growthdata;
    coverage(i,1) = length(enzymedata.kcat);
    
    %if length(enzymedata.proteins) == length(model.genes)
    res_enzyme(i,3) = length(model.genes);
    %else
   %     warning(['check enzymedata.proteins field num with model.genes field for species: ', species{i}]);
    %end
    a = join(enzymedata.enzyme, ' and ');
    a = split(a,' and ');
    allgene = unique(a); % enzyme that are included in the ecmodel
    res_enzyme(i,1) = length(allgene);
    
    cd ../model_dl
    z = load([species{i},'_dl.mat']);
    model = z.model;
    enzymedata = z.enzymedata;
    coverage(i,2) = length(enzymedata.kcat);
     a = join(enzymedata.enzyme, ' and ');
    a = split(a,' and ');
    allgene = unique(a);
    res_enzyme(i,2) = length(allgene);
    cd ../
end
coverage = coverage(:,1:2)./coverage(:,3);
result_final = num2cell(coverage(:));
result_final(:,2) = [repmat({'Classic'},length(species),1);repmat({'DL&Posterior'},length(species),1)];
writecell(result_final,'res_rxnCoverage.txt','Delimiter',',','QuoteStrings',false)
violin = violinplot(cell2mat(result_final(:,1)),result_final(:,2),'ShowNotches',false,'ShowMean' ,false,'ViolinAlpha',1,'EdgeColor',[0,0,0],'ShowData',false,'BoxColor',[1,1,1]);
violin(1).ViolinColor = [253,224,221]./255;
violin(2).ViolinColor = [197,27,138]./255;
violin(2).MedianPlot.SizeData = 1;
violin(1).MedianPlot.SizeData = 1;
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['{\itk}_c_a_t coverage for enzymatic rxns'],'FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
    'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
    'XAxisLocation','top', 'XTick', [],...
    'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])

save('res_rxnCoverage.mat','coverage','species')
res_enzyme = res_enzyme(:,1:2)./res_enzyme(:,3);
result_final = num2cell(res_enzyme(:));
result_final(:,2) = [repmat({'Classic'},length(species),1);repmat({'DL&Posterior'},length(species),1)];
writecell(result_final,'res_enzymeCoverage.txt','Delimiter',',','QuoteStrings',false)
save('res_enzymeCoverage.mat','res_enzyme','species')
violin = violinplot(cell2mat(result_final(:,1)),result_final(:,2),'ShowNotches',false,'ShowMean' ,false,'ViolinAlpha',1,'EdgeColor',[0,0,0],'ShowData',false,'BoxColor',[1,1,1]);
violin(1).ViolinColor = [253,224,221]./255;
violin(2).ViolinColor = [197,27,138]./255;
violin(2).MedianPlot.SizeData = 1;
violin(1).MedianPlot.SizeData = 1;
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['{\itk}_c_a_t coverage for enzymes'],'FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box off
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 10,...
    'Color','None','XColor','k','YColor','k', 'LineWidth', 1,...
    'XAxisLocation','top', 'XTick', [],...
    'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])
cd(current_path)