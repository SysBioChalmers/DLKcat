% This analysis is to test whether the model proteome is similar to the
% measured proteome
current_path = pwd;
% load proteome data
[~,~,raw_proteome] = xlsread('Proteome_ref.xlsx','alldata');
raw_proteome = raw_proteome(2:end,:); % the first is the head line
cond = unique(raw_proteome(1,:),'stable');
lable  = unique(raw_proteome(2,:),'stable');
speciesid = {'Sce','Kma','Kla','Yli'};
species = {'Saccharomyces_cerevisiae','Kluyveromyces_marxianus','Kluyveromyces_lactis','Yarrowia_lipolytica'};

states = {'auto','DL','Bayesian_DL_mean'}; % 'Bayesian_DL' and 'global' arealsooption that can ploted global means a global kcat would be used for all enzymes Bayesian_DL means that all posterior models will be used

cd ../../Results
for i = 1:length(cond)
    for k = 1:length(states)
        k
        state = states{k};
        disp([num2str(i),'/',num2str(length(cond))])
        %subplot(4,6,i)
        code = cond{i};
        code = split(code,'_');
        [~,idx] = ismember(code(1),speciesid);
        code(end+1) = species(idx);
        code(end+1) = {['model_bayesian_max/',code{7}]};
        % get flux
        [model,sol_result_mean] = getflux(code{7},code{2},code{3},code{4},code{5},code{8},state); % species chemostat/batch media o2 carbon inputpath
        sol_result{i,k} = sol_result_mean;
        models{i} = model;
        % get proteome
        idx = ismember(raw_proteome(1,:),cond(i));
        raw_proteome_tmp = raw_proteome(4:end,idx);
        if strcmp(code{1},'Sce') || strcmp(code{6},'protein')
            proteome.gene = raw_proteome_tmp(:,1); % geneID
            if strcmp(code{6},'protein')
                proteome.gene = convert2panID(proteome.gene); % change to S.cereviase gene
            end
        else
            proteome.gene = raw_proteome_tmp(:,2); % geneID in the model species@seq1
        end
        proteome.abun = cell2mat(raw_proteome_tmp(:,end)); % protein abundance in mmol
        proteome.abun(~cell2mat(cellfun(@ischar,proteome.gene,'UniformOutput',false))) = []; % delete empty line
        proteome.gene(~cell2mat(cellfun(@ischar,proteome.gene,'UniformOutput',false))) = []; % delete empty line
        proteome.gene(proteome.abun == 0|isnan(proteome.abun)) = []; % delete 0 abundance protein
        proteome.abun(proteome.abun == 0|isnan(proteome.abun)) = []; % delete 0 abundance protein
        if strcmp(code{1},'Sce') || strcmp(code{6},'protein')
            [~,idx] = ismember(model.proteins,proteome.gene);
        else
            [~,idx] = ismember(model.genes,proteome.gene);
        end
        proteome_abun = zeros(length(model.genes),1);
        proteome_abun(idx~=0) = proteome.abun(idx(idx~=0));
        proteome_abun = proteome_abun.*model.MWs(1:length(model.genes));
        enzUsage_abun = calculateproteome(model,sol_result{i,k},model.genes);
        
        
        abun(:,1) = proteome_abun;
        gene(:,1) = model.proteins;
        abun(:,k+1) = enzUsage_abun;
        clear proteome_abun enzUsage_abun
    end
    
    % get pan abundance by sum all paralog abundance
    protein = unique(model.proteins);
    allresult.gene(1:length(protein),1,i) = protein;
    for j = 1:length(protein)
        idx = ismember(model.proteins,protein(j));
        allresult.abun(j,:,i) = sum(abun(idx,:),1);
    end
    clear abun gene   
    
    
    % find congene and calculate the rmse
    tmp_result =  allresult.abun(:,:,i);
    tmp_result(tmp_result < 1E-6) = 0; % lower than the tolerance
    tmp_result_gene =  allresult.gene(:,:,i);
    tmp_result_gene(any(tmp_result==0,2),:) = [];
    tmp_result(any(tmp_result==0,2),:) = [];
    
    for m = 1:length(states)
        % plot for corr
        figure
        [RHO4,PVAL4] = corr(log10(tmp_result(:,1)),log10(tmp_result(:,1+m)),'Type','Pearson');
        scatter(log10(tmp_result(:,1)),log10(tmp_result(:,1+m)),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
        text(median(log10(tmp_result(:,1))),median(log10(tmp_result(:,1+m))),['r = ' num2str(round(RHO4,3))],'FontSize',6,'FontName','Helvetica');
        text(median(log10(tmp_result(:,1))),median(log10(tmp_result(:,1+m))-1),['p = ' num2str(round(PVAL4,10))],'FontSize',6,'FontName','Helvetica');
        title([strrep(cond{i},'_','\_'),' (n = ',num2str(length(tmp_result_gene)),')'],'FontSize',8,'FontName','Helvetica');
        result.pearsonr(i,m) = RHO4;
        result.pvalue(i,m) = PVAL4;
        result.rmse(i,m) = sqrt(immse(log10(tmp_result(:,1)),log10(tmp_result(:,1+m))));
        result.genenum(i,m) = length(tmp_result_gene);
    end
end

% figure 
colorveryhigh = [254,235,226]/255;
colorhigh = [251,180,185]/255;
%colormedium = [247,104,161]/255;
colorlow = [174,1,126]/255;
figure();
b = bar(1:22,result.rmse');
b(1).LineWidth = 0.5;
b(1).FaceColor = colorveryhigh;
b(2).LineWidth = 0.5;
b(2).FaceColor = colorhigh;
b(3).LineWidth = 0.5;
b(3).FaceColor = colorlow;
% b(4).LineWidth = 0.5;
% b(4).FaceColor = colormedium;
set(gca,'XTick',1:1:22);
set(gca,'XTickLabel',strrep(lable,'_','\_'));
xtickangle(90)
ylim([0 2.3]);
leg = legend({'Auto' 'DL' 'Bayesian_DL_mean'},'Location','northwest','Orientation','horizontal','FontSize',6,'FontName','Helvetica');
leg.ItemTokenSize = [10,2];
set(leg,'box','off');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('RMSE for proteome prediction','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[200 100 300 150]);
set(gca,'position',[0.1 0.2 0.75 0.5]);
% remove the tick from right axis and top axis
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
           'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
           'XAxisLocation','top', 'XTick', [],... 
           'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])