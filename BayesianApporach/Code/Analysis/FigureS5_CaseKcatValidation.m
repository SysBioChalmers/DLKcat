% Figure for validation of kcat

%% kcat for xylose ultilization

cd ../../Results/model_dl/

fid2 = fopen('SubstrateUsageInfo.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubBiologName = data(:,1);
SubModelName = data(:,2);
SubCondition = data(:,3);
Subtype = data(:,4);
data(:,1:4) = [];
fclose(fid2);

xyloseSpecies = strainlist(strcmp(data(strcmp(SubModelName,'D-xylose'),:),'1'));
NoxyloseSpecies = strainlist(strcmp(data(strcmp(SubModelName,'D-xylose'),:),'0'));

clearvars -except xyloseSpecies NoxyloseSpecies

% get rxnTargets
load('panmodel.mat')
model = panmodel;
        % sol the model first by maxmize ATP yield
        model_tmp = changeMedia(model,'D-xylose','MIN');
          model_tmp = changeRxnBounds(model_tmp,'r_4046',0,'l');
        model_tmp = changeRxnBounds(model_tmp,'r_4046',1000,'u');
        model_tmp = changeObjective(model_tmp,'r_4046',1);
        sol = optimizeCbModel(model_tmp,'max','one');
rxnTarget = model.rxns(find(sol.x));

inputpath = pwd;
[siginificantEnzyme,pvalue,kcatresult,pathwayEnzyme] = kcatSignAnalysis(xyloseSpecies,NoxyloseSpecies,inputpath,0.01,rxnTarget);


% box plot
 species_target = {'Meyerozyma_guilliermondii','Scheffersomyces_stipitis',...
'Spathaspora_arborariae','Spathaspora_passalidarum','yHMPu5000041693_Debaryomyces_nepalensis'};
[~,idx] = ismember([species_target,{'Yarrowia_lipolytica','Saccharomyces_cerevisiae'}],kcatresult.species);
idx2 = find(contains(kcatresult.rxns,{'r_1093','p_4816','r_1092'}));

a = kcatresult.kcat(idx2,idx);
h = bar(log10(a/3600),'FaceColor',[178, 24, 43]/255,'LineWidth',0.5,'FaceAlpha',0.3,'BarWidth',0.5,'EdgeColor',[178, 24, 43]/255)
[r,p] = ttest2(log10(a(1:5))/3600,log10(a(7)/3600));
text(1,1.2,['p = ' num2str(round(p,3))],'FontSize',6,'FontName','Helvetica','Color','k')
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Predicted {\itk}_c_a_t in log10 scale','FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
xticklabels({'mgu','pic','sar','spa','dne','yli','sce'})
set(gca,'FontSize',6,'FontName','Helvetica');
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])




%% fructose HGT case
%fructoseHGT = {'Starmerella_bombicola_JCM9596','Saitoella_complicata','Candida_apicola','Botrytis_cinerea','Sclerotinia_sclerotiorum','Stagonospora_nodorum','Wickerhamiella_domercqiae','Xylona_heveae','yHMPu5000026142_Citeromyces_matritensis','yHMPu5000034660_Diddensiella_caesifluorescens','yHMPu5000034865_Zygosaccharomyces_kombuchaensis','yHMPu5000034866_Zygosaccharomyces_bisporus','yHMPu5000034950_Citeromyces_hawaiiensis','Candida_versatilis','Zygosaccharomyces_bailii','Zygosaccharomyces_rouxii','Fusarium_graminearum','Aspergillus_nidulans','Aspergillus_nidulans'};
fructoseHGT = {'Zygosaccharomyces_bailii','Zygosaccharomyces_rouxii','Wickerhamiella_domercqiae','Starmerella_bombicola_JCM9596','yHMPu5000034865_Zygosaccharomyces_kombuchaensis'}

fructoseTrans = {'p_5228','r_1134'};
glcTrans = {'r_1166'};
for i = 1:length(fructoseHGT)

    load([fructoseHGT{i},'_dl.mat'])
    idx1 = find(contains(enzymedata.rxn_list,fructoseTrans));
    idx2 = find(contains(enzymedata.rxn_list,glcTrans));
    kcatFru(i) = max(enzymedata.kcat(idx1));
    kcatGlc(i) = max(enzymedata.kcat(idx2));
end
h = bar([log10(kcatFru/3600);log10(kcatGlc/3600)]','FaceColor','flat','LineWidth',0.5,'FaceAlpha',0.3,'BarWidth',0.5)
[r,p] = ttest2(log10(kcatFru/3600),log10(kcatGlc/3600))
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Predicted {\itk}_c_a_t [log10]','FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
xticklabels({'zba','zro','wdo','sbo','zko','yli','sce'})
set(gca,'FontSize',6,'FontName','Helvetica');
box off;
ax1 = gca;
ax2 = axes('Position', get(ax1, 'Position'), 'FontSize', 6,...
'Color','None','XColor','k','YColor','k', 'LineWidth', 0.5,...
'XAxisLocation','top', 'XTick', [],...
'YAxisLocation','right', 'YTick', []);
linkaxes([ax1, ax2])