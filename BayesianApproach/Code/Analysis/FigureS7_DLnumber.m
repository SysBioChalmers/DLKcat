current_path = pwd;
%this function is to plot the figure DL-ecGEMs number for Figure S6

fid2 = fopen('343_phenotype_clade.tsv');% DLKcat/BayesianApproach/Data
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    Strain_information(:,i) = temp{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
group = [];
clade_av = [];

strains_sortclade = [];
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    group = [group;i*ones(length(idx(idx~=0)),1)];
    strains_sortclade = [strains_sortclade;Strain_information(idx~=0,1)];
end
cd ../../Results/model_build_files/model_dl/
for k = 1:length(strains_sortclade)
    spec = strains_sortclade{k};
    z = load([spec,'_dl.mat']);
    idx = cellfun(@isempty,z.enzymedata.subunit);
    prot = unique(z.enzymedata.subunit(~idx));
    result(k,1) = length(prot);
    result(k,2) = length(z.model.genes); % not all genes are included in the prediction due to the substrate SMILES missing
    result(k,3) = length(z.enzymedata.rxn_list);
    subs = join(z.enzymedata.substrate,'; ');
    subs = unique(split(subs,'; '));
    result(k,4) = length(subs);
    
   k 
end


figure
hold on;
h = boxplot(result(:,1),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
set(gca,'xtick',[])
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Protein number','FontSize',7,'FontName','Helvetica','Color','k');
camroll(-90)
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
box on


figure
h = boxplot(result(:,3),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
set(gca,'xtick',[])
yticks([0:2000:6000])
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Enzymatic reaction number','FontSize',7,'FontName','Helvetica','Color','k');
camroll(-90)
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
box on

figure
h = boxplot(result(:,4),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
set(gca,'xtick',[])
ylim([800 1400])

ylabel('Substrate number','FontSize',7,'FontName','Helvetica','Color','k');
camroll(-90)

set(gca,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
box on

cd (current_path)