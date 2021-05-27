% subsystem kcat
load('strains.mat')
species = strains

fid2 = fopen('module_ec.txt');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
module_ec(:,i) = temp{i};
end
fclose(fid2);

    cd('model_dl')
    kcat_Ec_sorted = '';
    kcat_sorted = [];
    kcat_pathway_sorted = '';
for i = 1:length(species)
    % load data
    i
    load([species{i},'_dl.mat'])
    kcat = enzymedata.subunit_kcat(:);
    kcat_Ec = enzymedata.subunit_ec(:);
    idx = strcmp(kcat_Ec,'NA')| cellfun(@isempty,kcat_Ec);
    kcat = kcat(~idx);
    kcat_Ec = kcat_Ec(~idx);
    [~,idx] = ismember(kcat_Ec,module_ec(:,2));
    kcat_Ec_sorted = [kcat_Ec_sorted;kcat_Ec(idx~=0)];
    kcat_sorted = [kcat_sorted;kcat(idx~=0)];
   kcat_pathway_sorted = [kcat_pathway_sorted;module_ec(idx(idx~=0),3)];
    kcat_Ec = kcat_Ec(idx == 0);
    kcat = kcat(idx == 0);
    for j = 1:length(kcat_Ec)
        ec_split = strrep(split(kcat_Ec{j},';'),' ','');
        [~,idx] = ismember(ec_split,module_ec(:,2));
        path_tmp = unique(module_ec(idx(idx~=0),3));
        path_tmp = setdiff(path_tmp,'');
        if ~isempty(path_tmp)
            kcat_sorted = [kcat_sorted;repmat(kcat(j),length(path_tmp),1)];
            kcat_Ec_sorted = [kcat_Ec_sorted;repmat(kcat_Ec(j),length(path_tmp),1)];
            kcat_pathway_sorted = [kcat_pathway_sorted;path_tmp];
        end
    end
end
label = {'Primary - Carbohydrate & Energy Metabolism','Primary - amino acids, fatty acids and nucleotides','Intermediate','Secondary'}
figure
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;

for i = 1:5
    idx = strcmp(kcat_pathway_sorted,label{i});
    hold on
    h = cdfplot(log10(kcat_sorted(idx)./3600));
    h.Color = color(i,:);
end
box off
ylabel('Cumulative distribution','FontSize',8,'FontName','Helvetica','Color','k');
legend({'Primary-CE','Primary-AFN','Intermediate','Secondary'},'FontSize',6)
set(gcf,'position',[600 0 450 170]);
set(gca,'position',[0.1 0.4 0.27 0.56]);
saveas(z2,['growthfigure_',species{i},'.pdf']);