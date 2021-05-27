%% this function is for figure S6b

currentpath = pwd;
load('conserveMets.mat')
cd ../../Results/model_dl/
fid2 = fopen('343_phenotype_clade.tsv');
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
%% sort all OG groups for all species

load('ecdata.mat')
for k = 1:length(clades)
    species = Strain_information(ismember(Strain_information(:,2),clades(k)),1);
% define specialist & genralist
kcat_general = [];
kcat_special = [];

for i = 1:length(species)
    disp([num2str(i),'/343 species'])
    z = load([species{i},'_dl.mat']);
    model = z.model;
   enzymedata = z.enzymedata;
   res = tabulate(enzymedata.enzyme);
   spe = res(cell2mat(res(:,2)) == 1);
   gen = res(cell2mat(res(:,2)) > 1);
   idx_gen = ismember(enzymedata.enzyme,gen);
    kcat_general = [kcat_general;enzymedata.kcat(idx_gen)];
    idx_spe = ismember(enzymedata.enzyme,spe);
    kcat_special = [kcat_special;enzymedata.kcat(idx_spe)];
    % index the species
    kcat_general = kcat_general(kcat_general~=0);
    kcat_special = kcat_special(kcat_special~=0);
    clearvars data datanew datasorted model reducedModel
end
kcat_special_clade{k} = kcat_special;
kcat_general_clade{k} = kcat_general;
clear  kcat_general_splited_class_name kcat_general_splited_class kcat_special_splited_class kcat_special_splited_class_name
end


% to keep consistency of boxplot, we use python to plot the figure
res = [];
res_label1 = [];
res_lable2 = [];
for i = 1:13
    res = [res;kcat_special_clade{i}];
    res_label1 = [res_label1;repmat({'spe'},length(kcat_special_clade{i}),1)];
    res = [res;kcat_general_clade{i}];
    res_label1 = [res_label1;repmat({'gen'},length(kcat_general_clade{i}),1)];
    res_lable2 = [res_lable2;repmat(i,length(kcat_special_clade{i}) + length(kcat_general_clade{i}),1)];
end
t = cell2table([num2cell(log10(res/3600)),res_label1,num2cell(res_lable2)]);
t.Properties.VariableNames = {'kcat','type','clade'};
writetable(t,'kcat_gen_spe.txt','Delimiter','\t','QuoteStrings',false,'WriteRowNames',true)
t = cell2table([num2cell(log10(res/3600)),res_label1,num2cell(res_lable2)]);
writetable(t,'kcat_gen_spe.txt','Delimiter','\t','QuoteStrings',false,'WriteRowNames',true)