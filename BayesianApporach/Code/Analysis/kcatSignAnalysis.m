function [siginificantEnzyme,pvalue,kcatresult,pathwayEnzyme] = kcatSignAnalysis(group1,group2,inputpath,Pcutoff,rxnsTarget)

% This function is to find out the siginificant Enzyme for certain groups
% with certain types.

if nargin < 4
    Pcutoff = 0.05;
    rxnsTarget = {};
end
currentpath = pwd;

cd(inputpath)
species = [group1(:);group2(:)];

m = [];

for i = 1:length(species)
    disp([num2str(i),'/',num2str(length(species))])
    z = load([species{i},'_dl.mat']);
    enzymedata = z.enzymedata;
    pre = extractBefore(enzymedata.rxn_list,'_');
    tmp = extractAfter(enzymedata.rxn_list,'_');
    tmp = regexprep(tmp,'_(\d)','');
    enzymedata.rxn_list = strcat(pre,'_',tmp);
    if i == 1
        corxns = enzymedata.rxn_list;
    else
        corxns = intersect(corxns,enzymedata.rxn_list);
    end
    m(i) = length(corxns);
end
kcatresult.rxns = corxns;
for i = 1:length(species)
    disp([num2str(i),'/',num2str(length(species))])
    z = load([species{i},'_dl.mat']);
    enzymedata = z.enzymedata;
    tmp = strrep(enzymedata.rxn_list,'r_','');
    tmp = regexprep(tmp,'_(\d)','');
    enzymedata.rxn_list = strcat('r_',tmp);
    for j = 1:length(corxns)
        idx_rxn = ismember(enzymedata.rxn_list,corxns(j));
        kcatresult.kcat(j,i) = median(enzymedata.kcat(idx_rxn));
    end
end
kcatresult.species = species;

%pval = cell(length(corxns),1);
for i = 1:length(corxns)
[h(i),pval(i),~,~] = ttest2(log10(kcatresult.kcat(i,1:length(group1))),log10(kcatresult.kcat(i,length(group1)+1:end)),'Vartype','unequal');
end
pval = round(pval,2);
siginificantEnzyme = corxns(pval <= Pcutoff);
pvalue = pval(pval < Pcutoff);
pathwayEnzyme = intersect(rxnsTarget,siginificantEnzyme);

cd(currentpath)
end