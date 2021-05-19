function kcatCoverage(species,version1)
% this function alculate the kcat coverage foe the enzymes in the model

current_path = pwd;
if strcmp(version1,'dl')
    inputpath = '../KcatTuning/model_dl_new';
    cd(inputpath)
else
    inputpath = '../KcatTuning/model_auto';
    cd(inputpath)
end
% load model
for i = 1:length(species)
    if strcmp(version1,'dl')
        z = load([species{i},'_dl.mat']);
    else
      z = load([species{i},'_auto.mat']);
    end
    enzymedata = z.enzymedata;
    model = z.model;
    % enzyme coverage
    a = join(enzymedata.enzyme, ' and ');
    a = split(a,' and ');
    allgene = unique(a);
    tmp = setdiff(enzymedata.enzyme, allgene);
    res_enzyme(i,1) = length(tmp);
    if length(enzymedata.proteins) == length(z.model.genes)
    res_enzyme(i,2) = length(enzymedata.proteins);
    else
        warning(['check enzymedata.proteins field num with model.genes field for species: ', species{i}]);
    end
    clear a
    % rxn coverage
    a = cellfun(@isempty,model.rules);
    a = find(~a);
    res_rxn(i,1) = length(enzymedata.enzyme);
    res_rxn(i,2) = length(a);
end
clearvars -except species res_enzyme res_rxn current_path version1

cd(current_path)
save(['Results/res_kcatcoverage',version1,'.mat'])