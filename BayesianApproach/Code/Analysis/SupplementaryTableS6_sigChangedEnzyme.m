species = {'Saccharomyces_cerevisiae','Yarrowia_lipolytica'};
currentpath = pwd;
cd ../../Results/model_build_files/
load('Yeast8.mat')
org_model.subSystems = cellfun(@(x) x{1},org_model.subSystems,'UniformOutput',false);
for j = 1:length(species)
    cd ('model_Bayesian/')
    cd(species{j})
    nfound = length(dir('kcat_genra*.txt'));
    tmp = readmatrix(['kcat_genra',num2str(nfound),'.txt'],'FileType','text','Delimiter',',');
    kcat_posterior = tmp(1:end-2,:);
    ss = num2cell(kcat_posterior',1);
    [a,b] = arrayfun(@updateprior,ss);
    kcat_posterior_mean = a';
    cd ../../model_dl
    z = load([species{j},'_dl.mat']);
    enzymedata = z.enzymedata;
    model = z.model;
    kcat_prior = arrayfun(@getrSample, enzymedata.kcat,enzymedata.kcat_var,enzymedata.enzyme_ec_kcat_range(:,1),enzymedata.enzyme_ec_kcat_range(:,2),repmat(100,length(enzymedata.kcat),1),'UniformOutput',false);
    kcat_prior = cell2mat(kcat_prior);
    [~,p]= ttest2(log10(kcat_posterior),log10(kcat_prior),'Dim',2,'Vartype','unequal');
    p_adj = pval_adjust(p,'sidak');
    sig_enzyme(:,1) = enzymedata.rxn_list(p_adj < 0.01);
    [~,idx] = ismember(sig_enzyme(:,1),enzymedata.rxn_list);
    sig_enzyme(:,2) = num2cell(log10(enzymedata.kcat(idx)/3600));
    sig_enzyme(:,3) = num2cell(log10(kcat_posterior_mean(idx)/3600));
    sig_enzyme(:,4) = num2cell(~cellfun(@isempty,enzymedata.subunit(idx,2)));
    rxns = extractBefore(sig_enzyme(:,1),7);
    [~,idx] = ismember(rxns,org_model.rxns);
    sig_enzyme(idx~=0,5) = org_model.subSystems(idx(idx~=0));
    sig_enzyme(cellfun(@isempty,sig_enzyme(:,5)),5) = {''};
    sig_subsystem = tabulate(sig_enzyme(:,5));
    rxns_all = extractBefore(enzymedata.rxn_list(:,1),7);
    [~,idx] = ismember(rxns_all,org_model.rxns);
    enzyme_all = cell(length(rxns_all),1);
    enzyme_all(idx~=0,1) = org_model.subSystems(idx(idx~=0));
    enzyme_all(cellfun(@isempty,enzyme_all(:,1)),1) = {''};
    tot_subsystem = tabulate(enzyme_all);
    cd ../
    [~,idx] = ismember(sig_subsystem(:,1),tot_subsystem(:,1));
    sig_subsystem(:,3) = tot_subsystem(idx,2);
    % calculate enrichment fisher exact test
    x = cell2mat(sig_subsystem(:,2:3)); % sig_rxn_num tot_rxn_num
    for i = 1:length(sig_subsystem(:,1))
    y(1,1) = x(i,1);
    y(2,1) = x(i,2) -x(i,1);
    y(1,2) = sum(x(:,1)) - x(i,1);
    y(2,2) = sum(x(:,2))- sum(x(:,1)) - x(i,2) + x(i,1);
    [h,p,stats] = fishertest(y);
    p_res(i,1) = p;
    end
    p_adj = pval_adjust(p_res,'BH');
    sig_subsystem(:,4) = num2cell(p_res);
    sig_subsystem(:,5) = num2cell(p_adj);
    save([species{j},'_sig_enzyme.mat'],'sig_enzyme','p_adj','p_res','sig_subsystem')
    clear sig_enzyme sig_enzyme 
end



 