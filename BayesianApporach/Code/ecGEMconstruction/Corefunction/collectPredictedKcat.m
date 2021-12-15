function enzymedata = collectPredictedKcat(model,MWdata,ecdata,Protein_stoichiometry,inputpath,withmedian)
% This function is to load predicted kcat and to generate a enzymedata
% format

current_path = pwd;
species = strrep(model.id,' specific model genereted from panYeast','');

load('conserveMets.mat')
% load predicted kcat
cd(inputpath)
fileName = [species,'_PredictionResults.txt'];
fid  = fopen(fileName);
data = textscan(fid,[repmat('%s ',[1,7]) '%s'],'Delimiter','\t');
data = [data{1:end}];
fclose(fid);
cd(current_path)

% get enzyme list and corresponding rxn list
enzyme_list = model.grRules(~cellfun(@isempty,model.grRules)); % all enyzmes in the model
rxn_list = model.rxns(~cellfun(@isempty,model.grRules));

tf_tmp = count(enzyme_list,' and ');
max_subunit = max(tf_tmp) + 1;% max subunit for constrain the matrix nrxn*length subunit

enzymedata = struct();
enzymedata.proteins = MWdata.genes;
enzymedata.proteinMW = MWdata.MW;
enzymedata.enzyme = enzyme_list;
enzymedata.rxn_list = rxn_list;
enzymedata.substrate = cell(length(enzyme_list),1);
enzymedata.sub_prod = cell(length(enzyme_list),1);
enzymedata.subunit = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_stoichiometry = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_ec = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_conf = nan(length(enzyme_list),max_subunit);
enzymedata.kcat = zeros(length(enzyme_list),1);
enzymedata.kcat_var = zeros(length(enzyme_list),1);
enzymedata.kcat_conf = zeros(length(enzyme_list),1);
enzymedata.MW = zeros(length(enzyme_list),1);
enzymedata.kcat_max = zeros(length(enzyme_list),1);
enzymedata.subunit_MW = zeros(length(enzyme_list),max_subunit);
for i = 1:length(enzyme_list)
    if mod(i,1000) == 0
        disp(['Assigning kcats:' num2str(i) '/' num2str(length(enzyme_list))]);
    end
    enzyme_id = enzyme_list{i};
    
    % add substrates
    metrxn_id = rxn_list(i);
    idx_tmp = ismember(model.rxns,metrxn_id);
    s_tmp = model.S(:,idx_tmp);
    mets_tmp = model.metNames(s_tmp < 0);
    mets_tmp = cellfun(@(x) x(1:strfind(x,' [')-1),mets_tmp,'UniformOutput',false);
    mets_tmp = strjoin(mets_tmp,'; ');
    enzymedata.substrate(i,1) = {mets_tmp};
    mets_tmp = model.metNames(s_tmp ~= 0);
    mets_tmp = cellfun(@(x) x(1:strfind(x,' [')-1),mets_tmp,'UniformOutput',false);
    mets_tmp = strjoin(mets_tmp,'; ');
    enzymedata.sub_prod(i,1) = {mets_tmp};
    
    % add subunits
    idx_tmp = ismember(model.rxns,metrxn_id);
    subunits_tmp = model.genes(logical(model.rxnGeneMat(idx_tmp,:)));
    subunits_pan = model.proteins(logical(model.rxnGeneMat(idx_tmp,:)));
    enzymedata.subunit(i,1:length(subunits_tmp)) = subunits_tmp;
    
    % add stoichiometry of subunits
    enzymedata.subunit_stoichiometry(i,1:length(subunits_tmp)) = ones(length(subunits_tmp),1);
    [~,idx] = ismember(subunits_pan,Protein_stoichiometry.protein);
    enzymedata.subunit_stoichiometry(i,idx~=0) = Protein_stoichiometry.stoichiometry(idx(idx~=0));
    
    % add ec numbers of subunits
    for j = 1:length(subunits_pan) % using panID to index ec number
        subunit_tmp = subunits_pan(j);
        if ismember(subunit_tmp,ecdata.id)
            ec_tmp = unique(ecdata.ec(ismember(ecdata.id,subunit_tmp)));
            ec_tmp = strjoin(ec_tmp,'; ');
            kcatrange_tmp(1,:) = unique(ecdata.kcatrange(ismember(ecdata.id,subunit_tmp),:));
        else
            ec_tmp = 'NA';
            kcatrange_tmp = 0;
        end
        enzymedata.subunit_ec(i,j) = {ec_tmp};
        enzymedata.subunit_ec_kcat_range{i,j} = kcatrange_tmp;
        clearvars  kcatrange_tmp
    end
    kcat_range_tmp = [enzymedata.subunit_ec_kcat_range{i,:}];
    enzymedata.enzyme_ec_kcat_range(i,1) = min(kcat_range_tmp);
    enzymedata.enzyme_ec_kcat_range(i,2) = max(kcat_range_tmp);
    
    % add MW of subunits
    [~,idx] = ismember(subunits_tmp,MWdata.genes);
    enzymedata.subunit_MW(i,1:length(idx)) = MWdata.MW(idx);
    enzymedata.MW(i,1) = sum(enzymedata.subunit_MW(i,:).*enzymedata.subunit_stoichiometry(i,:));
    
    % assign kcats and confidence scores for subunits
    substrate_tmp = enzymedata.substrate(i,1);
    sub_prod_tmp = enzymedata.substrate(i,1);
    [finalkcat_tmp,kcat_subunit_coef_tmp,kcat_all_tmp] = searchkcat_predict(metrxn_id,subunits_tmp,substrate_tmp,sub_prod_tmp,data,conserveMets);
    enzymedata.subunit_kcat(i,1:length(subunits_tmp)) = finalkcat_tmp;
    enzymedata.subunit_kcat_conf(i,1:length(subunits_tmp)) = kcat_subunit_coef_tmp; % using the predictiong value
    
    % assign kcats and confidence scores for enzymes
    conf_tmp = max(enzymedata.subunit_kcat_conf(i,:));
    enzymedata.kcat_conf(i,1) = conf_tmp;
    finalkcat_list = enzymedata.subunit_kcat(i,:).*enzymedata.subunit_stoichiometry(i,:); %consider protein stoichiometry
    finalkcat_all = kcat_all_tmp.*enzymedata.subunit_stoichiometry(i,1:length(subunits_tmp));
    kcat_tmp = finalkcat_list(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    if ~isempty(kcat_tmp(kcat_tmp~=0))
        enzymedata.kcat(i,1) = max(kcat_tmp(kcat_tmp~=0)); %choose median/mean/max among subunits
    else
        enzymedata.kcat(i,1) = 0; %choose median among subunits
    end
    enzymedata.kcat_var(i,1) = std(log10(finalkcat_all(:)/3600));
end

% no kcat assigned enzyme assumed to be the median of the collected dataset


% 1st enzymedata median value
if withmedian
    medianvalue = median(enzymedata.kcat(enzymedata.kcat ~= 0));
    enzymedata.kcat_var(1:end) = 3;
    enzymedata.kcat_var(enzymedata.kcat == 0) = 4;
    enzymedata.kcat_conf(enzymedata.kcat == 0) = -2;
    enzymedata.kcat(enzymedata.kcat == 0) = medianvalue;
else
    % 2 st enzymedata
    idx = (enzymedata.kcat == 0 | isnan(enzymedata.kcat));
    enzymedata.enzyme(idx) = [];
    enzymedata.rxn_list(idx) = [];
    enzymedata.substrate(idx) = [];
    enzymedata.sub_prod(idx) = [];
    enzymedata.subunit(idx,:) = [];
    enzymedata.subunit_stoichiometry(idx,:) = [];
    enzymedata.subunit_ec(idx,:) = [];
    enzymedata.subunit_kcat(idx,:) = [];
    enzymedata.subunit_kcat_conf(idx,:) = [];
    enzymedata.kcat_var(idx,:) = [];
    enzymedata.kcat_conf(idx) = [];
    enzymedata.MW(idx) = [];
    enzymedata.kcat_max(idx,:) = [];
    enzymedata.subunit_MW(idx,:) = [];
    enzymedata.enzyme_ec_kcat_range(idx,:) = [];
    enzymedata.subunit_ec_kcat_range(idx,:) = [];
    enzymedata.kcat(idx,:) = [];
end

enzymedata.kcat_var(enzymedata.kcat_var <= 1) = 1;
enzymedata.kcat_var(isnan(enzymedata.kcat_var)) = 1; % 1 is from the deep_learning model
% we are assuming that variance from the kcat prediction are the same, or
% should we consider the total variance from the database?

% and also when we have a predicted high kcat but we need to change the
% kcat into the range
tmp = log10(enzymedata.kcat./3600) - enzymedata.enzyme_ec_kcat_range(:,1);
idx = (tmp < 0 & enzymedata.enzyme_ec_kcat_range(:,1) ~= 0);
enzymedata.enzyme_ec_kcat_range(idx,1) = log10(enzymedata.kcat(idx)./3600);
enzymedata.kcat_conf(idx) = -4;
tmp = log10(enzymedata.kcat./3600) - enzymedata.enzyme_ec_kcat_range(:,2);
idx = (tmp > 0 & enzymedata.enzyme_ec_kcat_range(:,2) ~= 0);
enzymedata.enzyme_ec_kcat_range(idx,2) = log10(enzymedata.kcat(idx)./3600);
enzymedata.kcat_conf(idx) = -4;
cd(current_path)
end