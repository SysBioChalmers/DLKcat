%% collectkcats
function enzymedata = collectkcats(model,MWdata,ecdata,Protein_stoichiometry,withmedian)
% contains model.proteins field which stores the same id with the ec file
% which we have
org_name = strrep(model.id,' specific model genereted from panYeast','');
org_name = strrep(org_name,'_',' ');
org_name = lower(org_name);
%% Load kcats data from BRENDA database

allkcats = loadkcats;
% generate a ec file for all pangenes %%%%%%important
%% Collect EC number
if nargin < 3
    [ecdata,Protein_stoichiometry] = ECprepPanGEM;
end
%% Assign kcat for each reaction with GPR

% collect enzyme data
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
enzymedata.subunit = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_stoichiometry = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_ec = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_conf = nan(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_median = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_var = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_max = zeros(length(enzyme_list),max_subunit);
enzymedata.kcat = zeros(length(enzyme_list),1);
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
    for j = 1:length(subunits_tmp)
        ec_tmp = enzymedata.subunit_ec(i,j);
        substrate_tmp = enzymedata.substrate(i,1);
        [finalkcat_tmp, conf_score_tmp,median_ec,std_ec,max_ec] = searchkcat(ec_tmp,substrate_tmp,org_name,allkcats);
        enzymedata.subunit_kcat(i,j) = finalkcat_tmp;
        enzymedata.subunit_kcat_conf(i,j) = conf_score_tmp;
        enzymedata.subunit_kcat_median(i,j) = median_ec;
        enzymedata.subunit_kcat_var(i,j) = std_ec;
        enzymedata.subunit_kcat_max(i,j) = max_ec;
    end
    
    % assign kcats and confidence scores for enzymes
    conf_tmp = max(enzymedata.subunit_kcat_conf(i,:));
    enzymedata.kcat_conf(i,1) = conf_tmp;
    finalkcat_list = enzymedata.subunit_kcat(i,:).*enzymedata.subunit_stoichiometry(i,:); %consider protein stoichiometry
    kcat_tmp = finalkcat_list(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
%   enzymedata.kcat(i,1) = median(kcat_tmp); %choose median among subunits
    enzymedata.kcat(i,1) = min(kcat_tmp); %choose minimum among subunits
    finalkcat_median = enzymedata.subunit_kcat_median(i,:).*enzymedata.subunit_stoichiometry(i,:); %consider protein stoichiometry
    median_tmp = finalkcat_median(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    enzymedata.kcat_median(i,1) = min(median_tmp); %choose minimum among subunits
    finalkcat_var = enzymedata.subunit_kcat_var(i,:); %consider protein stoichiometry
    var_tmp = finalkcat_var(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    enzymedata.kcat_var(i,1) = max(var_tmp); %choose max among subunits
    finalkcat_max = enzymedata.subunit_kcat_max(i,:).*enzymedata.subunit_stoichiometry(i,:); %consider protein stoichiometry
    max_tmp = finalkcat_max(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    enzymedata.kcat_max(i,1) = max(max_tmp); %choose max among subunits
end

% 1st no kcat assigned enzyme assumed to be the median of the collected dataset 
if withmedian
    medianvalue = median(enzymedata.kcat(~isnan(enzymedata.kcat)));
    medianvar = median(enzymedata.kcat_var(~isnan(enzymedata.kcat)));
    enzymedata.kcat_var(isnan(enzymedata.kcat)) = medianvar;
    enzymedata.kcat_var(enzymedata.kcat_var == 0) = medianvar;
    enzymedata.kcat_var(enzymedata.kcat_var < 3) = 3;
    %enzymedata.kcat_median(enzymedata.kcat_median == 0) = medianvar;
    enzymedata.kcat_median(isnan(enzymedata.kcat)) = medianvalue;
    enzymedata.kcat(isnan(enzymedata.kcat)) = medianvalue;
else
    % 2 nd no kcat assigned enzyme not taken into account
    idx = (enzymedata.kcat == 0 | isnan(enzymedata.kcat));
    enzymedata.kcat(idx) = 0;
    enzymedata.enzyme(idx) = [];
    enzymedata.rxn_list(idx) = [];
    enzymedata.substrate(idx) = [];
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
    enzymedata.kcat(idx,:) = [];
    enzymedata.enzyme_ec_kcat_range(idx,:) = [];
    enzymedata.subunit_ec_kcat_range(idx,:) = [];
    enzymedata.subunit_kcat_median(idx,:) = [];
    enzymedata.subunit_kcat_var(idx,:) = [];
    enzymedata.subunit_kcat_max(idx,:) = [];
    enzymedata.kcat_median(idx) = [];
end
enzymedata.kcat_var(enzymedata.kcat_var < 1) = 1;

% There are a few cases that multiple isozymes are assigned different kcats
% due to different EC numbers, e.g., EC4.2.1.- and EC4.2.1.36. This should
% be avoided by re-assigning the kcat with higher confidence score for the
% isozyme with lower confidence score?



end



    