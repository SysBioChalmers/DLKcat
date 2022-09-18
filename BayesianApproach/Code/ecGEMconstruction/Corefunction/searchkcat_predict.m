function [kcat_subunit,kcat_subunit_conf,kcat_all] = searchkcat_predict(rxn,subunit,sub,sub_prod,data,conserveMets)

% -1 in the conf means the kcat is from prediction
% -3 in the coef means the kcat is not found
% -2 in the coef means the kcat is from median value
% index the rxn in the kcatprediction file


% match rxn
[~,idx_rxn] = ismember(rxn,data(:,1));
if idx_rxn ~= 0
    % match substrate
    subunits_file = split(data(idx_rxn,6),';');
    sub_file = split(data(idx_rxn,3),';'); % sub_file is from the prediction result
    sub = split(sub,'; '); % sub is from the rxn
    tmp = setdiff(sub_file,sub);
    if ~isempty(tmp)
        warning(['predictedkcat file does not contain the sub in the model, please check the sub:',cell2mat(join(tmp,';')),' in rxn',cell2mat(rxn)]);
    end
    
    kcat_tmp = split(data(idx_rxn,8),';');
    kcat_tmp = strrep(kcat_tmp,'#','0');
    emptykcatidx = cellfun(@isempty,kcat_tmp);
    kcat_subunit = zeros([length(sub_file) length(subunits_file)]);
    kcat_all = zeros([length(sub_file) length(subunits_file)]);
    kcat_subunit(~emptykcatidx,:) = cellfun(@str2num,split(kcat_tmp(~emptykcatidx),','));
    kcat_all = cellfun(@str2num,split(kcat_tmp(~emptykcatidx),','));
    
    % filter unwanted kcats for currency mets usch as h2O ATP nadph nadh
    sub_new = sub_file;
    if length(intersect(sub_prod,conserveMets.pair)) >= 2
        sub_new = setdiff(sub_new,conserveMets.pair);% only remove conservedMets from sub if there are paired for example NADH and NAD, then it can save kcats prediction for NAD synthesis pathway
    end
    sub_new = setdiff(sub_new,conserveMets.currency);
    
    if isempty(sub_new)
        sub_new = sub_file;
        warning(['keep the sub '  cell2mat(join(sub_file,';')),' for the rxn ',cell2mat(rxn)])
    end
    [~,idx_sub] = ismember(sub_new,sub_file); % get idx for non-currency mets
    kcat_subunit = kcat_subunit(idx_sub(idx_sub~=0),:);
    kcat_subunit = max(kcat_subunit,[],1); % find max kcats for all subs for each subunit
    
    % get idx for subunit
    
    [~,idx_subunit] = ismember(subunit,subunits_file);
    if isempty(kcat_subunit)
        warning(['get the kcat for the rxn',cell2mat(rxn)])
        kcat_subunit = zeros(length(subunit));
        kcat_subunit_conf = repmat(-3,length(subunit),1);
        kcat_all = 0;
    else
        kcat_subunit = kcat_subunit(idx_subunit).*3600; % change from /s to /h
        kcat_subunit_conf = repmat(-1,length(subunit),1);
        kcat_all = kcat_all.*3600; % change from /s to /h
    end
else
    kcat_subunit = zeros(length(subunit));
    kcat_subunit_conf = repmat(-3,length(subunit),1);
    kcat_all = 0;
end

if isempty(kcat_all)
    kcat_all = 0;
end

if length(kcat_all(1,:)) == 1
    kcat_all = reshape(kcat_all,1,length(kcat_all(:,1)));
end

end

