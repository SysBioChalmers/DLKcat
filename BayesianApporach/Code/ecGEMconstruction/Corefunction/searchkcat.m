%% searchkcat
function [finalkcat, conf_score,median_ec,std_ec,max_ec] = searchkcat(ec,substrate,org_name,allkcats)

% Input:  ec, EC number, e.g., 'EC1.1.1.1' or 'EC1.1.2.4; EC1.1.99.40';
%         substrate, metabolite name, e.g., 'acetyl-CoA; oxaloacetate';
%         org_name, organism name, e.g., 'saccharomyces cerevisiae';
%         allkcats, kcats from BRENDA database.
% Output: finalkcat, kcat value in the unit of "/h", max for multiple EC numbers;
%         conf_score, confidence score of the searched kcat value.
% CONFIDENCE SCORE
% 4:    Both substrates and organism matched, then kcat = max of matched kcats.
% 3:    Organism matched but substrates not matched, then kcat = median of all substrates matched kcats.
% 2:    No organism matched, then kcat = median of all matched substrates in other organisms.
% 1:    No organism and substrates matched, then kcat = median of all reported kcats of the EC number.
% 0:    No kcat reported for the EC number in BRENDA, then kcat = nan.

% Note: if multiple EC numbers provided, then the kcat with highest
%       confidence score will be selected.

whole_ec = allkcats{1};
whole_substrate = allkcats{2};
whole_org = allkcats{3};
whole_kcat = allkcats{4};


eclist = split(ec,'; ');
kcatlist = zeros(length(eclist),1);
conflist = zeros(length(eclist),1);
std_ec = zeros(length(eclist),1);
median_ec = zeros(length(eclist),1);
max_ec = zeros(length(eclist),1);
substrate = lower(substrate);
sublist = split(substrate,'; ');
sublist_plus = cellfun(@(x) strcat(x,'+'),sublist,'UniformOutput',false);%nad/nadp and nad+/nadp+
sublist = [sublist;sublist_plus];

for i = 1:length(eclist)
    ec_tmp = eclist(i);
    if ~ismember(ec_tmp,whole_ec)
        kcat = nan;
        conf = 0;
        std_ec_tmp = 0;
        median_ec_tmp = 0;
        max_ec_tmp = 0;
    else
        idx_tmp = ismember(whole_ec,ec_tmp);
        sub_tmp = whole_substrate(idx_tmp);
        org_tmp = whole_org(idx_tmp);
        kcat_tmp = whole_kcat(idx_tmp);
        std_ec_tmp = std(log10(kcat_tmp/3600),0);
        median_ec_tmp = median(kcat_tmp);
        max_ec_tmp = max(kcat_tmp);
        
        match_idx_sub = ismember(sub_tmp,sublist);
        match_idx_org = ismember(org_tmp,org_name);
        match_idx_combined = match_idx_sub & match_idx_org;
        
        if any(match_idx_combined) %both organism and substrate matched
            kcat = max(kcat_tmp(match_idx_combined)); %choose max
            conf = 4;
        else
            if any(match_idx_org) %only organism matched
                kcat = median(kcat_tmp(match_idx_org)); %choose median
                conf = 3;
            else %if no organism matched
                if any(match_idx_sub) %if substrates matched, then choose median
                    kcat = median(kcat_tmp(match_idx_sub)); %choose median
                    conf = 2;
                else
                    kcat = median(kcat_tmp);
                    conf = 1;
                end
            end
        end
        
        % another round of fuzzy matching
        sublist = strrep(sublist,'_','-');
        match_idx_sub = ismember(sub_tmp,sublist);
        match_idx_combined = match_idx_sub & match_idx_org;
        if any(match_idx_combined) %both organism and substrate matched
            kcat = max(kcat_tmp(match_idx_combined)); %choose max
            conf = 4;
            
        end
        kcatlist(i,1) = kcat;
        conflist(i,1) = conf;
        std_ec(i,1) = std_ec_tmp;
        median_ec(i,1) = median_ec_tmp;
        max_ec(i,1) = max_ec_tmp;
    end
    
    finalkcat_tmp = kcatlist(conflist == max(conflist));
    
    finalkcat = max(finalkcat_tmp);
    conf_score = max(conflist);
    std_ec = max(std_ec);
    median_ec = max(median_ec);
    max_ec = max(max_ec);
end
