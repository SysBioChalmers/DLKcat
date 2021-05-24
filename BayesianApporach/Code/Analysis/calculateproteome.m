function enzUsage_abun = calculateproteome(model,sol_result,queryproteins)
% give protein abundance in mg/gDW and in the order as the proteinss
usages = find(~cellfun(@isempty,strfind(model.rxns,'prot_'))); % pool doesn't start with prot


usageProtein = strrep(model.rxns(usages),'prot_','');
[~,idx1] = ismember(usageProtein,model.enzymes);
usageMW = model.MWs(idx1);
usages = sol_result(usages,1);


enzUsage_abun = usages.*usageMW;  % mmol/gDW * MW = mg/gDW

[~,idx] = ismember(queryproteins,usageProtein);

enzUsage_abun(idx~=0) = enzUsage_abun(idx(idx~=0));

end