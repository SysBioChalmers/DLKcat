function sol_c= simulatecarbonsource(model,carbonsources,state,dilution)

if strcmpi(state,'batch')
    for i = 1:length(carbonsources)
        model_tmp = model;
        [~,idx] = ismember(strcat(carbonsources{i},' exchange'),model_tmp.rxnNames);
        model_tmp = changeRxnBounds(model_tmp,'r_1714',0,'l');
        model_tmp = changeRxnBounds(model_tmp,model_tmp.rxns(idx),-1000,'l');
        model_tmp = changeObjective(model_tmp,model_tmp.rxns(strcmp(model_tmp.rxnNames,'growth')),1);
        sol = optimizeCbModel(model_tmp,'max','one');
        sol_c(:,i) = sol.x;
    end
elseif  strcmpi(state,'chemostat')
    for i = 1:length(carbonsources)
        model_tmp = model;
        [~,idx] = ismember(strcat(carbonsources{i},' exchange'),model_tmp.rxnNames);
        model_tmp = changeRxnBounds(model_tmp,'r_1714',0,'l');
        model_tmp = changeRxnBounds(model_tmp,model_tmp.rxns(idx),-1000,'l');
        model_tmp = changeRxnBounds(model_tmp,model_tmp.rxns(strcmp(model_tmp.rxnNames,'growth')|strcmp(model_tmp.rxnNames,'biomass exchange')),dilution);
        model_tmp = changeObjective(model_tmp,model_tmp.rxns(idx),1);
        sol = optimizeCbModel(model_tmp,'max','one');
        if isempty(sol.f) | isnan(sol.f)
            sol = solveLP(model_tmp,1);
            disp('using raven')
        end
        sol_c(:,i) = sol.x;
    end
end
end

