function [model,sol_result] = getflux_ec(strain,chemostat,media,anox,carbon,inputpath,state)

if ~strcmp(chemostat,'max')
    dilutionrate = str2double(chemostat);
    chemostat = 'chemo';
elseif strcmp(chemostat,'max')
    chemostat = 'batch';
end

cd(inputpath)

if strcmp(state,'initial')
    load([strain,'_gecko.mat'])
else
    load([strain,'_gecko_improved.mat'])
end

model = ecModel_batch;
if strcmp(anox,'anaerobic')
    anox = 1;
else
    anox = 0;
end
[model,~] = changeMedia(model,carbon,media,anox);
if strcmpi(chemostat,'batch')
        model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'biomass pseudoreaction')),1000,'u');
    model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'biomass pseudoreaction')),0,'l');

    model = changeObjective(model,model.rxns(strcmp(model.rxnNames,'growth')),1);
    try
        sol = optimizeCbModel(model,'max','one');
    catch
        disp('skip this one')
        sol.x = zeros(length(model.rxns),1);
        sol.f = nan;
    end
else
    [~,idx] = ismember([carbon,' exchange'],model.rxnNames);
    model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'biomass pseudoreaction')),1000,'u');
    model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'biomass pseudoreaction')),0,'l');

    model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'growth')|strcmp(model.rxnNames,'biomass exchange')),dilutionrate*0.999,'l');
    model = changeRxnBounds(model,model.rxns(strcmp(model.rxnNames,'growth')|strcmp(model.rxnNames,'biomass exchange')),dilutionrate,'u');
   model = changeObjective(model,model.rxns(idx),1);
    try
        sol = optimizeCbModel(model,'max','one');
    catch
        disp('skip this one, due to no solution')
        sol.x = zeros(length(model.rxns),1);
        sol.f = nan;
    end
    if isempty(sol.f) | isnan(sol.f)
        sol = solveLP(model,1);
        disp('using raven')
    end
end
sol_result = sol.x;
end