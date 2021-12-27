%% fitGAM
function [model,GAM,NGAM] = fitGAM(model,growthdata,condition)

if nargin < 3
    condition = 'aerobic';
end

% filter chemostat data

idx = find(strcmpi(growthdata(:,14),condition) & strcmpi(growthdata(:,2),'(S)-lactate'));
exp_GAM = cell2mat(growthdata(idx,3:11)); % u sub ace eth gly pyr co2 o2
exp_GAM = exp_GAM.*[1,-1,1,1,1,1,1,1,-1];
[~,I] = sort(exp_GAM(:,1));
exp_GAM = exp_GAM(I,:);
ATPyield = zeros(length(exp_GAM(:,1)),1);


for i = 1:length(exp_GAM(:,1))
    ex_mets = {'biomass pseudoreaction','(S)-lactate exchange','acetate exchange','ethanol exchange','glycerol exchange','formate exchange','carbon dioxide exchange','oxygen exchange','EX_protein_pool'};
    [~,idx] = ismember(ex_mets,model.rxnNames);
    model_tmp = model;
    bioPos = strcmp(model_tmp.rxnNames,'biomass pseudoreaction');
    for j = 1:length(model_tmp.mets)
        S_ix  = model_tmp.S(j,bioPos);
        isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
            'H+ [cytoplasm]','phosphate [cytoplasm]'},model_tmp.metNames{j})) == 1;
        if S_ix ~= 0 && isGAM
            model_tmp.S(j,bioPos) = 0;
        end
    end
    
    % first minimize carbon source uptake
    model_tmp.lb(strcmp(model_tmp.rxns,model_tmp.rxns(idx(2)))) = exp_GAM(i,2); % not constrain the substrate usage
    model_tmp.lb(strcmp(model_tmp.rxns,model_tmp.rxns(idx(1)))) = exp_GAM(i,1);
    objrxn = strcmpi(model.rxnNames,'non-growth associated maintenance reaction');
    model_tmp = setParam(model_tmp,'lb',model_tmp.rxns(objrxn),0);
    model_tmp = setParam(model_tmp,'ub',model_tmp.rxns(objrxn),1000);
    model_tmp = setParam(model_tmp,'obj',model_tmp.rxns(objrxn),1);
    sol = optimizeCbModel(model_tmp,'max');
    ATPyield(i) = sol.f;
end
d = 1;
try
    r = polyfit(exp_GAM(1:4,1),ATPyield(1:4),d);
    GAM = r(1);
    NGAM = r(2);
    model = changeGAM(model,GAM,NGAM);
catch
    warning('not enough growth data')
    r = polyfit(exp_GAM(:,1),ATPyield(:), 1);
    GAM = r(1);
    NGAM = r(2);
end
end

function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >1
    if NGAM > 0
    pos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');%NGAM
    %model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub to be mu
    model.lb(pos) = NGAM;
    model.ub(pos) = NGAM;
    end
    
end

end
