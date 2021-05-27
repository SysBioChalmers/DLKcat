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
    pos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');%NGAM
    %model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub to be mu
    model.lb(pos) = NGAM;
    model.ub(pos) = NGAM;
    
end

end
