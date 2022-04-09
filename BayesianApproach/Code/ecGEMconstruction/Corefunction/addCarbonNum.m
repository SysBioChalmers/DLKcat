function model = addCarbonNum(model)

% get carbonnum for each exchange rxn to further calculation of error
    model.excarbon = zeros(1,length(model.rxns));
    [EXrxn, EXrxnIdx] = getExchangeRxns(model);
    [CarbonNum,EXfors] = getcarbonnum(model,EXrxn);
    model.excarbon(EXrxnIdx) = CarbonNum;
    model.excarbon(strcmp(model.rxnNames,'growth')) = 41;
end
