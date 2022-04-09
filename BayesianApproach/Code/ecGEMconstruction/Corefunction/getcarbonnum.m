function [CarbonNum,EXfors,changedmets] = getcarbonnum(model,exrxn)
% exrxn must be an exchange rxn and only contain one mets in the equation
[~,idx] = ismember(exrxn,model.rxns);
EXmets = model.S(:,idx);
% find all exchange mets
EXmetsIdx = zeros(length(exrxn),1);
for k = 1:length(EXmets(1,:))
    EXmetsIdx(k) = find(EXmets(:,k)); % 
end
EXfors = model.metFormulas(EXmetsIdx);
[Ematrix, elements] = getElementalComposition(EXfors,{'C'});
Ematrix = Ematrix(:,1);
changedmets = exrxn(isnan(Ematrix));
Ematrix(isnan(Ematrix)) = 1;
CarbonNum = Ematrix';
end