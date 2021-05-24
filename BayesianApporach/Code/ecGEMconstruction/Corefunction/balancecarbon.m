function [biomass_carbon] = balancecarbon(model,growthdata)

% This function is to scale growthdata to fix the carbon recovery to be 1.
[EXrxn, EXrxnIdx] = getExchangeRxns(model);
[CarbonNum,~,changedrxns] = getcarbonnum(model,EXrxn);
model.excarbon(EXrxnIdx,1) = CarbonNum;
model.excarbon(ismember(model.rxns,changedrxns)) = 0;
disp('carbon has been changed back to zeros:')
printRxnFormula(model,'rxnAbbrList',changedrxns,'metNameFlag',true);
% perform one simulation and calculate the biomass carbon
model = changeRxnBounds(model,'r_2111',0.1,'b'); % biomass
model = changeRxnBounds(model,'r_1714',-1000,'l'); % glc
sol = optimizeCbModel(model,'max','one');
sol_ex = sol.x(EXrxnIdx);
biomass_carbon = abs(sum(sol_ex.*model.excarbon(EXrxnIdx))/0.1);
if abs(biomass_carbon-40) > 5
    warning('biomass carbon number may be wrong, please check this')
end

model.excarbon(strcmp(model.rxnNames,'growth')) = biomass_carbon;
exp = cell2mat(growthdata(:,3:11));   
exp_new = exp.*[1,-1,1,1,1,1,1,1,-1];

for j = 1:length(growthdata(:,1))
    if numel(find(isnan(exp_new(j,:)))) == 0
        ex_mets = {'growth',[growthdata{j,2},' exchange'],'acetate exchange','ethanol exchange','glycerol exchange','pyruvate exchange','ethyl acetate exchange','carbon dioxide exchange','oxygen exchange'};
        [~,idx] = ismember(ex_mets,model.rxnNames);
        excarbon = model.excarbon(idx);
        leftout = sum(exp_new(j,:).*excarbon');
        if abs(leftout)/abs((exp_new(j,2).*excarbon(2)')) > 0.08
            warning(['check the carbon recovery for this growthdata:',num2str(j),', leftover:',num2str(abs(leftout)/abs((exp_new(j,2).*excarbon(2)')))])
            exp(j,end-1) = exp(j,end-1) + leftout;
            exp(j,end) = exp(j,end)/exp(j,end-1)*(exp(j,end-1) + leftout); % scale O2 with CO2
        else
            disp(['pass:',num2str(j),', leftover:',num2str(leftout)])
        end
     
    end
end
            
        
end
