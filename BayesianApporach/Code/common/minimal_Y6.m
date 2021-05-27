function model = minimal_Y6(model)
% change Y6 model media to minimal - ammonium, glucose, oxygen,
% phosphate, sulphate
% Bicarbonate production is blocked to get rid of conflict due to bicarbonate and carbon dioxide
% are regarded as equivalent
% the function is from:https://doi.org/10.1371/journal.pcbi.1004530

% start with a clean slate: set all exchange reactions to upper bound = 1000
% and lower bound = 0 (ie, unconstrained excretion, no uptake)

%Feiran Li 2018.09.05
protein_Content = model.lb(strcmp(model.rxns,{'EX_protein_pool'}));
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

model.lb(strcmp(model.rxns,{'EX_protein_pool'})) = protein_Content;

desiredExchanges = {'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };   % potassium exchange

blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange

glucoseExchange = {'r_1714'};     % D-glucose exchange

uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges);
glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
BlockedRxnIndex      = findRxnIDs(model,blockedExchanges);

if length(find(uptakeRxnIndexes~= 0)) ~= 15
    warning('Not all exchange reactions were found.')
end

model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0))     = -1000;
model.lb(glucoseExchangeIndex) = -1000;

model.lb(BlockedRxnIndex) = 0;
model.ub(BlockedRxnIndex) = 0;

end
