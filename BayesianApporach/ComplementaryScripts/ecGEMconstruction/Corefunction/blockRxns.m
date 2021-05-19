%% blockRxns 
function model = blockRxns(model)
% Some reactions should be block to avoid weird flux distributions.

% block Gcy1, an alternative glycerol dissimilation pathway that is active 
% under microaerobic conditions (PMID: 22979944)
idx = find(startsWith(model.rxns,'r_0487_'));
model = changeRxnBounds(model,model.rxns(idx),0,'b');

% block newly added isozyme
%model = changeRxnBounds(model,'r_0438_5',0,'b'); % ferrocytochrome-c:oxygen oxidoreductase

% block newly added isozyme for enolase as Eno1 and Eno2 show higher protein levels
% model = changeRxnBounds(model,'r_0366_1_fwd',0,'b'); % enolase
% model = changeRxnBounds(model,'r_0366_1_rvs',0,'b'); % enolase
% model = changeRxnBounds(model,'r_0366_4_fwd',0,'b'); % enolase
% model = changeRxnBounds(model,'r_0366_4_rvs',0,'b'); % enolase
% model = changeRxnBounds(model,'r_0366_5_fwd',0,'b'); % enolase
% model = changeRxnBounds(model,'r_0366_5_rvs',0,'b'); % enolase

% model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK
idx = regexp(cellstr(model.rxns),'r_4262[\w*]+rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b'); % citrate hydroxymutase

% block some reactions that done in PMID: 28779005.
% idx = regexp(cellstr(model.rxns),'r_2045[\w*]+rvs');
% idx = find(cellfun(@isempty,idx)==0);
% model = changeRxnBounds(model,model.rxns(idx),0,'b'); % serine transport from [m] to [c]

idx = find(startsWith(model.rxns,'r_0659_'));
model = changeRxnBounds(model,model.rxns(idx),0,'b'); % isocitrate dehydrogenase (NADP)

idx = regexp(cellstr(model.rxns),'r_0725[\w*]fwd');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b');% methenyltetrahydrofolate cyclohydrolase

idx = find(startsWith(model.rxns,'r_0918'));
model = changeRxnBounds(model,model.rxns(idx),0,'b'); % phosphoserine transaminase

idx = find(startsWith(model.rxns,'r_4235'));
model = changeRxnBounds(model,model.rxns(idx),0,'b'); % weird reaction from glc to g6p

idx = regexp(cellstr(model.rxns),'r_4216[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b'); % block the reaction to produce FMN without ATP

idx = regexp(cellstr(model.rxns),'r_4264[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b');  % succinate:NAD+ oxidoreductase

% The following two reactions account for the transport of pyruvate from
% [c] to [m] without enzyme cost, should be blocked.
model = changeRxnBounds(model,'r_1137',0,'b');
model = changeRxnBounds(model,'r_1138',0,'b');

% Block backward reactions for PLP synthesis.
idx = regexp(cellstr(model.rxns),'r_4211[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b');
idx = regexp(cellstr(model.rxns),'r_4212[\w*]rvs');
idx = find(cellfun(@isempty,idx)==0);
model = changeRxnBounds(model,model.rxns(idx),0,'b');

% for not producing urea as byproduct
printRxnFormula(model,'rxnAbbrList','r_4790','metNameFlag',true);
[~,idx] = ismember('r_4790',model.rxns);
model.lb(idx(idx~=0)) = 0;