%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eModel = convertToEnzymeModel(model,Genes,uniprots,kcats)
% Converts standard GEM to GEM accounting for enzymes as pseudo
% metabolites, with -(1/kcat) as the corresponding stoich. coeffs.
%
% INPUT:
% model             The GEM structure (1x1 struct)
% enzymedata        structure which contains kcat and MW fro each enzyme
% complex and also for subunit info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eModel = convertToGeckoModel(model,enzymedata,total_protein)

% expand matrix to include the enzymes in the model
model.id = ['emodel_',strrep(model.id,' specific model genereted from panYeast','')];
model = minimal_Y6(model);
nmet = length(model.mets);
nrxn = length(model.rxns);
ngene = length(model.genes);
model.mets(nmet+1:nmet+ngene,1) = strcat(model.genes,repmat({' [c]'},ngene,1));
model.rxns(nrxn+1:nrxn+ngene,1) = strcat(repmat({'prot_'},ngene,1),model.genes);
model.mets(nmet+ngene+1,1) = {'protein_pool [c]'};
model.rxns(nrxn+ngene+1,1) = {'EX_protein_pool'};
model.S(nmet+ngene+1,nrxn+ngene+1) = -1;
model.lb(nrxn+1:nrxn+ngene) = 0;
model.lb(nrxn+ngene+1) = -total_protein*1000;
model.ub(nrxn+1:nrxn+ngene+1) = 1000;
% add enzyme into metabolic rxn
for i = 1:length(enzymedata.rxn_list)
   [~,idx] = ismember(enzymedata.rxn_list(i),model.rxns);
   subunitlist = enzymedata.subunit(i,:);
   subunit_num = cell2mat(cellfun(@any,subunitlist,'UniformOutput',false));
   subunitlist = subunitlist(1,subunit_num);
   subunitlist_tmp = strcat(subunitlist',repmat({' [c]'},sum(subunit_num),1));
   subunitcoef = enzymedata.subunit_stoichiometry(i,subunit_num);
   [~,idx2] = ismember(subunitlist_tmp,model.mets);
   model.S(idx2,idx) = -subunitcoef./enzymedata.kcat(i);
end

% add all coef for enzyme conversion from protein pool
for i = 1:length(enzymedata.proteins)
    subunit_rxn = strcat('prot_',enzymedata.proteins{i});
    subunit_met = strcat(enzymedata.proteins{i},' [c]');
    [~,idx] = ismember(subunit_rxn,model.rxns);
    [~,idx2] = ismember(subunit_met,model.mets);
    if idx~=0
    model.S(nmet+ngene+1,idx) = -enzymedata.proteinMW(i);
    model.S(idx2,idx) = 1;
    end
end

% fix other field
model.enzymes = enzymedata.proteins;
model.MWs = enzymedata.proteinMW;
model.rxnNames(nrxn+1:nrxn+ngene+1) = model.rxns(nrxn+1:nrxn+ngene+1);
model.rxnECNumbers(nrxn+ngene+1) = {''};
model.rxnMetaNetXID(nrxn+ngene+1) = {''};
model.rxnECNumbers(nrxn+ngene+1) = {''};
model.rxnNotes(nrxn+ngene+1) = {''};
model.rxnKEGGID(nrxn+ngene+1) = {''};
model.rxnConfidenceScores(nrxn+ngene+1) = 0;

model.metChEBIID(nmet+ngene+1) = {''};
model.metFormulas(nmet+ngene+1) = {''};
model.metKEGGID(nmet+ngene+1) = {''};
model.b(nmet+ngene+1) = 0;
model.c(nrxn+ngene+1) = 0;
model.grRules(nrxn+ngene+1) = {''};
model.rules(nrxn+ngene+1) = {''};
model.metNames(nmet+1:nmet+ngene+1) = model.mets(nmet+1:nmet+ngene+1);
model.metCharges(nmet+ngene+1) = 0;
model.metNotes(nmet+ngene+1) = {''};
if isfield(model,'rxnGeneMat')
    model = rmfield(model,'rxnGeneMat');
end
model.csense = repmat('E',length(model.mets),1);
[m,n] = size(model.rxnNotes);
if m < n
    model.rxnNotes = model.rxnNotes';
end
[m,n] = size(model.rxnConfidenceScores);
if m < n
    model.rxnConfidenceScores = model.rxnConfidenceScores';
end
eModel = model;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
