%% splitModel 
function new_model = splitModel(model)

disp('Splitting reactions...');
% Change gene ID in the model.rules into that in the model.genes.
for i = 1:length(model.genes)
    old = strcat('x(',num2str(i),')');
    new = model.genes{i};
    model.rules = cellfun(@(x) strrep(x,old,new),...
                            model.rules,'UniformOutput',false);
end

% Generate a matrix for splitting isozymes.
matrix1 = struct();
matrix1.RxnList = cell(0,1);
matrix1.CompoList = cell(0,1);
matrix1.CoeffList = zeros(0,1);
matrix1.CatalystList = cell(0,1);
matrix1.LBList = zeros(0,1);
matrix1.UBList = zeros(0,1);

% Assuming there are the following types of GPRs in the M model.
% 'geneA'; '(geneA or geneB or ...)'; '(geneA and geneB and ...)';
% '((geneA and geneB and ...) or (geneC and geneD and ...) or ...)'.

% Main loop
for i = 1:length(model.rxns)

    idx_substrate = model.S(:,i) < 0;
    idx_product = model.S(:,i) > 0;
    CompoS = model.mets(idx_substrate);
    CompoP = model.mets(idx_product);
    CoeffS = model.S(idx_substrate,i);
    CoeffP = model.S(idx_product,i);
    CompoList = [CompoS;CompoP];
    CoeffList = [CoeffS;CoeffP];
    n = length(CompoList);
    Rxn = model.rxns(i);
    Catalyst = model.rules(i);
    lb_tmp = model.lb(i);
    ub_tmp = model.ub(i);
    RxnList = repmat(Rxn,n,1);
    CatalystList = repmat(Catalyst,n,1);
    LBList = repmat(lb_tmp,n,1);
    UBList = repmat(ub_tmp,n,1);
    
    z = strtrim(model.rules{i});%get GPR and delete leading and trailing whitespace
    x = length(matrix1.RxnList);
    if isempty(z)%if is spontaneous reaction
        matrix1.RxnList(x+1:x+n,1) = RxnList;
        matrix1.CompoList(x+1:x+n,1) = CompoList;
        matrix1.CoeffList(x+1:x+n,1) = CoeffList;
        matrix1.CatalystList(x+1:x+n,1) = CatalystList;
        matrix1.LBList(x+1:x+n,1) = LBList;
        matrix1.UBList(x+1:x+n,1) = UBList;
    else%if is enzymatic reaction
        z = strrep(z,'(','');
        z = strrep(z,')','');
        z = strtrim(z);%delete leading and trailing whitespace if it has
        if ~contains(z,' | ') && ~contains(z,' & ')%if has a single enzyme
            matrix1.RxnList(x+1:x+n,1) = RxnList;
            matrix1.CompoList(x+1:x+n,1) = CompoList;
            matrix1.CoeffList(x+1:x+n,1) = CoeffList;
            matrix1.CatalystList(x+1:x+n,1) = CatalystList;
            matrix1.LBList(x+1:x+n,1) = LBList;
            matrix1.UBList(x+1:x+n,1) = UBList;
        else%if has multiple genes, either isozymes or complexes
            if ~contains(z,' | ')%if only has complex
                matrix1.RxnList(x+1:x+n,1) = RxnList;
                matrix1.CompoList(x+1:x+n,1) = CompoList;
                matrix1.CoeffList(x+1:x+n,1) = CoeffList;
                matrix1.CatalystList(x+1:x+n,1) = repmat({z},n,1);
            	matrix1.LBList(x+1:x+n,1) = LBList;
                matrix1.UBList(x+1:x+n,1) = UBList;
            else%if has isozymes, and has or does not have complex
                if ~contains(z,' & ')%if only has isozymes
                    z = strsplit(z,' | ');%split isozymes
                    for j = 1:length(z)%change reaction ID
                        rxnname_tmp = strcat(model.rxns{i},'_',num2str(j));
                        catalyst_tmp = z{j};
                        catalyst_tmp = strtrim(catalyst_tmp);
                        matrix1.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
                        matrix1.CompoList(x+1:x+n,1) = CompoList;
                        matrix1.CoeffList(x+1:x+n,1) = CoeffList;
                        matrix1.CatalystList(x+1:x+n,1) = repmat({catalyst_tmp},n,1);
                    	matrix1.LBList(x+1:x+n,1) = LBList;
                        matrix1.UBList(x+1:x+n,1) = UBList;
                        x = length(matrix1.RxnList);
                    end
                else%if has both isozymes and complex
                    if contains(z,' | ')
                    % if GPR is made up of iso-complexes, e.g., '( x(5) & x(2) ) | ( x(5) & x(4) )'
                        z = strsplit(z,' | ');%split isozymes
                        for j = 1:length(z)%change reaction ID
                            rxnname_tmp = strcat(model.rxns{i},'_',num2str(j));
                            catalyst_tmp = z{j};
                            catalyst_tmp = strtrim(catalyst_tmp);
                            matrix1.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
                            matrix1.CompoList(x+1:x+n,1) = CompoList;
                            matrix1.CoeffList(x+1:x+n,1) = CoeffList;
                            matrix1.CatalystList(x+1:x+n,1) = repmat({catalyst_tmp},n,1);
                        	matrix1.LBList(x+1:x+n,1) = LBList;
                            matrix1.UBList(x+1:x+n,1) = UBList;
                            x = length(matrix1.RxnList);
                        end
                    else
                        error('There are complicated GPRs in M model.')
                    end
                end
            end
        end
    end
end

% Generate a matrix for splitting reversible reactions.
matrix2 = struct();
matrix2.RxnList = cell(0,1);
matrix2.CompoList = cell(0,1);
matrix2.CoeffList = zeros(0,1);
matrix2.CatalystList = cell(0,1);
matrix2.LBList = zeros(0,1);
matrix2.UBList = zeros(0,1);

UnqRxnList = unique(matrix1.RxnList);
for i = 1:length(UnqRxnList)
    id_tmp = UnqRxnList(i);
    idx = ismember(matrix1.RxnList,id_tmp);
    
    RxnList = matrix1.RxnList(idx);
    CompoList = matrix1.CompoList(idx);
    CoeffList = matrix1.CoeffList(idx);
    CatalystList = matrix1.CatalystList(idx);
    LBList = matrix1.LBList(idx);
    UBList = matrix1.UBList(idx);
    n = length(RxnList);
    
    x = length(matrix2.RxnList);%count rows in matrix2
    if isempty(CatalystList{1})%if is spontaneous reaction
        matrix2.RxnList(x+1:x+n,1) = RxnList;
        matrix2.CompoList(x+1:x+n,1) = CompoList;
        matrix2.CoeffList(x+1:x+n,1) = CoeffList;
        matrix2.CatalystList(x+1:x+n,1) = CatalystList;
        matrix2.LBList(x+1:x+n,1) = LBList;
        matrix2.UBList(x+1:x+n,1) = UBList;
    else%if is enzymatic reaction
        if LBList(1) >= 0%if is not reversible
            matrix2.RxnList(x+1:x+n,1) = RxnList;
            matrix2.CompoList(x+1:x+n,1) = CompoList;
            matrix2.CoeffList(x+1:x+n,1) = CoeffList;
            matrix2.CatalystList(x+1:x+n,1) = CatalystList;
            matrix2.LBList(x+1:x+n,1) = LBList;
            matrix2.UBList(x+1:x+n,1) = UBList;
        else%if is reversible
            %add forward rows
            if UBList(1) > 0
                rxnname_tmp = strcat(RxnList{1},'_fwd');
                matrix2.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
                matrix2.CompoList(x+1:x+n,1) = CompoList;
                matrix2.CoeffList(x+1:x+n,1) = CoeffList;
                matrix2.CatalystList(x+1:x+n,1) = CatalystList;
                matrix2.LBList(x+1:x+n,1) = zeros(n,1);
                matrix2.UBList(x+1:x+n,1) = UBList;
            end
            x = length(matrix2.RxnList);
            %add reverse rows
            rxnname_tmp = strcat(RxnList{1},'_rvs');
            matrix2.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
            matrix2.CompoList(x+1:x+n,1) = CompoList;
            matrix2.CoeffList(x+1:x+n,1) = -1*CoeffList;
            matrix2.CatalystList(x+1:x+n,1) = CatalystList;
            matrix2.LBList(x+1:x+n,1) = zeros(n,1);
            matrix2.UBList(x+1:x+n,1) = -1*LBList;
        end
    end
end

matrixSplit = matrix2;

% Converted to COBRA model
new_model = struct();
new_model.rxns = cell(0,1);
new_model.lb = zeros(0,1);
new_model.ub = zeros(0,1);
new_model.mets = cell(0,1);
new_model.S = sparse(0,0);
new_model.b = zeros(0,1);
new_model.rxnGeneMat = sparse(0,0);
new_model.c = zeros(0,1);
new_model.rules = cell(0,1);
new_model.genes = cell(0,1);
new_model.csense = char();

new_model.metFormulas = cell(0,1);
new_model.metNames = cell(0,1);
new_model.metKEGGID = cell(0,1);
new_model.metChEBIID = cell(0,1);
new_model.rxnNames = cell(0,1);

new_model.osenseStr = model.osenseStr;

if isfield(model,'compNames')
    new_model.compNames = model.compNames;
end
if isfield(model,'comps')
    new_model.comps = model.comps;
end

UnqRxnList = unique(matrixSplit.RxnList);
for i = 1:length(UnqRxnList)
    
    id_tmp = UnqRxnList(i);
    idx = ismember(matrixSplit.RxnList,id_tmp);
    
    RxnList = matrixSplit.RxnList(idx);
    CompoList = matrixSplit.CompoList(idx);
    CoeffList = matrixSplit.CoeffList(idx);
    CatalystList = matrixSplit.CatalystList(idx);
    LBList = matrixSplit.LBList(idx);
    UBList = matrixSplit.UBList(idx);
    
    rxn_tmp = RxnList{1};
    catalyst_tmp = CatalystList{1};
    lb_tmp = unique(LBList);
    ub_tmp = unique(UBList);
    
    % add metabolites
    [isInModelidx,~] = ismember(CompoList,new_model.mets);
    NotInModel = CompoList(~isInModelidx);
    idx_met_tmp = ismember(model.mets,NotInModel);
	metid_tmp = model.mets(idx_met_tmp);
	metname_tmp = model.metNames(idx_met_tmp);
	metformula_tmp = model.metFormulas(idx_met_tmp);
	ChEBIiD_tmp = model.metChEBIID(idx_met_tmp);
	KEGGid_tmp = model.metKEGGID(idx_met_tmp);
    new_model = addMetabolite(new_model,metid_tmp,...
                              'metName',metname_tmp,...
                              'metFormula',metformula_tmp,...
                              'ChEBIID',ChEBIiD_tmp,...
                              'KEGGId',KEGGid_tmp);
    
    % add reactions
    reactionname_tmp = model.rxnNames{ismember(model.rxns,rxn_tmp(1:6))};
    new_model = addReaction(new_model,rxn_tmp,...
                            'reactionName',reactionname_tmp,...
                            'metaboliteList',CompoList,...
                            'stoichCoeffList',CoeffList,...
                            'lowerBound',lb_tmp,...
                            'upperBound',ub_tmp,...
                            'geneRule',catalyst_tmp);
end

% add proteins
[~,idx] = ismember(new_model.genes,model.genes);
new_model.proteins = model.proteins(idx);
new_model = changeObjective(new_model,'r_2111');
new_model.id = model.id;
% sort protein MNX IDs
[~,idx] = ismember(new_model.mets,model.mets);
new_model.metMetaNetXID = model.metMetaNetXID(idx);
% new_model.subSystem = repelem({'Metabolism'},length(new_model.rxns))';