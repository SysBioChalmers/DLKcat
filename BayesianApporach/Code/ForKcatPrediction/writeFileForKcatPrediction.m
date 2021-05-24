function writeFileForKcatPrediction(species,inputpath)
% load('strains.mat') % strains
% species = strains;
% inputpath = 'downloaded model from https://github.com/SysBioChalmers/Yeast-Species-GEMs/tree/master/Reconstruction_script/ModelFiles/xml'
current_path = pwd;

load('MNXprop.mat') % generate from MNXpropExtract.m
load('conserveMets.mat')
% remove conservedmets from prediction substrates

tmp = {'Ala-tRNA(Ala);Arg-tRNA(Arg);Arg-tRNA(Arg);Asn-tRNA(Asn);Asn-tRNA(Asn);Asp-tRNA(Asp);Asp-tRNA(Asp);Cys-tRNA(Cys);fMet-tRNA(fMet);Gln-tRNA(Gln);Glu-tRNA(Glu);Glu-tRNA(Glu);Gly-tRNA(Gly);His-tRNA(His);His-tRNA(His);Ile-tRNA(Ile);Ile-tRNA(Ile);Leu-tRNA(Leu);Leu-tRNA(Leu);Lys-tRNA(Lys);Lys-tRNA(Lys);Met-tRNA(Met);Met-tRNA(Met);Phe-tRNA(Phe);Phe-tRNA(Phe);Pro-tRNA(Pro);Ser-tRNA(Ser);Thr-tRNA(Thr);Thr-tRNA(Thr);Trp-tRNA(Trp);Trp-tRNA(Trp);Tyr-tRNA(Tyr);Tyr-tRNA(Tyr);Val-tRNA(Val);Val-tRNA(Val);tRNA(Ala);tRNA(Arg);tRNA(Arg);tRNA(Asn);tRNA(Asn);tRNA(Asp);tRNA(Asp);tRNA(Cys);tRNA(Gln);tRNA(Glu);tRNA(Glu);tRNA(Gly);tRNA(His);tRNA(His);tRNA(Ile);tRNA(Ile);tRNA(Leu);tRNA(Leu);tRNA(Lys);tRNA(Lys);tRNA(Met);tRNA(Met);tRNA(Phe);tRNA(Phe);tRNA(Pro);tRNA(Ser);tRNA(Thr);tRNA(Thr);tRNA(Trp);tRNA(Trp);tRNA(Tyr);tRNA(Tyr);tRNA(Val);tRNA(Val);L-glutamyl-tRNA(Gln);Gln-tRNA(Gln);tRNA(Pro);Pro-tRNA(Pro);Glycyl-tRNA(Ala);D-tyrosyl-tRNA(Tyr);H2O;H+'};
tmp = split(tmp,';');
conserveMets.currency = tmp;

for i = 1:length(species)
    disp([num2str(i),'/',length(species),' for species',species{i}])
    cd(inputpath)
    load([species{i},'.mat'])
    cd(current_path)
    model = reducedModel;
    model = splitModel(model);
    save([species{i},'_split.mat'],'model')
    modelr = ravenCobraWrapper(model);
    a = cell(length(model.rxns),7);
    conservedidx_pair = find(ismember(modelr.metNames,conserveMets.pair));
    conservedidx_currency = find(ismember(modelr.metNames,conserveMets.currency));
    % map to the yeastGEM_MNX
    [~,idx] = ismember(modelr.mets,MNXYeast(:,2));
    model.metMetaNetXID(idx~=0) = MNXYeast(idx(idx~=0),1);
    
    % map to missingID
    for m = 1:length(missingID(:,1))
        idx = find(ismember(modelr.metNames,missingID(m,1)));
        model.metMetaNetXID(idx(idx~=0)) = missingID(m,2);
    end
    model.metSmiles = cell(length(model.mets),1);
    model.metMNXname = cell(length(model.mets),1);
    model.metMNXref = cell(length(model.mets),1);
    [~,idx] = ismember(model.metMetaNetXID,MNXdepr(:,1));
    model.metMetaNetXID(idx~=0) = MNXdepr(idx(idx~=0),2);
    [~,idx] = ismember(model.metMetaNetXID,MNXprop(:,1));
    model.metSmiles(idx~=0) = MNXprop(idx(idx~=0),9);
    model.metMNXname(idx~=0) = MNXprop(idx(idx~=0),2);
    model.metMNXref(idx~=0) = MNXprop(idx(idx~=0),3);
    model.metSmiles(cellfun(@isempty,model.metSmiles)) = {''};
    model.metMNXname(cellfun(@isempty,model.metMNXname)) = {''};
    model.metMNXref(cellfun(@isempty,model.metMNXref)) = {''};
    for j = 1:length(model.rxns)
        a(j,1) = model.rxns(j);
        mets = find(full(model.S(:,j)) < 0);
%       mets_prods = find(full(model.S(:,j)));
%         if length(intersect(mets_prods,mets)) >= 2
%             mets = setdiff(mets,conservedidx_pair);% only remove conservedMets from sub if there are paired for example NADH and NAD, then it can save kcats prediction for NAD synthesis pathway 
%         end
         mets = setdiff(mets,conservedidx_currency); % remove the currency mets from the prediction
        metsMNX = model.metMetaNetXID(mets);
        metsName = modelr.metNames(mets); % useing the raven model metname without comp
        metsSmiles = model.metSmiles(mets);
        metsMNXname = model.metMNXname(mets);
        metsMNXref = model.metMNXref(mets);
        a(j,2) = join(metsMNX,';');
        a(j,3) = join(metsName,';');
        a(j,4) = join(metsMNXname,';');
        a(j,5) = join(metsSmiles,';');
        a(j,6) = join(model.genes(find(model.rxnGeneMat(j,:))),';');
        a(j,7) = join(metsMNXref,';');

    end
    a(cellfun(@isempty,a(:,6)),:) = [];
    columnID = {'# rxnID','MNXID','met_model','met_standard_name','Smiles','genes','refID'};
    writetable(cell2table(a,'VariableNames',columnID),[species{i},'ForKcatPrediction.txt'],'Delimiter','\t')
i
end
end