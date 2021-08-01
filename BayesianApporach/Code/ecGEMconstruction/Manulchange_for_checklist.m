% manual curation for checklist
% k.lactis Neurospora_crassa use NADPH
inputpath = '../../Result/ssGEMs'; % this can be downloaded from the GitHub/Yeast-Species-GEMs
current_path = pwd;
cd(inputpath)
strains = {'Kluyveromyces_lactis','Neurospora_crassa'};
for i = 1:length(strains)
    load([strains{i},'.mat'])
    model = reducedModel;
    if i == 2
        model = addrxnBack(model,model_original,{'r_0773'},{'Neurospora_crassa@Seq_4559'});
    end
    rxnformula = cell2mat(printRxnFormula(model,'rxnAbbrList','r_0773','metNameFlag',true)); % r_0773 is the alternative complex I
    rxnformula = strrep(rxnformula,'NADH [','NADPH [');
    rxnformula = strrep(rxnformula,'NAD [','NADP(+) [');
    model = changerxn(model,'r_0773',rxnformula);
    reducedModel = model;
    save([strains{i},'.mat'],'reducedModel')
end

%% K.lactis accept NADPH for r_0486
strains = {'Kluyveromyces_lactis'};
for i = 1:length(strains)
    load([strains{i},'.mat'])
    model = reducedModel;
    rxnformula = cell2mat(printRxnFormula(model,'rxnAbbrList','r_0486','metNameFlag',true)); % r_0773 is the alternative complex I
    rxnformula = strrep(rxnformula,'NADH [','NADPH [');
    rxnformula = strrep(rxnformula,'NAD [','NADP(+) [');
    rxnformula = strrep(rxnformula,'glyceraldehyde 3-phosphate','glyceraldehyde&3-phosphate');
    rxnformula = strrep(rxnformula,'<=>','=>');
    cd(current_path)
    model = changerxn(model,'r_0486',rxnformula); 
    [~,idx] = ismember('r_0486',model.rxns);
    model.rxns(idx) = {'r_6000'};
    model.rxnNames{idx} = [model.rxnNames{idx},'NADP'];
    model = addrxnBack(model,model_original,{'r_0486'},model.grRules(idx));
    reducedModel = model;
    cd(inputpath)
    save([strains{i},'.mat'],'reducedModel')
end

%% Ogataea_parapolymorpha and Ogataea_polymorpha have alternative mito NDH1
cd(inputpath)
strains = {'Ogataea_polymorpha','Ogataea_parapolymorpha','Komagataella_pastoris'}; % for Komagataella_pastoris@Seq_333 is used with the published latest model to search PAS_chr3_0792
GPR = {'Ogataea_polymorpha@Seq_1455 or Ogataea_polymorpha@Seq_1941','Ogataea_parapolymorpha@Seq_624 or Ogataea_parapolymorpha@Seq_1666','Komagataella_pastoris@Seq_333 or Komagataella_pastoris@Seq_4347'};
for i = 1:length(strains)
    load([strains{i},'.mat'])
    model = reducedModel;
    cd(current_path)
    model = addrxnBack(model,model_original,{'r_0773'},GPR(i));
    [~,idx] = ismember(split(GPR(i),' or '),model.genes);
    model.proteins(idx) = {'YMR145C'};
    reducedModel = model;
    cd(inputpath)
    save([strains{i},'.mat'],'reducedModel')
end

%% add lactose exchange for k.lcatis
strains = {'Kluyveromyces_lactis','Kluyveromyces_marxianus'};
for i = 1:length(strains)
    load([strains{i},'.mat'])
    model = reducedModel;
    model = addExchangeRxn(model,'m_4652[e]');
    model.rxnNames(strcmp(model.rxns,'EX_m_4652[e]')) = {'lactose exchange'};
    reducedModel = model;
    save([strains{i},'.mat'],'reducedModel')
end
cd(current_path)