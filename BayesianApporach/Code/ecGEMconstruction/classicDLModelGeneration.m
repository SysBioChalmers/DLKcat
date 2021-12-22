function classicDLModelGeneration(strains,modelpath,dbpath,kcatpredictionPath)
% This function is to analyze why some cell can grow very fast to 0.8h-1,
% while others can only grow to 0.18h-1
% STEP 1 split reactions reversible reactions and reactions that
% catalyzed by isoenzymes
% STEP 2 collect growth rates phenotype data
% STEP 3 collect growth rates phenotype data

% dbpath = '/Users/feiranl/Documents/GitHub/DLKcat/BayesianApporach/Results/Proteinfasta'
% modelpath = '/Users/feiranl/Documents/GitHub/DLKcat/BayesianApporach/Results/ssGEMs'
% kcatpredictionPath = '/Users/feiranl/Documents/GitHub/DLKcat/BayesianApporach/Results/PredcitedKcat343species'
% load('strains.mat') % strains
%% step 1
% Reformulate the orginal model
%   1. Reactions with isozymes (i.e., 'OR' case in GPR) will be copied,
%   each isozyme will be assigned to each copy.
%   2. Reversible reactions will be split into forward and reverse ones.
current_path = pwd;
mkdir('../../Results/model_build_files/model_classic/')
mkdir('../../Results/model_build_files/model_dl/')
 for m = 1:length(strains)
    strain = strains{m};
    cd(modelpath)
%     model = load([strain,'.mat']);
%     model = model.reducedModel;
%    model_split = splitModel(model);
%    model = model_split;
    model = load([strain,'_split.mat']);
    tmp = fieldnames(model);
    model = eval(['model.',cell2mat(tmp)]);
    cd(current_path)

    disp('finish splitting model')
    
    % read fasta file
    cd(dbpath);
    File = [strain,'.fasta'];
    %File = lower(File);
    [Header, Sequence] = fastaread(File);
    Header = Header';
    cd(current_path)
    
    % calculate MW for the sequence
    [~,idx] = ismember(model.genes,Header);
    clear MWdata
    for i = 1:length(model.genes)
        seq = Sequence{idx(i)};
        seq = strrep(seq,'X','N'); % fix the uncertainties by replacing x to the median mass aa asn
        seq = strrep(seq,'x','n');
        seq = strrep(seq,'*','');
        MWdata.MW(i,1) = getProteinMW(seq);
        MWdata.genes(i,1) = model.genes(i);
    end
    
    model = blockRxns(model);
    model = addCarbonNum(model);
%     % modify two stuff
%     model.S(strcmp(model.mets,'s_0796[e]'),startsWith(model.rxns,'r_1134')) = 0;
%     model.S(strcmp(model.mets,'s_0794[c]'),startsWith(model.rxns,'r_1134')) = 0;
%     model.S(strcmp(model.mets,'s_0796[e]'),startsWith(model.rxns,'r_1139')) = 0;
%     model.S(strcmp(model.mets,'s_0794[c]'),startsWith(model.rxns,'r_1139')) = 0;
    
    disp('add growth data')
    % add growth data
    [~,~,growthrates] = xlsread('growthratedata.xlsx','growthrates');
    % {species,'carbon','µmax(h?1)','rsub(mmolg?1h?1)','racetate(mmolg?1h?1)','rethanol(mmolg?1h?1)','rglycerol(mmol/gDW/h)','rpyruvuate(mmol/gDW/h)','rethyl acetate(mmol/gDW/h)','qCO2(mmolg?1h?1)','qO2(mmolg?1h?1)','temp','ph','aerobic','condition','media'}
     growthrates = growthrates(2:end,[4,7,10,12,13,14,15,16,17,18,19,23,24,25,2,6]);
%     %filter this species at least contains growth rate and substrate uptake
%     %rate, filter that whhether aerobic/anaerobic, filter lower temp than 20
%     idx = ismember(lower(growthrates(:,1)),lower(strrep(strain,'_',' '))) & ~isnan(cell2mat(growthrates(:,3))) & ~isnan(cell2mat(growthrates(:,4))) & ~strcmp(growthrates(:,14),'limited') & (cell2mat(growthrates(:,12))>20|isnan(cell2mat(growthrates(:,12)))) & (cell2mat(growthrates(:,13))>4|isnan(cell2mat(growthrates(:,13))));
%     growthdata = growthrates(idx,:);
%     idx = ismember(lower(growthrates(:,1)),lower(strrep(strain,'_',' ')))  & (strcmpi(growthrates(:,16),'minimal')|strcmpi(growthrates(:,16),'NA'));
%     tmp = growthrates(idx,:);
%     sub = unique(tmp(:,2));
%     max_growth = {};
%     for j = 1:length(sub)
%         idx = find(strcmpi(tmp(:,2),sub(j)));
%         [~,b] = max(cell2mat(tmp(idx,3)));
%         max_growth = [max_growth;tmp(idx(b),:)];
%     end
%     
    
%     % optional
%     exp = cell2mat(growthdata(:,3:11)); % u sub ace eth gly pyr co2 o2
%     testsize = length(exp(:,1)); % only test 30 data
%     if testsize > 40
%         testsize = 40;
%     end
%     [~,I] = mink(sum(isnan(exp),2),testsize);
%     growthdata = growthdata(I,:);
%     
%     
    idx = ismember(lower(growthrates(:,1)),lower(strain))  & ~isnan(cell2mat(growthrates(:,3))) & ~isnan(cell2mat(growthrates(:,4))) & strcmpi(growthrates(:,15),'Continuous');
    growthdata = growthrates(idx,:);
    cond = {'aerobic','anaerobic'};
    for x = 1:length(cond)
        idx = ismember(lower(growthrates(:,1)),lower(strain))  & strcmpi(growthrates(:,15),'batch') & ~isnan(cell2mat(growthrates(:,3))) & ~isnan(cell2mat(growthrates(:,4))) & strcmpi(growthrates(:,14),cond{x});
        tmp = growthrates(idx,:);
        sub = unique(tmp(:,2));
        for j = 1:length(sub)
            idx = find(strcmpi(tmp(:,2),sub(j)));
            if length(idx) > 1
                [~,b] = max(cell2mat(tmp(idx,3)));
                median_tmp = tmp(idx(b),:);
            else
                median_tmp = tmp(idx,:);
            end
            if cell2mat(median_tmp(1,3)) > 0.1
                growthdata = [growthdata;median_tmp];
            end
        end
    end
   idx = ismember(lower(growthrates(:,1)),lower(strain));
    tmp = growthrates(idx,:);
    sub = unique(tmp(:,2));
    max_growth = {};
    for j = 1:length(sub)
        idx = find(strcmpi(tmp(:,2),sub(j)));
        [~,b] = max(cell2mat(tmp(idx,3)));
        if max(cell2mat(tmp(idx,3))) > 0.1
            max_growth = [max_growth;tmp(idx(b),:)];
        end
    end
    
    clearvars tmp mode_split Header Sequence
    
    disp('add enzyme data')
    % go to the EC file
    load('Protein_stoichiometry.mat') % generate in the function ECprepPanGEM.m should stored at pdbe % defines the subunit stoi for the complex downloaded from ComplexPortal, can be found in the Data/dabase/PDBe
    load('ecdata.mat') % generate in the function ECprepPanGEM.m should stored at ecGEMconstruction
    load('rxn2block.mat') % generate from the function blockbyproduct.m stored at ecGEMconstruction
    enzymedata = collectkcats(model,MWdata,ecdata,Protein_stoichiometry,false);
    save(['../../Results/model_build_files/model_classic/',strain,'_classic.mat'],'rxn2block','enzymedata','max_growth','growthrates','max_growth','model','MWdata','Protein_stoichiometry','strain','growthdata')
    
    % store enzymedata_classic
    enzymedata_classic = enzymedata;
    enzymedata_DL_org = collectPredictedKcat(model,MWdata,ecdata,Protein_stoichiometry,kcatpredictionPath,0); % 0 means that not with median
    
    % replace the DL predicted kcats with in vitro mesurement if available
    enzymedata = updatekcats(enzymedata_DL_org,enzymedata_classic);
    % enzymedata = matchkapp2kcats(enzymedata)
    save(['../../Results/model_build_files/model_dl/',strain,'_dl.mat'],'rxn2block','enzymedata','max_growth','growthrates','max_growth','model','MWdata','Protein_stoichiometry','strain','growthdata','enzymedata_DL_org')
    disp(['finish',strain])
end
