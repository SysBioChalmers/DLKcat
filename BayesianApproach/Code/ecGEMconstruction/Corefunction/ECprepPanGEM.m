function [ecdata,Protein_stoichiometry] = ECprepPanGEM
load('panmodel.mat')  %MLKcat/BayesianApproach/Results/ssGEMs
model = model_original;
% Reformulate the orginal model
%   1. Reactions with isozymes (i.e., 'OR' case in GPR) will be copied,
%   each isozyme will be assigned to each copy. 
%   2. Reversible reactions will be split into forward and reverse ones.

model_split = splitModel(model);
[grRules,rxnGeneMat]   = standardizeGrRules(model_split,true);

org_name = 'saccharomyces cerevisiae';

% generate a ec file for all pangenes %%%%%%important
%% Collect EC number
org_name_tmp = strrep(org_name,' ','+');
address = ['https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:',org_name_tmp,'&columns=genes(OLN),ec&format=tab'];
websave('uniprot.tab',address);
%Build Swissprot table (uniprot code - protein name - gene names - EC number - sequence):
fileID_uni        = fopen('uniprot.tab');
swissprot         = textscan(fileID_uni,'%s %s','delimiter','\t');
swissprot         = [swissprot{1} swissprot{2}];
swissprot(1,:)    = [];
fclose(fileID_uni);
tmp = cellfun(@isempty,swissprot); % filter out empty cell
swissprot = swissprot(~sum(tmp,2),:);
id_Uniprot = swissprot(:,1);
ec_Uniprot = swissprot(:,2);
ec_Uniprot = cellfun(@(x) ['EC',x], ec_Uniprot, 'UniformOutput',false);
ec_Uniprot = cellfun(@(x) strrep(x,'; ','; EC'), ec_Uniprot, 'UniformOutput',false);
%Download KEGG data:
websave('keggorg.txt','http://rest.kegg.jp/list/organism');
fID         = fopen('keggorg.txt');
Data    = textscan(fID,'%s%s%s%s','Delimiter','\t','HeaderLines',1);
fclose(fID);
mapping(:,1) = Data{2};
mapping(:,2) = Data{3};
tmp = regexprep(mapping(:,2),'\([^)]*\)','');
tmp = strtrim(tmp);
idx = find(ismember(lower(tmp),lower(org_name)));
if idx ~=0
    keggID = mapping{idx,1};
else
    error('check org_name to see whether it is in the kegg list:http://rest.kegg.jp/list/organism ')
end
clear Data mapping tmp;

base = 'http://rest.kegg.jp/link/ec/';
websave('ec_number_kegg.txt',[base keggID]);
fID         = fopen('ec_number_kegg.txt');
Data    = textscan(fID,'%s%s','Delimiter','\t','HeaderLines',1);
fclose(fID);
mapping(:,1) = Data{1};
mapping(:,2) = Data{2};
mapping(:,1) = strrep(mapping(:,1),[keggID,':'],'');
mapping(:,2) = strrep(mapping(:,2),'ec:','EC');

id_KEGG = mapping(:,1);
ec_KEGG = mapping(:,2);

id_list = unique([id_Uniprot;id_KEGG]);

ecdata = struct();
ecdata.id = cell(0,1);
ecdata.ec = cell(0,1);

for i = 1:length(id_list)
    id = id_list(i);
    if ismember(id,id_Uniprot)
        ec_U_tmp = ec_Uniprot(ismember(id_Uniprot,id));
        ec_U_tmp = split(ec_U_tmp);
        id_U_tmp = repelem(id,length(ec_U_tmp))';
        ecdata.id = [ecdata.id;id_U_tmp];
        ecdata.ec = [ecdata.ec;ec_U_tmp];
    end
    if ismember(id,id_KEGG)
        ec_K_tmp = ec_KEGG(ismember(id_KEGG,id));
        id_K_tmp = repelem(id,length(ec_K_tmp))';
        ecdata.id = [ecdata.id;id_K_tmp];
        ecdata.ec = [ecdata.ec;ec_K_tmp];
    end
end
ecdata.ec = strrep(ecdata.ec,';','');
%% collect for each reaction with GPR
% update the subunit collection based on 1201 orthologs, the assumption is
% that orthologs should have the same copies in the complex formation
load('Protein_stoichiometry.mat');% obtained from pdbe for stoich for subunit of yeast ComplementaryData/databases/pcdb/
load('ortholog_changedtoPanID.mat') % all orthologs in yeast8 ComplementaryData/

[~,idx] = ismember(ortholog(:,1),Protein_stoichiometry.protein);
orth_stoi = Protein_stoichiometry.stoichiometry(idx(idx~=0));
orth = ortholog(idx~=0,2);
Protein_stoichiometry.protein = [Protein_stoichiometry.protein;orth];
Protein_stoichiometry.stoichiometry = [Protein_stoichiometry.stoichiometry;orth_stoi];

[~,idx] = ismember(ortholog(:,1),ecdata.id);
orth_stoi = ecdata.ec(idx(idx~=0));
orth = ortholog(idx~=0,2);
ecdata.id = [ecdata.id;orth];
ecdata.ec = [ecdata.ec;orth_stoi];

%% get module info for each enzyme
ecdata.class(1:length(ecdata.ec),1) = {''};
fileName = 'module_ec.txt';
fid  = fopen(fileName);
data = textscan(fid,'%s %s %s','Delimiter','\t');
datanew(:,1) = data{1};
datanew(:,2) = data{2};
datanew(:,3) = data{3};
clear data;
fclose(fid);
datanew = strrep(datanew,'ec:','EC');
[~,idx_sub] = ismember(ecdata.ec,datanew(:,2));
ecdata.class(idx_sub~=0,1) = datanew(idx_sub(idx_sub~=0),3);

%% get kcat range for each ec number
ecdata.kcatrange(1:length(ecdata.ec),1:2) = 0;
load('kcatrange.mat')
ecdatatmp = strrep(ecdata.ec,'EC','');
[~,idx_ec] = ismember(ecdatatmp,kcatrange(:,1));
ecdata.kcatrange(idx_ec~=0,:) = cell2mat(kcatrange(idx_ec(idx_ec~=0),2:3));
