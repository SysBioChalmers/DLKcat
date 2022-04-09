%% pdbe_processing
%  Process data obtained from pdbe.

% The pdbe dataset was obtained 2020.04.26. modified from coFactorYeast without the cofactor extraction

% 1. Different pdb ids even with the same corresponding protein id could
% have different cofactor and protein stoichiometry information.
% Regarding cofactor, all types of cofactors would be selected for each
% protein from all its pdb ids, and the max copy number would be selected.
% Regarding protein stoichiometry, the pdb id with the highest coverage
% would be selected if a protein has multiple pdb ids. If this still leads
% to multiple stoichiometry result, then the highest stoichiometry would be
% selected.
% 2. Given that the collected pdbe information do not include all yeast
% genes, especially genes collected in the dataset but their homologue
% genes not. Therefore, yeast homologue genes information should be used to
% assign pdbe information for the genes that are not in the collected
% dataset but their homologues are.
% 3. The remaining genes are assumed to be monomer with no factors bound.

% Timing: ~ 1700 s

tic;

%% Load data

% load pdbe dataset
load('pdb.mat');%obtained 2020.04.26

% load yeast homologue genes
[~,txt,~] = xlsread('Yeast_homologue_genes.xlsx');%obtained 2020.01.02
list_gene = txt(2:end,2);
list_homologuegene = txt(2:end,5);
clear txt;

% load blastp result
[num,txt,~] = xlsread('blastp_res.xlsx');%obtained 2020.01.02
blastp = struct;
blastp.protein = txt(:,1);
blastp.pdb = txt(:,2);
blastp.identity = num(:,1);
blastp.coverage = num(:,2);
blastp.evalue = num(:,3);
blastp.bit = num(:,4);
clear num txt;

%% 1. filter gene id
% split multiple genes, remove 'NA' and non-yeast genes
pdb_processing_1 = struct();
pdb_processing_1.pdbid = cell(0,1);
pdb_processing_1.geneid = cell(0,1);
pdb_processing_1.name = cell(0,1);
pdb_processing_1.cofactor = cell(0,1);
pdb_processing_1.peptide_entityid = zeros(0,1);
pdb_processing_1.prot_stoic = zeros(0,1); % number of peptides per complex
pdb_processing_1.cofactor_stoic = cell(0,1); % number of cofactors in the complex
pdb_processing_1.pdb_chains = cell(0,1);
pdb_processing_1.all_pdb_chains = cell(0,1);

for i = 1:length(pdb.pdbid)
    disp(['Processing raw data: ' num2str(i) '/' num2str(length(pdb.pdbid))]);
    gene_tmp = pdb.geneid{i};
    gene_list_tmp = strsplit(gene_tmp,'; ')';
    gene_list_tmp = unique(gene_list_tmp);
    for j = 1:length(gene_list_tmp)
        gene_tmp_tmp = gene_list_tmp{j};
        if ~strcmp('NA',gene_tmp_tmp) && (strcmp('Q',gene_tmp_tmp(1))||strcmp('R',gene_tmp_tmp(1))||strcmp('Y',gene_tmp_tmp(1)))
            pdb_processing_1.pdbid = [pdb_processing_1.pdbid;pdb.pdbid(i)];
            pdb_processing_1.geneid = [pdb_processing_1.geneid;{gene_tmp_tmp}];
            pdb_processing_1.name = [pdb_processing_1.name;pdb.name(i)];
            pdb_processing_1.cofactor = [pdb_processing_1.cofactor;pdb.cofactor(i)];
            pdb_processing_1.peptide_entityid = [pdb_processing_1.peptide_entityid;pdb.peptide_entityid(i)];
            pdb_processing_1.prot_stoic = [pdb_processing_1.prot_stoic;pdb.prot_stoic(i)];
            pdb_processing_1.cofactor_stoic = [pdb_processing_1.cofactor_stoic;pdb.cofactor_stoic(i)];
            pdb_processing_1.pdb_chains = [pdb_processing_1.pdb_chains;pdb.pdb_chains(i)];
            pdb_processing_1.all_pdb_chains = [pdb_processing_1.all_pdb_chains;pdb.all_pdb_chains(i)];
        end
    end
end


%% 2. Collect protein stoichiometry information

genes_collected = unique(pdb_processing_1.geneid);
blastp_all_pdb = strcat(blastp.pdb,':',blastp.protein);

pdb_protein_stoichiometry = struct();
pdb_protein_stoichiometry.pdbid = cell(0,1);
pdb_protein_stoichiometry.protein = cell(0,1);
pdb_protein_stoichiometry.prot_stoic = zeros(0,1); % number of peptides per complex

for i = 1:length(genes_collected)
    disp(['Collecting protein stoichiometry: ' num2str(i) '/' num2str(length(genes_collected))]);
    gene_tmp = genes_collected(i);
    idx_tmp = ismember(pdb_processing_1.geneid,gene_tmp);
	loc_tmp = find(idx_tmp);

    cmp_num_tmp = zeros(0,4);
    cmp_txt_tmp = cell(0,1);
    for j = 1:length(loc_tmp)
        pdbid_tmp = pdb_processing_1.pdbid(loc_tmp(j));
        pdbchain_tmp = pdb_processing_1.all_pdb_chains(loc_tmp(j));
        pdbchain_tmp = split(pdbchain_tmp,'; ');
        for k = 1:length(pdbchain_tmp)
            str_tmp = strcat(pdbid_tmp,'_',pdbchain_tmp(k),':',gene_tmp);
            if ismember(str_tmp,blastp_all_pdb)
                idx_blastp_tmp = ismember(blastp_all_pdb,str_tmp);
                n_tmp = length(find(idx_blastp_tmp));
                identity_tmp = blastp.identity(idx_blastp_tmp);
                coverage_tmp = blastp.coverage(idx_blastp_tmp);
                evalue_tmp = blastp.evalue(idx_blastp_tmp);
                bit_tmp = blastp.bit(idx_blastp_tmp);
                row_tmp = [identity_tmp coverage_tmp evalue_tmp bit_tmp];
                cmp_num_tmp = [cmp_num_tmp;row_tmp];
                cmp_txt_tmp = [cmp_txt_tmp;repelem(str_tmp,n_tmp)'];
            else
                cmp_num_tmp = [cmp_num_tmp;[0 0 1000 0]];
                cmp_txt_tmp = [cmp_txt_tmp;str_tmp];
            end
        end
    end

    % select max coverage
    idx_cmp_tmp = cmp_num_tmp(:,2) == max(cmp_num_tmp(:,2));
    strlist_tmp = cmp_txt_tmp(idx_cmp_tmp);
    strlist_tmp = unique(cellfun(@(x) x(1:4),strlist_tmp,'UniformOutput',false));
    idx_tmp_tmp = ismember(pdb_processing_1.pdbid,strlist_tmp) & ismember(pdb_processing_1.geneid,gene_tmp);
    tmp_pdbid = pdb_processing_1.pdbid(idx_tmp_tmp);
    tmp_prot_stoic = pdb_processing_1.prot_stoic(idx_tmp_tmp);

    unique_prot_stoic = unique(tmp_prot_stoic);
    unique_pdbid = cell(length(unique_prot_stoic),1);
    for j = 1:length(unique_prot_stoic)
        tmp_pdbid_list = tmp_pdbid(tmp_prot_stoic == unique_prot_stoic(j));
        tmp_pdbid_list = unique(tmp_pdbid_list);
        unique_pdbid(j) = {strjoin(tmp_pdbid_list,'/')};
    end

    pdb_protein_stoichiometry.protein = [pdb_protein_stoichiometry.protein;repelem(gene_tmp,length(unique_prot_stoic))'];
    pdb_protein_stoichiometry.prot_stoic = [pdb_protein_stoichiometry.prot_stoic;unique_prot_stoic];
    pdb_protein_stoichiometry.pdbid = [pdb_protein_stoichiometry.pdbid;unique_pdbid];

end
save('raw_pdb_protein_stoichiometry.mat','pdb_protein_stoichiometry');

%% 4. Process the collected pdb data
Protein_stoichiometry = struct;
Protein_stoichiometry.protein = cell(0,1);
Protein_stoichiometry.stoichiometry = zeros(0,1);
Protein_stoichiometry.note = cell(0,1);
% For multiple pdb ids, the max stoichiometry would be chosen.
protein_list = unique(pdb_protein_stoichiometry.protein);
for i = 1:length(protein_list)
    proteinid = protein_list(i);
    idx_tmp = ismember(pdb_protein_stoichiometry.protein,proteinid);
    loc_tmp = find(idx_tmp);
    stoic_list = pdb_protein_stoichiometry.prot_stoic(idx_tmp);
    pdbid_list = pdb_protein_stoichiometry.pdbid(idx_tmp);
    max_stoic = max(stoic_list);
    max_pdbid = pdbid_list(stoic_list == max_stoic);
    pdbid = strcat('PDBid:',max_pdbid);
    Protein_stoichiometry.protein = [Protein_stoichiometry.protein;proteinid];
    Protein_stoichiometry.stoichiometry = [Protein_stoichiometry.stoichiometry;max_stoic];
    Protein_stoichiometry.note = [Protein_stoichiometry.note;pdbid];
end
% Add stoichiometry for homologue proteins
protein_list = list_gene(~ismember(list_gene,Protein_stoichiometry.protein));
for i = 1:length(protein_list)
    proteinid = protein_list(i);
    homol_protein = list_homologuegene(ismember(list_gene,proteinid));
    idx_tmp = ismember(homol_protein,Protein_stoichiometry.protein);
    if any(idx_tmp)
        list_homol_protein = homol_protein(idx_tmp);
        list_homol_stoic = Protein_stoichiometry.stoichiometry(ismember(Protein_stoichiometry.protein,list_homol_protein));
        max_stoic = max(list_homol_stoic); % choose max if multiple homologue proteins
        max_homol_protein = list_homol_protein(list_homol_stoic == max_stoic);
        homol_protein = strcat('Paralog:',max_homol_protein);
        Protein_stoichiometry.protein = [Protein_stoichiometry.protein;proteinid];
        Protein_stoichiometry.stoichiometry = [Protein_stoichiometry.stoichiometry;max_stoic];
        Protein_stoichiometry.note = [Protein_stoichiometry.note;homol_protein];
    end
end
save('Protein_stoichiometry.mat','Protein_stoichiometry');

toc;
