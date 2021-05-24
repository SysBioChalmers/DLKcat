function gene = convert2panID(gene)


fileName    = ['../../../Multi_scale_evolution/pan_genome/result/id_mapping/Saccharomyces_cerevisiae.tsv'];
fID         = fopen(fileName);
protData    = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
fclose(fID);
geneID_core = protData{2}; % geneID in 343 yeast species with @seq
panID_final = protData{5}; % panID
draftgeneID      = protData{6}; % mRNAID geneID CDS
draftgeneID      = extractBefore(draftgeneID,' ');
panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');

[~,idx] = ismember(gene,draftgeneID);
gene(idx ~= 0) = panID_final(idx(idx ~= 0));
end