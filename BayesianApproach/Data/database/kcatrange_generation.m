function kcatrange_generation
fileName = 'Kcat_combination_41559.tsv';
fid  = fopen(fileName);
data = textscan(fid,[repmat('%s ',[1,7]) '%s'],'Delimiter','\t','Headerlines',1);
data = [data{1:end}];
fclose(fid);


% extract all wildtype
idx = strcmp(data(:,2),'wildtype');
data = data(idx,:);
tmp = tabulate(data(:,1));
ecnumber = tmp(cell2mat(tmp(:,2))>10,:);
kcatrange = zeros(length(ecnumber),2);
for i = 1:length(ecnumber)
    idx = strcmp(data(:,1),ecnumber(i));
    kcat_tmp = cellfun(@str2num,data(idx,7));
    kcatrange(i,1) = min(log10(kcat_tmp));
    kcatrange(i,2) = max(log10(kcat_tmp));
end
kcatrange = [ecnumber(:,1),num2cell(kcatrange)];
% save file
save('kcatrange.mat','kcatrange')
    
end