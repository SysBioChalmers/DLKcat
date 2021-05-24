
% extract sabiork data
     KCAT_file      = 'Kcat_sabio.txt';

     %Extract sabiork DATA from files information
     scallingFactor = 1;   %keep 1/s
     data_cell = openDataFile_SB(KCAT_file,scallingFactor);
     data_cell = [data_cell{1},lower(data_cell{2}),lower(data_cell{3}),num2cell(data_cell{4}),data_cell{5},data_cell{6}];
     idx = ismember(data_cell(:,5),'wildtype');
     data_cell = data_cell(idx,[1:4,6]);
     idx = contains(data_cell(:,2),';');
     new_data_cell = data_cell(~idx,:);
     idx = find(idx);
     for j = 1:length(idx)
         j
         sub = split(data_cell(idx(j),2),';');
         tmp(:,1) = repmat(data_cell(idx(j),1),length(sub),1);
         tmp(:,2) = sub;
         tmp(:,3) = repmat(data_cell(idx(j),3),length(sub),1);
         tmp(:,4) = repmat(data_cell(idx(j),4),length(sub),1);
         tmp(:,5) = repmat(data_cell(idx(j),5),length(sub),1);
         new_data_cell = [new_data_cell;tmp];
         clear tmp;
     end
     


     KCAT_file      = 'max_KCAT.txt';

     %Extract BRENDA DATA from files information
     scallingFactor = 1;   %[1/s] -> [1/h]
     allkcats       = openDataFile(KCAT_file,scallingFactor); 
     allkcats = [allkcats{1},allkcats{2},allkcats{3},num2cell(allkcats{4}),allkcats{5}];
depreceted   =[];  
     % combine BRENDA and Sabiork
     for i = 1:length(new_data_cell(:,1))
         i
         if any(ismember(allkcats(:,1),new_data_cell(i,1)) & ismember(allkcats(:,2),new_data_cell(i,2)) & ismember(allkcats(:,3),new_data_cell(i,3))) 
             idx = find(ismember(allkcats(:,1),new_data_cell(i,1)) & ismember(allkcats(:,2),new_data_cell(i,2)) & ismember(allkcats(:,3),new_data_cell(i,3)));
             if length(idx) > 1
                 error('multi mapping')
             else
                [~,I] = max(cell2mat([allkcats(idx,4),new_data_cell(i,4)]));
                if I == 2
                    allkcats(idx,:) = new_data_cell(i,:);
                    depreceted = [depreceted;allkcats(idx,:);new_data_cell(i,:)];
                end
             end
         else
             allkcats = [allkcats;new_data_cell(i,:)];
         end
     end
     writetable(cell2table(allkcats),'Max_combined.txt','Delimiter','\t')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string_cells = stringSplit(cell_array)
         string_cells = {strsplit(cell_array,'//')};
         string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_cell = openDataFile(fileName,scallingFactor)
     fID          = fopen(fileName);
     data_cell    = textscan(fID,'%s %s %s %f %s','delimiter','\t');
     fclose(fID);
     data_cell{4} = data_cell{4}*scallingFactor;
     %Split string for each organism in the BRENDA data 
     %{name, taxonomy, KEGG code}
     data_cell{3}  = cellfun(@stringSplit, data_cell{3});
     data_cell{5} = repmat({'brenda'},length(data_cell{1}),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_cell = openDataFile_SB(fileName,scallingFactor)
     fID          = fopen(fileName);
     %EntryID	Type	ECNumber	Substrate	EnzymeType	PubMedID	Organism	UniprotID	Value	Unit
     data_cell    = textscan(fID,'%s %s %s %s %s %s %s %s %f %s','delimiter','\t','Headerlines',1);
     fclose(fID);
     data_cell = [data_cell(3),data_cell(4),data_cell(7),data_cell(9),data_cell(5)];
     data_cell{4} = data_cell{4}*scallingFactor;
     %add EC in the EC number 
     data_cell{1}  = cellfun(@(x) ['EC',x], data_cell{1}, 'UniformOutput',false);
     data_cell{6} = repmat({'sabiork'},length(data_cell{1}),1);
     data_cell{2}  = strrep(data_cell{2},'H2O2','%hyperoxide%');% remove H2O from substrate
     data_cell{2}  = strrep(data_cell{2},'H2O;','');% remove H2O from substrate
     data_cell{2}  = strrep(data_cell{2},';H2O','');% remove H2O from substrate
     data_cell{2}  = strrep(data_cell{2},'%hyperoxide%','H2O2');
     
     % split the substrate
     
end