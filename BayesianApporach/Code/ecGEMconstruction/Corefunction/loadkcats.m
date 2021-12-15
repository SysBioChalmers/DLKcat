function allkcats = loadkcats

     KCAT_file      = 'max_KCAT.txt';

     %Extract BRENDA DATA from files information
     scallingFactor = 3600;   %scaling from [1/s] to [1/h]
     allkcats       = openDataFile(KCAT_file,scallingFactor); 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string_cells = stringSplit(cell_array)
         string_cells = {strsplit(cell_array,'//')};
         string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_cell = openDataFile(fileName,scallingFactor)
     fID          = fopen(fileName);
     data_cell    = textscan(fID,'%s %s %s %f %s','delimiter','\t','HeaderLines',1);
     fclose(fID);
     data_cell{4} = data_cell{4}*scallingFactor;
     %Split string for each organism in the BRENDA data 
     %{name, taxonomy, KEGG code}
     data_cell{3}  = cellfun(@stringSplit, data_cell{3});
end

