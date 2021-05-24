%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% table = truncateValues(table,cols)
%   table   cell array or table where some columns have values that should
%           be truncated
%   cols    index or indices of columns with values to be truncated
%
% Benjamin Sanchez    Last edited: 2018-05-25
% Eduard Kerkhoven    Last edited: 2018-12-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table = truncateValues(table,cols)

[m,n] = size(table);
for i = 1:m
    for j = cols
        orderMagn  = max([ceil(log10(abs(table{i,j}))),0]);
        table{i,j} = round(table{i,j},6-orderMagn);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
