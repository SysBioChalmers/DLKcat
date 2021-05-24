function pc = pval_adjust(p, method)
% PVAL_ADJUST Adjust p-values for multiple comparisons. Given a set of
% p-values, returns p-values adjusted using one of several methods. This is
% an implementation of the p.adjust R function, the documentation of which
% can be found at http://www.inside-r.org/r-doc/stats/p.adjust.
%
%   pc = PVAL_ADJUST(p, method)
%
% Parameters:
%        p - Numeric vector or matrix of p-values. Contrary to the R
%            function, this function does not handle missing values.
%   method - Correction method, one of: 'holm', 'hochberg', 'hommel', 
%            'bonferroni', 'BH', 'BY', 'fdr', 'sidak' or 'none'.
%
% Outputs:
%       pc - A numeric vector or matrix of corrected p-values, with the
%            same shape/dimensions of p.
%
% Notes:
%
% This function has two main differences regarding the R implementation:
%
%   1. Contrary to the R function, this function does not handle missing
%      values.
%   2. It adds one additional correction method, 'sidak', as described
%      here: https://en.wikipedia.org/wiki/%C5%A0id%C3%A1k_correction
%
%
% Copyright (c) 2016 Nuno Fachada
% Distributed under the MIT License (See accompanying file LICENSE or copy 
% at http://opensource.org/licenses/MIT)
%

% Number of p-values
np = numel(p);

% Reshape input into a row vector, keeping original shape for later
% converting results into original shape
pdims = size(p);
p = reshape(p, 1, np);

% Method 'hommel' is equivalent to 'hochberg' of np == 2
if (np == 2) &&  strcmp(method, 'hommel')
    method = 'hochberg';
end;

% Just one p-value? Return it as given.
if np <= 1

    pc = p;
    
% What method to use?
elseif strcmp(method, 'holm')
    
    % Sort p-values from smallest to largest
    [pc, pidx] = sort(p);
    [~, ipidx] = sort(pidx);
    
    % Adjust p-values
    pc = min(1, cmax((np - (1:np) + 1) .* pc));
    
    % Put p-values back in original positions
    pc = pc(ipidx);

elseif strcmp(method, 'hochberg')

    % Descendent vector
    vdec = np:-1:1;
    
    % Sort p-values in descending order
    [pc, pidx] = sort(p, 'descend');
    
    % Get indexes of p-value indexes
    [~, ipidx] = sort(pidx);
    
    % Hochberg-specific transformation
    pc = ((np + 1) - vdec) .* pc;

    % Cumulative minimum
    pc = cmin(pc);
    
    % Reorder p-values to original order
    pc = pc(ipidx);

elseif strcmp(method, 'hommel')

    % Sort p-values from smallest to largest
    [pc, pidx] = sort(p);
    
    % Get indexes of p-value indexes
    [~, ipidx] = sort(pidx);
    
    % Generate vectors for cycle
    pa = repmat(min(np * pc ./ (1:np)), size(p));
    q = pa;
    
    % Begin cycle
    for i = (np - 1):-1:2
        i1 = 1:(np - i + 1);
        i2 = (np - i + 2):np;
        q1 = min(i * pc(i2) ./ (2:i));
        q(i1) = min(i * pc(i1), q1);
        q(i2) = q(np - i + 1);
        pa = max(pa, q);
    end;
    
    % Finalize result
    pa = max(pa, pc);
    pc = pa(ipidx);

elseif strcmp(method, 'bonferroni')
    
    % Simple conservative Bonferroni
    pc = p * numel(p);

elseif strcmp(method, 'BH') || strcmp(method, 'fdr')

    % Descendent vector
    vdec = np:-1:1;
    
    % Sort p-values in descending order
    [pc, pidx] = sort(p, 'descend');

    % Get indexes of p-value indexes
    [~, ipidx] = sort(pidx);

    % BH-specific transformation
    pc = (np ./ vdec) .* pc;

    % Cumulative minimum
    pc = cmin(pc);    
    
    % Reorder p-values to original order
    pc = pc(ipidx);

elseif strcmp(method, 'BY')

    % Descendent vector
    vdec = np:-1:1;
    
    % Sort p-values in descending order
    [pc, pidx] = sort(p, 'descend');

    % Get indexes of p-value indexes
    [~, ipidx] = sort(pidx);

    % BY-specific transformation
    q = sum(1 ./ (1:np));
    pc = (q * np ./ vdec) .* pc;

    % Cumulative minimum
    pc = cmin(pc);    
    
    % Reorder p-values to original order
    pc = pc(ipidx);

elseif strcmp(method, 'sidak')

    % Sidak correction
    pc = 1 - (1 - p) .^ np;
    
elseif strcmp(method, 'none')
    
    % No correction
    pc = p;
    
else
    
    % Unknown method
    error('Unknown p-value adjustment method');
    
end;

% Can't have p-values larger than one
pc(pc > 1) = 1;    

% Reshape result vector to original form
pc = reshape(pc, pdims);

% Helper function to determine the cumulative maximum
function p = cmax(p)

for i = 2:numel(p)
    if p(i) < p(i - 1)
        p(i) = p(i - 1);
    end;
end;

% Helper function to determine the cumulative minimum
function p = cmin(p)

for i = 2:numel(p)
    if p(i) > p(i - 1)
        p(i) = p(i - 1);
    end;
end;

