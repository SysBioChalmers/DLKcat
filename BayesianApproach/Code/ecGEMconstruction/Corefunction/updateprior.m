function [a,b] = updateprior(x)
x = x{:};

%  try
%     [~, pd, phat_ci]  = fitdist_ntrunc(log10(x./3600),[0,8]); % fix the normal distribution stuff with trucate
% 
%     if pd(1) < 0
%         pd(1) = 0;
%         pd(2) = 1;
%        
%         warning('reach limit')
%     elseif pd(1) > 8
%         pd(1) = 8;
%         pd(2) = 1;
%     end
% 
%     a = 10^(pd(1))*3600;
%     b= 10^(pd(2))*3600;
%  
% catch
    pd = fitdist(log10(x./3600),'Normal');
%     if pd.mu < 0
%        pd.mu = 0;
%         pd.sigma = 1;
%        
%         warning('reach limit')
%     elseif pd.mu > 8
%         pd.mu = 8;
%         pd.sigma = 1;
%     end
    a = 10^(pd.mu)*3600;
   % b= 10^(pd.sigma)*3600;
   b = pd.sigma;
%    warning('turncate not available')
    
end