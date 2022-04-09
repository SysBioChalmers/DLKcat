function r = getrSample(mu,sigma,lb,ub,step,method)
if nargin < 6
    method = 'normal';
end
if lb == ub && lb == 0
    lb = -2;
    ub = 9;
end
if strcmp(method,'normal')
    mutmp = log10(mu/3600);
    %sigmatmp = log10(sigma/3600);
    sigmatmp = sigma;
    pd = makedist('normal','mu',mutmp,'sigma',sigmatmp);
    %t = truncate(pd,-3,8);
    r = random(pd,1,step);
    r(r < lb) = lb;
    r(r > ub) = ub;
    r = 10.^(r).*3600;
elseif strcmp(method,'uniform')
    mutmp = log10(mu/3600);
    sigmatmp = sigma;
    pd = makedist('uniform','lower',mutmp-sigmatmp,'upper',mutmp + sigmatmp);
    t = truncate(pd,-2,8);
    r = random(t,1,step);
    r = 10.^(r).*3600;
end
end