function [param ffit varaccount] = Expfit3(domain,f,baseline)

%3 is an extension of 2, which only varies the decay constant and
%amplitude, while assuming it asymptotes to 'base'.

global yy xx base


%%%search%%%
xx = domain;
yy = f;
base = baseline;
param = expfitter3;

%%%%%%%%%%%


domu = unique(xx);

ffit = param(1)*exp(-param(2)*domu) + base;


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = param(1)*exp(-param(2)*xx) + base;
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%