function [param ffit varaccount domu] = Expfit2(domain,f,G0)

%2 allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = expfitter2(G0);
%%%%%%%%%%%

domu = unique(xx);

ffit = param(1)*exp(-domu*param(2)) + param(3);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = param(1)*exp(-xx*param(2)) + param(3);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
