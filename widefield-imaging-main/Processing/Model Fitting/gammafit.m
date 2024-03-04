function [param ffit varaccount] = gammafit(domain,f,G0)

%2 allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = gammafitter(G0);
%%%%%%%%%%%

domu = unique(xx);

ffit = domu.^param(1) + param(2);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
%expect = param(1)*phi(xx-param(2)).^param(3) + param(4);
expect = xx.^param(1) + param(2);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
