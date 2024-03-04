function [param ffit varaccount] = NL(domain,f,G0)

%2 allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain;
yy = f;
param = NLfitter(G0);
%%%%%%%%%%%

domu = unique(xx);

ffit = param(1)*phi(domu-param(2)).^param(3) + param(4);
%ffit = param(1)*phi(domu-param(2)).^param(3);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = param(1)*phi(xx-param(2)).^param(3) + param(4);
%expect = param(1)*phi(xx-param(2)).^param(3);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
