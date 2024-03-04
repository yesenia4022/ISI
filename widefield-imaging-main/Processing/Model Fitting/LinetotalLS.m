function [param ffit dom varaccount] = LinetotalLS(domain,f,G0)

%allows 'domain' vs. 'f' to be a scatter plot (i.e. domain is not
%monotonic)

global yy xx


%%%search%%%
xx = domain';
yy = f;
param = LinetotalLSfitter(G0);
%%%%%%%%%%%

dom = linspace(0,max(xx),100);

ffit = param(1)*dom + param(2);


%%%%
id = find(isnan(yy.*xx));
yy(id) = [];
xx(id) = [];
expect = param(1)*xx + param(2);
varaccount = (var(yy)-var(yy-expect))/var(yy);
%%%%
