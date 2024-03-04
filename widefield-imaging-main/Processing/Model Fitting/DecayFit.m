function [param ffit varaccount] = DecayFit(domain,f,guess)

global yy xx initguess

initguess = guess;

%%%search%%%
xx = domain;
yy = f;
param = decayfitter;
%%%%%%%%%%%

A = param(1);
B = param(2);
ffit = 1./(A./xx + B);

%%%%
% id = find(isnan(yy.*xx));
% yy(id) = [];
% xx(id) = [];
% expect = param(1)./(1 + exp(-xx*param(2)));
varaccount = (var(yy)-var(yy-ffit))/var(yy);

%%%%
