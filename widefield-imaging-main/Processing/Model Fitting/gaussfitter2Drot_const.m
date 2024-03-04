function [x,f] = gaussfitter2Drot_const

global RF

x0 = gaussfitguess2Drot_const;

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-8,'TolX',.000004);  %this really helps for this function

% [ygp xgp sigguess sigguess origuess];
lb = [-2 -2 1 1 0];
ub = [2 2 1.5 1.5 180];

[x,f] = fmincon('gaussfitter_handle2Drot_const',x0,[],[],[],[],lb,ub,[],options);
%[x,f] = fminsearch('gaussfitter_handle2Drot',x0);
