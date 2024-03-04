function [x,f] = gaussfitter2Drot

global RF

x0 = gaussfitguess2Drot;

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',1e-8,'TolX',.000004);  %this really helps for this function
[x,f] = fminsearch('gaussfitter_handle2Drot',x0,options);
%[x,f] = fminsearch('gaussfitter_handle2Drot',x0);
