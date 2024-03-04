function [x,f] = Expfitter2Drot

global RF

x0 = Expfitguess2Drot;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle2Drot',x0,options);
[x,f] = fminsearch('Expfitter_handle2Drot',x0);
