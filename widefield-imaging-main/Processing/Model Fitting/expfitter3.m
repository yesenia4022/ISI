function param = expfitter3

global yy xx

x0 = expfitguess3;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('expfitter_handle',x0,options);
param = fminsearch('expfitter3_handle',x0);


