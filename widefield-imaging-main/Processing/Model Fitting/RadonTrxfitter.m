function x = RadonTrxfitter

global RF

x0 = RadonTrxguess;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('RadonTrxfitter_handl',x0,options);
[x f] = fminsearch('RadonTrx_handle',x0);

