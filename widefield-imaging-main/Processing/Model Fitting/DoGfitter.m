function x = DoGfitter

x0 = DoGfitguess;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);
[x,f] = fminsearch('DoGfitter_handle',x0);

