function x = SoGfitter

x0 = SoGfitguess;

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('SoGfitter_handle',x0,options);

[x] = fminsearch('SoGfitter_handle',x0);
