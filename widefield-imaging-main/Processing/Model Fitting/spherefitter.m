function [p er] = spherefitter

p0 = spherefitguess;

 %options = optimset('MaxFunEvals',600000,'MaxIter',600000,'TolFun',.000004,'TolX',.000004);
 %[p,er] = fminsearch('spherefitter_handle',p0,options);
[p,er] = fminsearch('spherefitter_handle',p0);

