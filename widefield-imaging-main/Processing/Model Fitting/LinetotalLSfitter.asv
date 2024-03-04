function param = LinetotalLSfitter(x0)

global yy xx

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
%param = fminsearch('expfitter2_handle',x0,options);
param = fminsearch('LinetotalLSfitter_handle',x0);


