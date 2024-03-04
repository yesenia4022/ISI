function param = saturationfitter(x0)

global yy xx

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
%param = fminsearch('expfitter2_handle',x0,options);
param = fminsearch('saturationfitter_handle',x0);


