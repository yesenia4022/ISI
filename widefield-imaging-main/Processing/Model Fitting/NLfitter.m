function param = NLfitter(x0)

%options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
%param = fminsearch('expfitter2_handle',x0,options);
param = fminsearch('NLfitter_handle',x0);


