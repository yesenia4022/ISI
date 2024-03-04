function [ro] = intersectfitter(ro)

global S D err iter
 
iter = 0;
%ro = intersectfitguess;

options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
[ro] = fminsearch('intersectfitter_handle',ro,options);
%[ro] = fminsearch('intersectfitter_handle',ro);

%figure,plot(err)
