function [x,f,g] = expfitter

global RF x0 xx

%dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('expfitter_handle',x0,options);
[x,f] = fminsearch('expfitter_handle',x0);

A = x(1);
alp = x(2);

B = x(3);

g = A*exp(-alp*xx)+B;
