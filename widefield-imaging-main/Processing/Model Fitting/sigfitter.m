function [x,f,g] = sigfitter

global RF xx initguess

if ~isempty(initguess)
    x0 = initguess;
else
    x0 = sigfitguess;
end
    
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('sigfitter_handle',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

xc = x(1);
sig = x(2);

A = x(3);
B = x(4);

d = xx-xc;

g = A./(1 + exp(-d*sig)) + B;
