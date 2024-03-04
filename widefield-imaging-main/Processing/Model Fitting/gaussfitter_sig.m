function [x,f,g] = gaussfitter_sig

global RF RF_sig

x0 = gaussfitguess;
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('gaussfitter_sig_handle',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

xc = x(1);
sig = x(2);

A = x(3);
B = x(4);

xx = 1:dim;

d = (xx-xc).^2;
g = A*exp(-d/(2*sig^2))+B;
