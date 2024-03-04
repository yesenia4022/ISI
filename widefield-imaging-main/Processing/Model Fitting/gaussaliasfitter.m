function [x,f,g] = gaussaliasfitter

global RF dom

x0 = gaussaliasfitguess;
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('gaussaliasfitter_handle',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

xc = x(1);
sig = x(2);

A = x(3);

xx = 1:dim;
d = (xx-xc).^2;
g = A*exp(-d/(2*sig^2));

ddom = dom(2)-dom(1);
shift360 = 360/ddom;

d = (xx-xc-shift360).^2;
g = g + A*exp(-d/(2*sig^2));

d = (xx-xc+shift360).^2;
g = g + A*exp(-d/(2*sig^2));
