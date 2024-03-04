function [x,f,g] = gaussfitter2D

global RF

x0 = gaussfitguess2D;
dim = size(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);
[x,f] = fminsearch('gaussfitter_handle2D',x0);

xc1 = x(1);
sig1 = x(2);

xc2 = x(3);
sig2 = x(4);

B = x(5);
A = x(6);

[xx yy] = meshgrid(1:dim(2),1:dim(1));
d1 = (yy-xc1).^2;
d2 = (xx-xc2).^2;
g = A*exp(-d1/(2*sig1^2)).*exp(-d2/(2*sig2^2))+B;

%g = A/(sqrt(2*pi)*sig)*exp(-d/(2*sig^2))+B;
