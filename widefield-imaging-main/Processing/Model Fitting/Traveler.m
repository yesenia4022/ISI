function [param ffit varaccount] = Traveler(f)

global RF

%%%%%%%%%%%
mi = min(f(:));
f = f-mi;
E = max(f(:));
f = f/E;
%%%%%%%%%%%%

%%%search%%%
RF = f;
param = Travelerfitter;
%%%%%%%%%%%

[xdom tdom] = meshgrid(1:length(f(1,:)),1:length(f(:,1)));

Xweight = exp(-(xdom-param(1)).^2/(2*param(2).^2));
%Xweight = exp(-abs(xdom-param(1))/param(2));

Delay = param(3)*abs(xdom-param(1)) + param(4);
Tfunc = exp(-(tdom-Delay).^2/(2*param(5).^2));

ffit = Xweight.*Tfunc;
%ffit = param(6)*ffit + param(7);

varaccount = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

mudum = param(1);

