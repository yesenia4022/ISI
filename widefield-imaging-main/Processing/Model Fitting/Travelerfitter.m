function [x,f,g] = Travelerfitter

global RF

x0 = Travelerfitguess;
dim = size(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);
%[x,f] = fminsearch('Travelerfitter_handle',x0);

% G(1) = ma_x;  %Spatial center
% G(2) = length(f(1,:))/5;  %spatial sigma
% G(3) = 1.5;   %slope (wave velocity)
% G(4) = ma_y;  %time to peak of center
% G(5) = 2;   %temporal sigma


LB = [1 .5 -1 1 2];
UB = [dim(2) dim(2)/2 3.5 dim(1) 10];
[x,f] = fmincon('Travelerfitter_handle',x0,[],[],[],[],LB,UB);



[xdom tdom] = meshgrid(1:length(f(1,:)),1:length(f(:,1)));

Xweight = exp(-(xdom-x(1)).^2/(2*x(2).^2));
%Xweight = exp(-abs(xdom-x(1))/x(2));

Delay = x(3)*abs(xdom-x(1)) + x(4);
Tfunc = exp(-(tdom-Delay).^2/(2*x(5).^2));

g = Xweight.*Tfunc;

%g = x(6)*g + x(7);

