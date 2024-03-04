function [param] = Spherefit(S)

global f

f = S;
%%%search%%%
[param er] = spherefitter
%%%%%%%%%%%

% dz = f(3,:)-param(3);
% dy = f(2,:)-param(2);
% dx = f(1,:)-param(1);
% dR = sqrt(dz.^2 + dy.^2 + dx.^2); %Distance of data points from center of sphere
% xyproj = sqrt(dx.^2 + dy.^2);  %xy projection
% phi = acos(xyproj./dR);
% theta = acos(dx./xyproj);
% 
% x = param(1) + param(4)*cos(phi)*cos(theta);
% y = param(2) + param(4)*cos(phi)*sin(theta);
% z = param(3) + param(4)*sin(phi);

%MSE = mean((ffit-f).*(ffit-f));

