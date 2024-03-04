function err = spherefitter_handle(param)

global f;

x0 = param(1);
y0 = param(2);
z0 = param(3);
R = param(4);

dz = f(3,:)-param(3);
dy = f(2,:)-param(2);
dx = f(1,:)-param(1);
dR = sqrt(dz.^2 + dy.^2 + dx.^2); %Distance of data points from center of sphere


err = mean(abs(dR - R));
