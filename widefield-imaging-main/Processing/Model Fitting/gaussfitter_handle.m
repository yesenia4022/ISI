function err = gaussaliasfitter_handle(param)

global RF;

dim = length(RF);

xc = param(1);
sx = param(2);

A = param(3);
base = param(4);

xx = 1:dim;

xx = xx-xc;

img = A*exp(-xx.^2./(2*sx^2))  +  base;

err = sum((img-RF).^2);

% err = (img-RF).^2;
% err = sort(err);
% err = sum(err(1:round(length(err)*.9)));
