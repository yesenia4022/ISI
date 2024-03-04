function err = CircGaussFit_handle(param)

global RF;

dim = length(RF);

xc = param(1);
kappa = param(2);

A = param(3);
base = param(4);

xx = linspace(0,2*pi,dim+1);
xx = xx(1:end-1);

xx = xx-xc;

img = exp(kappa*cos(xx));

img = img-min(img);
img = img/max(img);

img = A*img + base;

err = sum((img-RF).^2);