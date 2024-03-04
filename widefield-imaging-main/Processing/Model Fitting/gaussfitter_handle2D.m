function err = gaussfitter_handle2D(param)

global RF;

dim = size(RF);

xc1 = param(1);
sx1 = param(2);
xc2 = param(3);
sx2 = param(4);

base = param(5);
A = param(6);

[xx yy] = meshgrid(1:dim(2),1:dim(1));
yy = yy-xc1;
xx = xx-xc2;

img = A*exp(-xx.^2./(2*sx2^2)).*exp(-yy.^2./(2*sx1^2))  +  base;

err = sum((img(:)-RF(:)).^2);