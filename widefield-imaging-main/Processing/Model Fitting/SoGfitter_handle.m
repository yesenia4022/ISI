function err = SoGfitter_handle(param)

global RF dom;

xc = param(1);
sx = param(2);

A = param(3);
base = param(4);

xx = dom-xc;

f = A*exp(-xx.^2./(2*sx^2))  +  base;
f = f+fliplr(f);

err = sum((f-RF).^2);
