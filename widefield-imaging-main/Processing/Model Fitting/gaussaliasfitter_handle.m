function err = gaussaliasfitter_handle(param)

global RF dom;

dim = length(RF);

xc = param(1);
sx = param(2);

A = param(3);

xx = 1:dim;
xx = xx-xc;
img = A*exp(-xx.^2./(2*sx^2));

ddom = dom(2)-dom(1);
shift360 = 360/ddom;

xx = 1:dim;
xx = xx-xc-shift360;
img = img + A*exp(-xx.^2./(2*sx^2));

xx = 1:dim;
xx = xx-xc+shift360;
img = img + A*exp(-xx.^2./(2*sx^2));

err = sum((img-RF).^2);

% err = (img-RF).^2;
% err = sort(err);
% err = sum(err(1:round(length(err)*.9)));
