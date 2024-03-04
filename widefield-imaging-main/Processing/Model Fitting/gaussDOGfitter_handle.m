function err = gaussDOGfitter_handle(param)

global RF oridom sfdom;

Gfit = exp(-(oridom-param(1)).^2/(2*param(2).^2)); 
G1 = param(3)*exp(-sfdom.^2/(2*param(4).^2));
G2 = exp(-sfdom.^2/(2*param(5).^2));
DoG = G1 - G2;
ffit = param(7) + param(6)*(Gfit' * DoG);

err = sum((ffit(:)-RF(:)).^2);