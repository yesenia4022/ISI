function err = sigfitter_handle(param)

global yy xx

xc = param(1);
sig = param(2);

A = param(3);
B = param(4);

d = xx-xc;

ffit = A./(1 + exp(-d*sig)) + B;

err = mean((ffit-yy).^2);




