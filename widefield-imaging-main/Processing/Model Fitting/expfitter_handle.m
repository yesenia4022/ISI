function err = expfitter_handle(param)

global RF xx;

A = param(1);
alp = param(2);

B = param(3);

ffit = A*exp(-alp*xx) + B;

err = sum((ffit-RF).^2);

