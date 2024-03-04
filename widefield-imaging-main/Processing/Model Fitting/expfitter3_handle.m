function err = expfitter3_handle(param)

global yy xx base;

A = param(1);
alp = param(2);

ffit = A*exp(-alp*xx) + base;

err = nanmean((ffit-yy).^2);

