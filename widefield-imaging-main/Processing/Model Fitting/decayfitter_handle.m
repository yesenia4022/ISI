function err = decayfitter_handle(param)

global yy xx

A = param(1);
B = param(2);

ffit = 1./(A./xx + B);

err = mean((ffit-yy).^2);



