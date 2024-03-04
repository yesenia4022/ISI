function err = gammafitter_handle(param)

global yy xx;

gam = param(1);
base = param(2);


ffit = xx.^gam + base;

err = nanmean((ffit-yy).^2);


