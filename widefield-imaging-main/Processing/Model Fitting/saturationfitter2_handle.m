function err = saturationfitter2_handle(param)

global yy xx;

A = param(1);
alp = param(2);

B = param(3);

ffit = A*(1-exp(-xx*alp)) + B;

%id = find(~isnan(yy.*ffit));
%err = trimmean((ffit(id)-yy(id)).^2,20);

err = nanmean((ffit-yy).^2);


