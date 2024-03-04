function err = NLfitter_handle(param)

global yy xx;

k = param(1);
T = param(2);
n = param(3);
B = param(4);


ffit = k*phi(xx-T).^n + B;
%ffit = k*phi(xx-T).^n;

%id = find(~isnan(yy.*ffit));
%err = trimmean((ffit(id)-yy(id)).^2,20);

err = nanmean((ffit-yy).^2);


