function err = saturationfitter_handle(param)

global yy xx errvec;

A = param(1);
B = param(2);

lindom = linspace(0,max(xx),200);
ffit1 = A*(1 - exp(-lindom*B));

temp = ones(1,length(ffit1));
E = zeros(1,length(yy));
for i = 1:length(yy)
    Ds = sqrt((yy(i)*temp-ffit1).^2 + (xx(i)*temp-lindom).^2);
    E(i) = min(Ds);
end
id = find(~isnan(E));
err = trimmean(E(id),10);

%errvec = [errvec err];

%err = nanmean(abs(ffit1-yy));



