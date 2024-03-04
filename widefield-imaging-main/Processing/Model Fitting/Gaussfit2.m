function [param ffit mudum] = Gaussfit2(domain,f,circflag)

global RF
global pepORI

%%%%%%%%%%%
mi = min(f);
f = f-mi;
E = norm(f);
f = f/E;
%%%%%%%%%%%%

if circflag == 1
    [ma ind] = max(f);
    f = circshift(f,[0 round(length(f)/2)-ind]);
end

%%%search%%%
pepORI = domain;
RF = f;
param = gaussfitter;
%%%%%%%%%%%

ffit = exp(-((1:length(f))-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit + param(4);

MSE = mean((ffit-f).*(ffit-f));

mudum = param(1);

if circflag == 1
    ffit = circshift(ffit,[0 ind-round(length(ffit)/2)]);
    param(1) = param(1) -  (round(length(f)/2)-ind);
end

ffit = ffit*E + mi;
param(3) = param(3)*E;
param(4) = param(4)*E+mi;
