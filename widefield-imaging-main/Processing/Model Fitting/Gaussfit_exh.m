function [param ffit varacc ffitI domI] = Gaussfit_exh(f,dom,mudom,sigdom,ampdom,basedom,circflag)

%%%%%%%%%%%
mi = min(f);
f = f-mi;
E = max(f);
f = f/E;
%%%%%%%%%%%%

if circflag == 1
    [ma ind] = max(f);
    f = circshift(f,[0 round(length(f)/2)-ind]);
end

e = zeros(length(mudom),length(sigdom),length(ampdom),length(basedom));
for m = 1:length(mudom)
    for s = 1:length(sigdom)
        for a = 1:length(ampdom)
            for b = 1:length(basedom)
                
                ffit = exp(-(dom-mudom(m)).^2/(2*sigdom(s).^2));
                ffit = ampdom(a)*ffit + basedom(b);
                
                e(m,s,a,b) = sum((ffit-f).^2);
                
            end
        end
    end
end

%Find the best parameters with the error block

dim = size(e);

[dum id] = min(e(:));
baseid = ceil(id/(dim(1)*dim(2)*dim(3)));
base = basedom(baseid);

id = id - (baseid-1)*(dim(1)*dim(2)*dim(3));
ampid = ceil(id/(dim(1)*dim(2)));
amp = ampdom(ampid);

id = id - (ampid-1)*(dim(1)*dim(2));
sigid = ceil(id/dim(1));
sig = sigdom(sigid);

id = id - (sigid-1)*dim(1);
muid = id;
mu = mudom(muid);

param = [mu sig amp base];


%%%%%%%%%%%%

ddom = dom(2)-dom(1);

ffit = exp(-(dom-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit + param(4);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

ffit = ffit*E + mi;
param(3) = param(3)*E;
param(4) = param(4)*E+mi;

if circflag == 1
    ffit = circshift(ffit,[0 ind-round(length(ffit)/2)]);
    param(1) = param(1) -  (round(length(f)/2)-ind)*ddom;
end

%%%
N = 100;
domI = linspace(0,length(f)*ddom,N+1);
domI = domI(1:end-1);

ffitI = exp(-(domI-param(1)).^2/(2*param(2).^2));
ffitI = param(3)*ffitI + param(4);
