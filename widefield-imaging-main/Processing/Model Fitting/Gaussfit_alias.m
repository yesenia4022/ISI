function [param ffit varacc sigma] = Gaussfit_alias(dom,f,circflag)

global RF dom

%%%%%%%%%%%
E = max(f);
f = f/E;
%%%%%%%%%%%%

if circflag == 1
    [ma ind] = max(f);
    f = circshift(f,[0 round(length(f)/2)-ind]);
end

%%%search%%%
RF = f;
param = gaussaliasfitter;

%%%%%%%%%%%

ddom = dom(2)-dom(1);



ffit = exp(-((1:length(f))-param(1)).^2/(2*param(2).^2));
shift360 = 360/ddom;
ffit = ffit + exp(-((1:length(f)) - param(1) - shift360).^2/(2*param(2).^2));
ffit = ffit + exp(-((1:length(f)) - param(1) + shift360).^2/(2*param(2).^2));
ffit = param(3)*ffit;

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

if circflag == 1
    ffit = circshift(ffit,[0 ind-round(length(ffit)/2)]);
    param(1) = param(1) -  (round(length(f)/2)-ind);
end

ffit = ffit*E;
param(3) = param(3)*E;


param(1) = (param(1)-1)*ddom;
if circflag
    if param(1) < 0
        param(1) = 180+param(1);
    end
end
param(2) = param(2)*ddom;

%% get sigma as 61% of max

xxI = linspace(0,180,1001);
xxI = xxI(1:end-1);

ffitI = exp((-xxI.^2)/(2*param(2)^2));

ffitI = ffitI-min(ffitI);
ffitI = ffitI/max(ffitI);

thresh = exp(-1/2);
thresh = 1/sqrt(2);
[dum id] = min(abs(ffitI-thresh));

sigma = (xxI(2)-xxI(1))*(id-1);



%%


