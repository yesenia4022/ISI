function [mu sigma ffit varacc param] = CircGaussFit(f)

global RF

%%%%%%%%%%%
% mi = min(f);
% f = f-mi;
% E = max(f);
% f = f/E;
%%%%%%%%%%%%

RF = f;
 
%% make initial guess

xx = linspace(0,2*pi,length(f)+1);
xx = xx(1:end-1);

[dum idx] = max(f);

mi = prctile(f,.1);
ma = max(f);

x0 = [xx(idx) 2 ma-mi mi];

%% search
param = fminsearch('CircGaussFit_handle',x0);

%% Make a highly sampled fit to produce the sigma
xxI = linspace(0,2*pi,1001);
xxI = xxI(1:end-1);

ffitI = exp(param(2)*cos(xxI-param(1)));

[dum idma] = max(ffitI);
ffitI = circshift(ffitI,[0 1-idma]);  %shift peak to first element
ffitI = ffitI-min(ffitI);
ffitI = ffitI/max(ffitI);

thresh = exp(-1/2);
[dum id] = min(abs(ffitI(1:length(ffitI)/2)-thresh));

sigma = (xxI(2)-xxI(1))*(id-1);
sigma = sigma/2*180/pi;  %put it in the orientation domain
mu = param(1)/2*180/pi;

%% Make fit with same sample rate as the original

ffit = exp(param(2)*cos(xx-param(1)));

ffit = ffit-min(ffit);
ffit = ffit/max(ffit);

ffit = param(3)*ffit + param(4);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));


