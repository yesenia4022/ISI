function [param ffit varacc] = Gaussfit2D(domain,f,circflag)

global RF
global pepORI

%first dimension (rows) is circular variable

if circflag == 1
    ma = max(f(:));
    [indy indx] = find(f == ma);
    f = circshift(f,[round(length(f(:,1))/2)-indy 0]);
end

%%%search%%%
pepORI = domain;
RF = f;
param = gaussfitter2D;
%%%%%%%%%%%

ffitD1 = exp(-((1:length(f(:,1)))-param(1)).^2/(2*param(2).^2)); 
ffitD2 = exp(-((1:length(f(1,:)))-param(3)).^2/(2*param(4).^2));
ffit = param(5) + param(6)*(ffitD1' * ffitD2);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

if circflag == 1
    ffit = circshift(ffit,[indy-round(length(ffit(:,1))/2) 0]);
    param(1) = param(1) -  (round(length(f(:,1))/2)-indy);
end
 
