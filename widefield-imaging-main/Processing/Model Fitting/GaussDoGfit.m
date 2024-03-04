function [param ffit varacc] = GaussDoGfit(f,circflag,ydom,xdom)

global RF sfdom oridom

%first dimension (rows) is circular variable

if circflag == 1
    ma = max(f(:));
    [indy indx] = find(f == ma);
    f = circshift(f,[round(length(f(:,1))/2)-indy 0]);
end

RF = f;
oridom = ydom;
sfdom = xdom;

%Make guess


[idx1 idx2] = find(RF == max(RF(:)));
x0 = [idx1 idx2 4 5 min(RF(:)) max(RF(:))-min(RF(:))];

%Search
% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% param = fminsearch('gaussfitter_handle',x0,options);
param = fminsearch('gaussDOGfitter_handle',x0);

%%%%%%%%%%%
Gfit = exp(-(ydom-param(1)).^2/(2*param(2).^2)); 
G1 = param(3)*exp(-xdom.^2/(2*param(4).^2));
G2 = exp(-xdom.^2/(2*param(5).^2));
DoG = G1 - G2;
ffit = param(7) + param(6)*(Gfit' * DoG);


varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

if circflag == 1
    ffit = circshift(ffit,[indy-round(length(ffit(:,1))/2) 0]);
    param(1) = param(1) -  (round(length(f(:,1))/2)-indy);
end
 
