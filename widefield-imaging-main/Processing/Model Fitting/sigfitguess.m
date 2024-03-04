function G = sigfitguess

%Double check Initial guesses

global yy xx;


G = [median(xx) .2 (max(yy)-min(yy)) median(yy(:))];



% xc = param(1);
% sig = param(2);
% A = param(3);
% B = param(4);
% ffit = A./(1 + exp(-(xx-xc)*sig)) + B;