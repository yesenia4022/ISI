
function G = expfitguess3

%Double check Initial guesses

global yy base;

G = [max(yy)-base .1];
%G = [nanstd(yy) .1];
