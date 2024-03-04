
function G = expfitguess

%Double check Initial guesses

global RF;

f = RF;

G = [RF(1)-RF(end) .5 RF(end)];
