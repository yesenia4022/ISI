function G = gaussfitguess

%Double check Initial guesses

global RF;

f = RF;
[M,idx] = max(f);
sig = 30/(180/length(f)); %approximately 30 deg
%sig = length(f)/8;

mi = prctile(f,.1);
ma = max(f);

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
G = [idx sig ma-mi mi];
%G = [length(f)/2 1 0 0];