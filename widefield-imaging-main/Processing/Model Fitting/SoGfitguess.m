function G = SoGfitguess

%Double check Initial guesses

global RF dom;

f = RF;
[M,idx] = max(f((length(f)/2+1):end));
idx = idx+length(f)/2;
sig = 1; 

mi = prctile(f,.1);
ma = max(f);

%G = [idx length(f)/6 max(RF)-min(RF) min(RF)];
G = [dom(idx) sig ma-mi mi];
%G = [length(f)/2 1 0 0];