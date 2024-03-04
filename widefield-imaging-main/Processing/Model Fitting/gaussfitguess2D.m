
function G = gaussfitguess2D

%Double check Initial guesses

global RF;

[idx1 idx2] = find(RF == max(RF(:)));

G = [idx1 idx2 4 5 min(RF(:)) max(RF(:))-min(RF(:))];
