
function G = RadonTrxguess

%Double check Initial guesses

global RF pdomglob;

[dum,pori] = find(RF == max(RF(:)));

%[Sinwave amp; Sinwave phase; StdDev over position; Gain of funtion, baseline]

base = median(RF(:));

G = [.3*max(pdomglob) pi max(pdomglob)/10 max(RF(:))-base base];
