
function G = expfitguess2

%Double check Initial guesses

global yy;

id = find(~isnan(yy));

%G = [prctile(yy(id),5)-prctile(yy(id),90) .1 prctile(yy(id),15)];
G = [max(yy)-min(yy) .1 min(yy)];
%G = [-(max(yy)-min(yy)) .5 max(yy)];
