
function G = gaussfitguess2Dpolar

%Double check Initial guesses

global RF;

r = length(RF(:,1));
c = length(RF(1,:));
ma = max(max(RF(:,ceil(c/2):end)));
[y x] = find(flipud(RF(:,ceil(c/2):end)) == ma);
y = y(1) - ceil(r/2);
x = x(1)-1;
r = sqrt(y.^2 + x.^2);
theta = atan(y/(x+eps))*180/pi;
if sign(theta) == -1
    theta = theta + 180;
end

G = [theta r 10 8 min(RF(:)) max(RF(:))-min(RF(:))];
