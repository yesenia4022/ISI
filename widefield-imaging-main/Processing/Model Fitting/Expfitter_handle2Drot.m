function err = Expfitter_handle2Drot(param)

global RF cent origuess;

imW = length(RF(:,1));
domx = -floor(imW/2):ceil(imW/2)-1;
domy = ceil(imW/2)-1:-1:-floor(imW/2);

xg = domx(cent(2));
yg = domy(cent(1));

[x y] = meshgrid(domx,domy);

xgp = xg*cos(origuess*pi/180) + yg*sin(origuess*pi/180);
ygp = yg*cos(origuess*pi/180) - xg*sin(origuess*pi/180);

xp = x*cos(param(5)*pi/180) + y*sin(param(5)*pi/180);
yp = y*cos(param(5)*pi/180) - x*sin(param(5)*pi/180);

ffitD1 = exp(-abs(yp-ygp)*param(1));
ffitD2 = exp(-abs(xp-xgp)*param(2));
ffit = param(3)*ffitD1.*ffitD2 + param(4);

err = nanmean((ffit(:)-RF(:)).*(ffit(:)-RF(:)));

% id = find(~isnan(ffit(:).*RF(:)));
% ffit2 = ffit(id);
% RF2 = RF(id);
% 
% err = corrcoef(ffit2(:),RF2(:));
% err = -err(1,2);

