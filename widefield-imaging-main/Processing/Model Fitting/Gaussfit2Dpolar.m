function [param ffit MSE] = Gaussfit2Dpolar(domain,f)

global RF
global pepORI

cols = length(f(1,:));
rows = length(f(:,1));
[xx yy] = meshgrid(1:cols,1:rows);
xx = xx-ceil(cols/2);
yy = yy-ceil(rows/2);
yy = flipud(yy);
tt = atan(yy./(xx+eps))*180/pi;
tt(:,1:floor(cols/2)) = 180-fliplr(tt(:,ceil(cols/2)+1:end));
rr = sqrt(yy.^2 + xx.^2);

%%%search%%%
pepORI = domain;
RF = f;
param = gaussfitter2Dpolar;
%%%%%%%%%%%

d1 = (tt-param(1)).^2;  %ori domain
d2 = (rr-param(3)).^2;  %sf domain
ffit = param(6)*exp(-d2./(2*param(4)^2)).*exp(-d1./(2*param(2)^2))  +  param(5);
ffit = ffit+fliplr(flipud(ffit));

MSE = mean((ffit(:)-f(:)).*(ffit(:)-f(:)));


