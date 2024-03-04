function err = gaussfitter_handle2Dpolar(param)

global RF;

dim = size(RF);

xc1 = param(1);
sx1 = param(2);
xc2 = param(3);
sx2 = param(4);

base = param(5);
A = param(6);

[xx yy] = meshgrid(1:dim(2),1:dim(1));
xx = xx-ceil(dim(2)/2);
yy = yy-ceil(dim(1)/2);
yy = flipud(yy);
tt = atan(yy./(xx+eps))*180/pi;
tt(:,1:floor(dim(2)/2)) = 180-fliplr(tt(:,ceil(dim(2)/2)+1:end));
rr = sqrt(yy.^2 + xx.^2);

tt = tt-xc1;
rr = rr-xc2;

img = A*exp(-rr.^2./(2*sx2^2)).*exp(-tt.^2./(2*sx1^2))  +  base;
img = img+fliplr(flipud(img));

err = sum((img(:)-RF(:)).^2);