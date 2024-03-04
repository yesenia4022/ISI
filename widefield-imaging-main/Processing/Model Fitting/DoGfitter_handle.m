function err = DoGfitter_handle(param)

global RF logdom;

dim = length(RF);

sx1 = param(1);
A1 = param(2);

sx2 = param(3);
A2 = param(4);

base = param(5);

img1 = A1*exp(-logdom.^2./(2*sx1^2));
img2 = A2*exp(-logdom.^2./(2*sx2^2));

img = img1-img2+base;

err = sum((img-RF).^2);
