function [param ffit varaccount] = Expfit2Drot(f,varargin)

global RF errall
errall = [];

f = f-nanmedian(f(:));
f = f/max(f(:));

if ~isempty(varargin{1})
    global origuess 
    origuess = varargin{1};
end

%%%search%%%
RF = f;
param = gaussfitter2Drot;
%%%%%%%%%%%

imW = length(f(:,1));
domx = -floor(imW/2):ceil(imW/2)-1;
domy = ceil(imW/2)-1:-1:-floor(imW/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(5)*pi/180) + y*sin(param(5)*pi/180);
yp = y*cos(param(5)*pi/180) - x*sin(param(5)*pi/180); 

ffitD1 = exp(-abs(yp-param(1))*param(3));
ffitD2 = exp(-abs(yp-param(2))*param(4));

ffit = ffitD1.*ffitD2;

varaccount = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%Put the means back on the original axis

p2 = param(2)*cos(-param(5)*pi/180) + param(1)*sin(-param(5)*pi/180);
p1 = param(1)*cos(-param(5)*pi/180) - param(2)*sin(-param(5)*pi/180);
param(1) = p1;
param(2) = p2;

%figure,plot(errall)