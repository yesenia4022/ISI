function [param ffit varaccount] = Gaussfit2Drot_const(f,varargin)

global RF errall
errall = [];

f = f-nanmedian(f(:)) + .1;
f = f/max(f(:));

if ~isempty(varargin{1})
    global origuess 
    origuess = varargin{1};
end

%%%search%%%
RF = f;
param = gaussfitter2Drot_const;
%%%%%%%%%%%

imWy = length(f(:,1));
imWx = length(f(1,:));
domx = -floor(imWx/2):ceil(imWx/2)-1;
domy = ceil(imWy/2)-1:-1:-floor(imWy/2);

[x y] = meshgrid(domx,domy);

xp = x*cos(param(5)*pi/180) + y*sin(param(5)*pi/180);
yp = y*cos(param(5)*pi/180) - x*sin(param(5)*pi/180); 

ffitD1 = exp(-(yp-param(1)).^2/(2*param(3).^2));
ffitD2 = exp(-(xp-param(2)).^2/(2*param(4).^2));
ffit = ffitD1.*ffitD2;

varaccount = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

%Put the means back on the original axis

p2 = param(2)*cos(-param(5)*pi/180) + param(1)*sin(-param(5)*pi/180);
p1 = param(1)*cos(-param(5)*pi/180) - param(2)*sin(-param(5)*pi/180);
param(1) = p1;
param(2) = p2;

%figure,plot(errall)