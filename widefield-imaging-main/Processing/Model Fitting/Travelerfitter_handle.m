function err = Travelerfitter_handle(param)

global RF;

dim = size(RF);


[xdom tdom] = meshgrid(1:dim(2),1:dim(1));

Xweight = exp(-(xdom-param(1)).^2/(2*param(2).^2));
%Xweight = exp(-abs(xdom-param(1))/param(2));

Delay = param(3)*abs(xdom-param(1)) + param(4);
Tfunc = exp(-(tdom-Delay).^2/(2*param(5).^2));

ffit = Xweight.*Tfunc;

%ffit = param(6)*ffit + param(7);

%err = sum((ffit(:)-RF(:)).^2);

err = -sum(zscore(ffit(:)).*zscore(RF(:)))/(length(RF)-1);
