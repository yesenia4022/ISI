function [param ffit varacc ffitI domI] = SoGfit(domdum,f)

global RF dom

%%%%%%%%%%%
% mi = min(f);
% f = f-mi;
% E = max(f);
% f = f/E;
%%%%%%%%%%%%

dom = [-fliplr(domdum) domdum];
f = [fliplr(f) f];

%%%search%%%
RF = f;
param = SoGfitter;

%%%%%%%%%%%

ffit = exp(-(dom-param(1)).^2/(2*param(2).^2));
ffit = param(3)*ffit + param(4);
ffit = ffit+fliplr(ffit);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

% ffit = ffit*E + mi;
% param(3) = param(3)*E;
% param(4) = param(4)*E+mi;

%%%Now make the interpolated version%%%%%%%%%%% 
domI = linspace(dom(1),dom(end),100);
ffitI = exp(-(domI-param(1)).^2/(2*param(2).^2));
ffitI = param(3)*ffitI + param(4);
ffitI = ffitI+fliplr(ffitI);

%ffitI = ffitI*E + mi;

%%%%

dom = dom((length(dom)/2+1):end);
ffit = ffit((length(ffit)/2+1):end);

[dum id] = min(abs(dom(1)-domI));
domI = domI(id:end);
ffitI = ffitI(id:end);
