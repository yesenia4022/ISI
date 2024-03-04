function [param ffit varacc ffitI domI pk BW] = DoGfit(f,domain)

%Ian Nauhaus

global RF logdom

%%%%%%%%%%%
mi = min(f);
f = f-mi;
E = norm(f);
f = f/E;
%%%%%%%%%%%%

%%%search%%%
logdom = domain;
RF = f;
param = DoGfitter;
%%%%%%%%%%%

ffit1 = exp(-logdom.^2/(2*param(1).^2))*param(2);
ffit2 = exp(-logdom.^2/(2*param(3).^2))*param(4);
ffit = ffit1 - ffit2 + param(5);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

ffit = ffit*E + mi;
param(2) = param(2)*E;
param(4) = param(4)*E;
param(5) = param(5)*E+mi;

%%%
%%%
N = 50;

domI = logspace(log10(logdom(1)),log10(logdom(end)),N);

ffit1 = exp(-domI.^2/(2*param(1).^2))*param(2);
ffit2 = exp(-domI.^2/(2*param(3).^2))*param(4);
ffitI = ffit1 - ffit2 + param(5);

[dum idpk] = max(ffitI);
pk = domI(idpk);
if idpk>1 %if bandpass
    tc1 = ffitI(1:idpk);
    tc2 = ffitI(idpk:end);
    [dum sflco] = min(abs(tc1-pk/2));
    sflco = domI(sflco);
    [dum sfhco] = min(abs(tc2-pk/2));
    sfhco = sfhco + idpk - 1;
    sfhco = domI(sfhco);
    BW = log2(sfhco/sflco); 
else %if lowpass
    %pk = 0;
    BW = inf;
end