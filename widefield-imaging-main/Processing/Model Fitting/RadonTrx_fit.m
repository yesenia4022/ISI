function [param ffit varacc] = RadonTrx_fit(f,oridom,pdom)

%orientation domain on y axis
%position domain on x axis

global odomglob pdomglob RF Gsfit

[Gparam Gsfit] = Gaussfit(1:length(f(:,1)),mean(f,2)',1);

Gsfit = [Gsfit Gsfit];
f = [f; fliplr(f)];
oridom = [oridom oridom+180];

odomglob = oridom;
pdomglob = pdom;

%[Pref orientation; Tuning curve sig; Sinwave amp; Sinwave phase; StdDev
%over position; Gain of funtion, baseline]

f = f-min(f(:));
f = f/sum(f(:));

%%%search%%%
RF = f;
param = RadonTrxfitter;
%%%%%%%%%%%

ffitDori = Gsfit'*ones(1,length(pdomglob));

[xpg thetag] = meshgrid(pdomglob,odomglob);

prefpos = param(1)*cos(thetag*pi/180-param(2));

ffitDpos = exp( -(xpg-prefpos).^2/(2*param(3).^2) );

ffit = param(4)*ffitDpos.*ffitDori + param(5);

varacc = (var(f(:))-var(f(:)-ffit(:)))/var(f(:));

