function err = LinetotalLSfitter_handle(param)

global yy xx errvec;

A = param(1);
B = param(2);

% lindom = linspace(0,max(xx),200);
% 
% ffit1 = A*lindom + B;

% m1 = A;
% b1 = B;
% m2 = -1/m1;
% b2 = yy - m2*xx;
% xintersect = (b2-b1)/(m1-m2);
% yintersect = m1*xintersect + b1;
% 
% E = sqrt((xx-xintersect).^2 + (yy-yintersect).^2);
yhat = A*xx+B;
xhat = xx;

dx = xhat-xx;
dy = (yhat-yy);
E = abs(dx.*sin(atan(dy./dx)));

%E = abs(yy - yhat) + abs(xx-xhat);



% temp = ones(1,length(ffit1));
% E = zeros(1,length(yy));
% for i = 1:length(yy)
%     Ds = sqrt((yy(i)*temp-ffit1).^2 + (xx(i)*temp-lindom).^2);
%     E(i) = min(Ds);
% end
id = find(~isnan(E));
err = trimmean(E(id),0);

%errvec = [errvec err];

%err = nanmean(abs(ffit1-yy));



