function xout = myBigWeiner(alpha,sp,x)

%alpha is 1/tau
%implements 1D weiner filter for an exponential decay as the blur

% no = 1;
% 
% tdom = (0:length(x)-1)*sp;
% 
% Cp = .5*(1/(no*sqrt(1+no^2)) + 1/(1+no^2));
% Cm = .5*(1/(no*sqrt(1+no^2)) - 1/(1+no^2));
% 
% Beta = alpha*sqrt(1+1/no^2);
% 
% h1 = Beta*exp(-tdom(1:10)*Beta);
% h2 = fliplr(h1);
% 
% xout = Cp*conv(x,h1) - Cm*conv(x,h2);

no =1;
tdom = (0:length(x)-1)*sp;

h = alpha*exp(-tdom*alpha);

H = fft(h)./(abs(fft(h)).^2 + no^2);

xout = real(ifft(H.*fft(x)));


%%%






