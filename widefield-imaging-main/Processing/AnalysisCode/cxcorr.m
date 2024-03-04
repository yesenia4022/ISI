function r = cxcorr(a,b,varargin)

%% Circular correlation coefficient based on Fisher & Lee, Biometrika (1983).

if(nargin>2)  %% set the period
    per = varargin{1};
else
    per = 2*pi;
end

a = a/per*2*pi;
b = b/per*2*pi;

r = (res(a-b)^2-res(a+b)^2)/ sqrt((1-res(2*a)^2)*(1-res(2*b)^2));

end

function y = res(x)
    y = abs(mean(exp(1i*x)));
end