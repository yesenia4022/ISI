function [param ffit MSE] = Expfit(domain,f,guess)

global RF x0 xx

x0 = guess;
xx = domain;

%%%search%%%
RF = f;
param = expfitter;
%%%%%%%%%%%

ffit = param(1)*exp(-param(2)*domain) + param(3);

MSE = mean((ffit-f).*(ffit-f));

