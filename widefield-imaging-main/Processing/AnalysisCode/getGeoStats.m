function [mu sig] = getGeoStats(Y)

mu = geomean(Y);
sig = exp(sqrt(sum((log(Y)-log(mu)).^2)/length(Y)));