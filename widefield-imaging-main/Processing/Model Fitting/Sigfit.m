function [param ffit varaccount] = Sigfit(domain,f,varargin)

global yy xx initguess

if ~isempty(varargin)
    initguess = varargin{1};
else 
    initguess = [];
end

%%%search%%%
xx = domain;
yy = f;
param = sigfitter;
%%%%%%%%%%%

xc = param(1);
sig = param(2);
A = param(3);
B = param(4);
d = xx-xc;
ffit = A./(1 + exp(-d*sig)) + B;

%%%%
% id = find(isnan(yy.*xx));
% yy(id) = [];
% xx(id) = [];
% expect = param(1)./(1 + exp(-xx*param(2)));
varaccount = (var(yy)-var(yy-ffit))/var(yy);

%%%%
