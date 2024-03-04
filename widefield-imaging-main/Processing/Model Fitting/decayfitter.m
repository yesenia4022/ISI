function [x,f,g] = decayfitter

global RF xx initguess


x0 = initguess;

    
dim = length(RF);

% options = optimset('MaxFunEvals',6000,'MaxIter',6000,'TolFun',.00004,'TolX',.00004);
% [x,f] = fminsearch('gaussfitter_handle',x0,options);

[x,f] = fminsearch('decayfitter_handle',x0);

% lb = x0-[1.5 1 .1 .1];
% ub = x0+[1.5 1 .2 .1];
% [x,f] = fmincon('gaussfitter_handle',x0,[],[],[],[],lb,ub);

A = x(1);
B = x(2);

g = 1./(A./xx +  B);
