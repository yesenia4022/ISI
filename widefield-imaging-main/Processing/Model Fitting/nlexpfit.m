function yy = nlexpfit(param,xx)

B = param(1);
alp = param(2);

A = param(3);

yy = B*exp(-alp*xx) + A;


