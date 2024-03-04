function N = getnotrials

N = 0;
nc = getnoconditions;
for c = 1:nc
    
    nr = getnorepeats(c);
    
    N = N + nr;

end