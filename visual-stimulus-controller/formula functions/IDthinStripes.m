function dum = IDthinStripes(colormod,sfvec,tpvec)

if colormod == 1
    
    dum(1) = sfvec(1);
    dum(2) = tpvec(1);
    
elseif colormod == 5
    
    dum(1) = sfvec(2);
    dum(2) = tpvec(2);
    
end