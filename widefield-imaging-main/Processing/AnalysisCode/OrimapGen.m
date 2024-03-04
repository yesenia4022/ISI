function orimap = OrimapGen(f0dum)

k = 1;
for(i=0:length(f0dum)-1)
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks
        v = pepgetvalues;
        ori(k) = v(1);
        f0{k} = f0dum{i+1};
        k = k+1;
    end
end

orimap = 0;
for i = 1:length(f0)
    orimap = orimap + f0{i}*exp(1i*2*ori(i)*pi/180);    %Linear combination
end