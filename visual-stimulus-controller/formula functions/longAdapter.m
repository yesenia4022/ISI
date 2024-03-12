function STout = longAdapter(ori,oriA,STnorm,STadapt)

%Could put this in the looper
%STnorm = 1; STadapt = 20; oriA = 1; STout = longAdapter(ori,oriA,STnorm,STadapt); stim_time = STout;

if ori == oriA
    STout = STadapt;
else
    STout = STnorm;
end
