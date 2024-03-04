function F1image = CondF1(Tens,Tlim)

global Analyzer

Flim = getframeidx(Tlim,1);

tdom = linspace(Tlim(1),Tlim(2),length(Flim(1):Flim(2)));
try
    STIMfrate = Analyzer.framerate; %ms/period
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end
T = getParamVal('t_period',0)/STIMfrate*1000; %ms/cycle
ce = exp(1i*2*pi*tdom/T);



for c = 1:length(Tens)
    muImage = mean(Tens{c},3);
    
    k = 1;
    F1image{c} = 0;
    for i = Flim(1):Flim(2)
        F1image{c} = F1image{c} + (Tens{c}(:,:,i)-muImage)*ce(k);
        k = k+1;
    end
    F1image{c} = abs(F1image{c})/length(Flim(1):Flim(2)); %peak amplitude

end
