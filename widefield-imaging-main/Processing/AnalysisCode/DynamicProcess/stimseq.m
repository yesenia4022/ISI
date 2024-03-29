function [oriseq sfreqseq phaseseq] = stimseq(conds);

for k = 1:length(conds)
%for k = 1:5

    cond = conds(k);
    %pepsetcondition(cond);
    
    %%%%%%%%%%%%%%%%%%%
    stimpars = pepgetgratingparam';  %This has to be done before pepgetsequence for unknown reasons

    %v = pepGetValues;
    %symidx = pepsearchsymbol('rseed');
    %seq = pepGetSequence(v(symidx));
    seq = pepGetSequence(k);

    Ngrats = length(stimpars(:,1));  %no of grating types
    %Nind = length(unique(seq))   %no. of grating indices
    Nind = max(seq)+1;  %Sometimes they aren't all presented
    Nblanks = Nind - Ngrats;  %no. of indices for blanks
    stimpars = [stimpars; 999*ones(Nblanks,length(stimpars(1,:)))]; %flag blanks as 999
    
    phasevals = stimpars(:,3);
    sfreqvals = stimpars(:,2);
    orivals = stimpars(:,1);
    
    ida = find(orivals ~= 999);  
    orivals(ida) = angle(exp(1i*2*orivals(ida)*pi/180))*180/pi/2;  %makes values -90 to 90
    orivals(ida) = orivals(ida) + (1-sign(orivals(ida)-eps))*90;  %make 0 to 180
    
    sfreqvals(ida) = 1./sfreqvals(ida);
    

    oriseq{k} = orivals(seq+1);
    phaseseq{k} = phasevals(seq+1);
    sfreqseq{k} = sfreqvals(seq+1);

end