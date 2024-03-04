function Flim = getframeidx(Tlim,tcr)

%%%Find frame limits from time limits

%Tlim is relative to stimulus onset, not trial onset.

%if varargin has 2 elements, it assumes you have entered
%(condition,repeat).  If there is 1 element, it assumes you have entered
%the (timetag).

%Tlim is a 2 element vector in ms

global Analyzer

P = getParamStruct(0);

try
    %STIMfrate = Analyzer.framerate; %ms/period
    STIMfrate = Analyzer.syncInfo{1}.frameRate;
catch
    STIMfrate = 60;
    'Frame rate does not exist in Analyzer. Setting to 60Hz'
end


preScreenFlips = round(P.predelay*STIMfrate);
preDelayEnd = preScreenFlips/STIMfrate;

if length(tcr) ==  2  
    cond = tcr(1);
    rep = tcr(2);
elseif length(tcr) == 1
	ttag = tcr;
    [cond rep] = getcondrep(ttag);
end

tno = Analyzer.loops.conds{cond}.repeats{rep}.trialno;


%Flim(1) = id(1);  %closest frame to requested start time

aS = Analyzer.syncInfo{tno}.acqSyncs;
aS = aS - aS(1);

dS(1) = 0; %Display sync time: beginning of pre delay
dS(2) = preDelayEnd;

sttime = dS(2) + Tlim(1)/1000;
[dum id] = min(abs(aS-sttime));
Flim(1) = id(1);  %closest acquisition frame to requested start time

endtime = dS(2) + Tlim(2)/1000;
[dum id] = min(abs(aS-endtime));
Flim(2) = id(1);  %closest acquisition frame to requested end time


