function Flim = getframeidx(Tlim,varargin)

%%%Find frame limits from time limits


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

if length(varargin) ==  2  
    cond = varargin{1};
    rep = varargin{2};
elseif length(varargin) == 1
	ttag = varargin{1};
    [cond rep] = getcondrep_offL(ttag);
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



% dS = Analyzer.syncInfo{tno}.dispSyncs;
% 
% aS = Analyzer.syncInfo{tno}.acqSyncs;
% 
% sttime = dS(2) + Tlim(1)/1000;
% [dum id] = min(abs(aS-sttime));
% Flim(1) = id(1);  %closest frame to requested start time
% 
% endtime = dS(2) + Tlim(2)/1000;
% [dum id] = min(abs(aS-endtime));
% Flim(2) = id(1);  %closest frame to requested end time

