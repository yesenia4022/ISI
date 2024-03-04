function Tens = GetTrialData(Tlim,tcr)

%if varargin has 2 elements, it assumes you have entered
%(condition,repeat).  If there is 1 element, it assumes you have entered
%the (timetag).

%Tlim is a 2 element vector in ms

global datadir ACQinfo

if length(tcr) == 1
    ttag = tcr-1;
elseif length(tcr) == 2
    cond = tcr(1);
    rep = tcr(2);
    ttag = gettrial(cond,rep)-1;
end

ue = datadir(end-8:end-1);

fname = [datadir ue  '_' sprintf('%03d',ttag)];

if Tlim(1) == -inf & Tlim(2) == inf %get the entire trial
    Flim(1) = 1;
    Flim(2) = ACQinfo.numberOfFrames;
else
    Flim = getframeidx(Tlim,tcr);
end
    
%Flim = [50 130];
%Use this if each frame is a file
meanimage = 0;

Fdom = Flim(1):Flim(2);
Tens = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine,length(Fdom)); %pre-Allocate
for i = 1:length(Fdom)
    var = ['f' num2str(Fdom(i))];
    fnamedum = [fname '_' var];
    load(fnamedum)    
    Tens(:,:,i) = im;
end

%meanimage = meanimage/length(Flim(1):Flim(2));

