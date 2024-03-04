function meanimage = getTrialMean(Tlim,varargin)

%if varargin has 2 elements, it assumes you have entered
%(condition,repeat).  If there is 1 element, it assumes you have entered
%the (timetag).

%Tlim is a 2 element vector in ms

global datadir Analyzer trialno

global F1imagemu

if length(varargin) ==  2  
    cond = varargin{1};
    rep = varargin{2};
    ttag = gettimetag(cond,rep)-1; %Saved as 0 to N-1
elseif length(varargin) == 1
	ttag = varargin{1};
    [cond rep] = getcondrep_offL(ttag);
end

ue = datadir(end-8:end-1);

fname = [datadir ue  '_' sprintf('%03d',ttag)];


Flim = getframeidx(Tlim,cond,rep);  %Find frame limits from time limits
%Flim = [50 130];
%Use this if each frame is a file

Nobj = findobj('Tag','negativeSignalFlag');
negsignal = get(Nobj,'value');

Fobj = findobj('Tag','F1flag');
F1flag = get(Fobj,'value');


if ~F1flag %Compute F0 (mean)
    
    meanimage = 0;
    
    for i = Flim(1):Flim(2)
        var = ['f' num2str(i)];
        fnamedum = [fname '_' var];
        load(fnamedum)
        if negsignal
            meanimage = meanimage + 4096-(double(im));
        else
            meanimage = meanimage + double(im);
        end
        
    end
    meanimage = meanimage/length(Flim(1):Flim(2));
    
else %compute |F1|
    
    tdom = linspace(Tlim(1),Tlim(2),length(Flim(1):Flim(2)));
    try
        STIMfrate = Analyzer.framerate; %ms/period
    catch
        STIMfrate = 60;
        'Frame rate does not exist in Analyzer. Setting to 60Hz'
    end
    T = getParamVal('t_period',0)/STIMfrate*1000; %ms/cycle
    ce = exp(1i*2*pi*tdom/T);
    
    var = ['f' num2str(10)];
    fnamedum = [fname '_' var];
    load(fnamedum);
    baseIm = double(im);
    
    k = 1;
    F1image = 0;
    for i = Flim(1):Flim(2)
        var = ['f' num2str(i)];
        fnamedum = [fname '_' var];
        load(fnamedum);
        im = double(im)-baseIm; %Theoretically, not necessary.  But helps with "leakage"
       %im = double(im);

        if negsignal
            %F1image = F1image + (double(im)-meanimage)*ce(k);
            F1image = F1image + (4096-im)*ce(k);
        else
            %F1image = F1image + (double(im)-meanimage)*ce(k);
            F1image = F1image + im*ce(k);
        end
        k = k+1;
    end
    
    
    
    F1image = abs(F1image)/length(Flim(1):Flim(2)); %peak amplitude
    meanimage = F1image;
    
    %[c r] = getcondrep(trialno);  %get cond and rep for this trialno
    if rep == 1
        F1imagemu{cond} = zeros(size(F1image));
    else
        F1imagemu{cond} = F1imagemu{cond}+F1image/getnoreps(cond);
    end
end




%
%Use this if each frame is saved as a variable within the .mat file
% im = 0;
% for i = Flim(1):Flim(2)
%     var = ['f' num2str(i)];
%     load(fname,var)
%     eval(['im = im + double(' var ');']);
%     eval(['clear ' var]) %might run out of memory if you don't clear
% end
 
