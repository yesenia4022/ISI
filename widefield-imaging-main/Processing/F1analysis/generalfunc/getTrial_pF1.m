function Tens = getTrial(Tlim,varargin)

%if varargin has 2 elements, it assumes you have entered
%(condition,repeat).  If there is 1 element, it assumes you have entered
%the (timetag).

%Tlim is a 2 element vector in ms

global datadir

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
meanimage = 0;

for i = Flim(1):Flim(2)
    var = ['f' num2str(i)];
    fnamedum = [fname '_' var];
    load(fnamedum)
    if i == 1
        Tens = zeros(size(im,1),size(im,2),length(Flim(1):Flim(2))); %preA
    end
    
    Tens(:,:,i) = im;
   % meanimage = meanimage + (4096-double(im));
end

%meanimage = meanimage/length(Flim(1):Flim(2));

%
%Use this if each frame is saved as a variable within the .mat file
% im = 0;
% for i = Flim(1):Flim(2)
%     var = ['f' num2str(i)];
%     load(fname,var)
%     eval(['im = im + double(' var ');']);
%     eval(['clear ' var]) %might run out of memory if you don't clear
% end
 
