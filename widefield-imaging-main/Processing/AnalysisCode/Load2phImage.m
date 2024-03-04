function im = Load2phImage(fno,trialID)

%if trialID has 2 elements, it assumes you have entered
%(condition,repeat).  If there is 1 element, it assumes you have entered
%the (timetag).

%fno is the frame number

global datadir

if length(trialID) ==  2
    cond = trialID(1);
    rep = trialID(2);
    ttag = gettimetag(cond,rep)-1;  %They are saved as 0 to N-1
elseif length(trialID) == 1
	ttag = trialID-1;
end

ue = datadir(end-8:end-1);

fname = [datadir ue  '_' sprintf('%03d',ttag)];

var = ['f' num2str(fno)];

fname = [fname '_' var];
load(fname)

% load(fname,var)
% eval(['im = ' var ';'])

im = double(im);


