function [frameid lineid synctimes] = getsyncFrameTime(Tensdum)

%This gives the acquisition frame and line for each stimulus update.  It
%also gives the time that is happens relative to the first sync.
%%
global ACQinfo

msPerFrame = ACQinfo.linesPerFrame*ACQinfo.msPerLine;

Tens = squeeze(nanmean(Tensdum,2));
dim = size(Tens);
Tensvec = medfilt1(Tens(:),3);

W = round(msPerFrame/9.2)*5;
hh = zeros(size(Tensvec));
hh(1:W) = 1/W;
hh = abs(fft(hh));
Tensvec = ifft(hh.*fft(Tensvec));

mi = prctile(Tensvec,5); Tensvec = Tensvec-mi;
ma = prctile(Tensvec,95); Tensvec = Tensvec/ma;

Tensvec = (1+sign(Tensvec-.5))/4;
Tensvec = [0; abs(diff(Tensvec))];
Tensvec = reshape(Tensvec,dim(1),dim(2));
TensvecIm = mean(Tensvec);
frameid = find(TensvecIm);

id = find((frameid)*msPerFrame < getparam('predelay')*1000);
firstsync = frameid(id)*msPerFrame;
frameid(id) = [];
if ~isempty(id)  %if we "see" the first sync
    id = find((frameid-1)*msPerFrame > (getparam('predelay')+getparam('stim_time'))*1000 + firstsync(end));
else
    id = find((frameid-1)*msPerFrame > (getparam('predelay')+getparam('stim_time'))*1000);
end
frameid(id) = [];

clear lineid
for i = 1:length(frameid)
 
    im = squeeze(Tensvec(:,frameid(i)));
    idx = find(im);
    lineid(i) = idx(1);
    
end

synctimes = (frameid-1)*msPerFrame + lineid*ACQinfo.msPerLine;
synctimes = synctimes/1000;

%%%
hper = getparam('h_per');
msupdate = hper/80;  %just an approximation

dsync = diff(synctimes);
if ~isempty(find(dsync < msupdate/2 | dsync > msupdate*1.5))
    'Messed up calculation of synctimes!!!!!!!!'
end
    
    
