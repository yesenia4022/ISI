function makeGratingTexture_flashCartesian

%This one builds only the images that get played for this trial.
%Importantly, this one doesn't allow for drift or contrast reversal like
%the old version in the 'archive' file.

global Mstate screenNum loopTrial %movieBlock

global Gtxtr TDim  %'playgrating' will use these

Screen('Close')  %First clean up: Get rid of all textures/offscreen windows

Gtxtr = []; TDim = [];  %reset

P = getParamStruct;

screenRes = Screen('Resolution',screenNum); %sometimes this reports 0Hz
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

%Replace it with this (comes from flipInterval)
screenRes.hz = Mstate.refresh_rate;

%The following gives inaccurate spatial frequencies
% xN = 2*Mstate.screenDist*tan(P.x_size/2*pi/180);  %grating width in cm
% xN = round(xN*pixpercmX);  %grating width in pixels
% yN = 2*Mstate.screenDist*tan(P.y_size/2*pi/180);  %grating height in cm
% yN = round(yN*pixpercmY);  %grating height in pixels

%The following assumes the screen is curved
xcm = 2*pi*Mstate.screenDist*P.x_size/360;  %stimulus width in cm
xN = round(xcm*pixpercmX);  %stimulus width in pixels
ycm = 2*pi*Mstate.screenDist*P.y_size/360;   %stimulus height in cm
yN = round(ycm*pixpercmY);  %stimulus height in pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Create Sequence for this trial%%%%%

%Make spatial phase domain
phasedom = linspace(0,360,P.n_phase+1);
phasedom = phasedom(1:end-1);
%Make orientation domain
orirange = 180;
oridom = linspace(P.ori,P.ori+orirange,P.n_ori+1);
oridom = oridom(1:end-1);
%Make spatial frequency domain
if strcmp(P.sf_domain,'log')
    sfdom = logspace(log10(P.min_sf),log10(P.max_sf),P.n_sfreq);
elseif strcmp(P.sf_domain,'lin')
    sfdom = linspace(P.min_sf,P.max_sf,P.n_sfreq);
end
sfdom = unique(sfdom);

colordom = getColorDomain(P.colorspace);

N_Im = round(P.stim_time*screenRes.hz/P.h_per); %number of images to present

%s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed);


%%
[phasemesh orimesh sfmesh colormesh] = ndgrid(1:length(phasedom),1:length(oridom),1:length(sfdom),1:length(colordom));

seqdum = [];
for l = 1:ceil(N_Im/length(sfmesh(:)))
    s = RandStream('mt19937ar','Seed',P.rseed*100 + (l-1)*100);
    %Used this one until 8/10/18
    %s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed+(l-1)*100);
    [dumval dumid] = sort(rand(s,1,length(sfmesh(:))));
    seqdum = [seqdum dumid];
end

phaseseq = phasemesh(seqdum(1:N_Im));  %The sequence index, w/o replacement
oriseq = orimesh(seqdum(1:N_Im));  %The sequence index, w/o replacement
sfseq = sfmesh(seqdum(1:N_Im));  
colorseq = colormesh(seqdum(1:N_Im)); 

%%

%How I used to do it, pre 3/23/16
% phaseseq = randi(s,[1 length(phasedom)],1,N_Im); %N_Im random indices for the "mixed bag"
% oriseq = randi(s,[1 length(oridom)],1,N_Im); %N_Im random indices for the "mixed bag"
% sfseq = randi(s,[1 length(sfdom)],1,N_Im); %N_Im random indices for the "mixed bag"
% colorseq = randi(s,[1 length(colordom)],1,N_Im); %N_Im random indices for the "mixed bag"

% blankflag = zeros(1,N_Im);
% if P.blankProb > 0
%     nblanks = round(P.blankProb*N_Im);
%     bidx = zeros(1,N_Im);    
%     bidx(1:nblanks) = 1;
%     dum = randn(1,N_Im);
%     [dum id] = sort(dum);
%     bidx = find(bidx(id));  %randomly shuffle
%     
%     %blank condition is identified with the following indices
%     phaseseq(bidx) = 1;
%     oriseq(bidx) = 1;
%     sfseq(bidx) = length(sfdom) + 1;
%     colorseq(bidx) = 1;
%     
%     blankflag(bidx) = 1;
% end


%% Replace the last values of the sequence with blanks and then reshuffle.
%The last gratings of the sequence are the ones most likely to have multiple
%presentations.
blankflag = zeros(1,N_Im);
nblanks = round(P.blankProb*N_Im);

sfseq((end-nblanks+1):end) = length(sfdom) + 1;
phaseseq((end-nblanks+1):end)  = 1;
oriseq((end-nblanks+1):end)  = 1;
colorseq((end-nblanks+1):end)  = 1;

s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed);
[dumval dumid] = sort(rand(s,1,length(sfseq)));
sfseq = sfseq(dumid);    
oriseq = oriseq(dumid); 
phaseseq = phaseseq(dumid); 
colorseq = colorseq(dumid); 

blankflag(find(sfseq == length(sfdom)+1)) = 1;

% figure,hist(phaseseq,length(phasedom))
% figure,hist(oriseq,length(oridom))
% figure,hist(sfseq,length(sfdom))
% figure,hist(colorseq,length(colordom))

%%




for i = 1:N_Im    %loop through each image in the sequence
    if ~blankflag(i)
        %Compression: Change zoom based on spatial frequency in x/y dimensions
        %independently.  Spatial frequency along x/y dimensions are modulated
        %by the grating's orientation.
        Xcycperimage = sfdom(sfseq(i))*P.x_size * abs(cos(oridom(oriseq(i))*pi/180));
        Ycycperimage = sfdom(sfseq(i))*P.y_size * abs(sin(oridom(oriseq(i))*pi/180));
        xNp = min([25*Xcycperimage+1 xN]); %No point in having the resolution higher than the screen
        yNp = min([25*Ycycperimage+1 yN]);
        
        x_ecc = P.x_size/2;
        y_ecc = P.y_size/2;
        
        x_ecc = single(linspace(-x_ecc,x_ecc,xNp));  %deg
        y_ecc = single(linspace(-y_ecc,y_ecc,yNp));
        
        [x_ecc y_ecc] = meshgrid(x_ecc,y_ecc);  %deg
        
        Im = buildImage(oridom(oriseq(i)),sfdom(sfseq(i)),phasedom(phaseseq(i)),x_ecc,y_ecc,P); %Make the shape
        putinTexture(Im,colordom,colorseq(i),P,i); %Put in texture as RGB
    else
        putinTexture(0,colordom,colorseq(i),P,i); %Blank
    end
end

TDim = size(Im);

%Save it if 'running' experiment
if Mstate.running
    Pseq = struct;
    Pseq.phaseseq = phaseseq;
    Pseq.oriseq = oriseq;
    Pseq.sfseq = sfseq;
    Pseq.colorseq = colorseq;
    
    domains = struct;
    domains.oridom = oridom;
    domains.sfdom = sfdom;
    domains.phasedom = phasedom;
    domains.colordom = colordom;
    
    saveLog_Hart(domains,Pseq)
    
    %This version of 'saveLog' comes from the old code, where I would call
    %it from 2 places (the make/play files).  It no longer needs the
    %varargin input structure, but I just kept it out of consistency (and laziness).
%     if loopTrial == 1
%         saveLog(domains)
%     else
%         saveLog(Pseq,P.rseed)  %append log file with the latest sequence
%     end
    
end


function temp = buildImage(ori,sfreq,phase,x_ecc,y_ecc,P)


sdom = x_ecc*cos(ori*pi/180) - y_ecc*sin(ori*pi/180);    %deg
sdom = sdom*sfreq*2*pi + pi; %radians
temp = cos(sdom - phase*pi/180);  

switch P.s_profile
    
    case 'sin'
        temp = temp*P.contrast/100;
        
    case 'square'
        thresh = cos(P.s_duty*pi);
        temp = sign(temp-thresh);
        temp = temp*P.contrast/100;
        
    case 'pulse'
        thresh = cos(P.s_duty*pi);
        temp = (sign(temp-thresh) + 1)/2;
        temp = temp*P.contrast/100;
        
end


function putinTexture(Im,colordom,colorID,P,i)

global Gtxtr screenPTR

%%%%%%%%%%%%%%%%%%%%%%%
%This is a total hack%%
if strcmp(P.colorspace,'DKL')
    switch colordom(colorID)
        case 4 %S
            Im = Im*.16/.86 * 2;
        case 5 %L-M
            Im = Im;
        case 6 %L+M
            Im = Im*.16/1.05;
    end
elseif strcmp(P.colorspace,'LMS')
    switch colordom(colorID)
        case 2 %L
            Im = Im;
        case 3 %M
            Im = Im*.21/.24;
        case 4 %S
            Im = Im*.21/.86 * 2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imcell{1} = Im;

Idraw = ImtoRGB(Imcell,colordom(colorID),P,[]);
Gtxtr(i) = Screen(screenPTR, 'MakeTexture', Idraw);
